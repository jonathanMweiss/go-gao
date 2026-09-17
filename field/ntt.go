// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"errors"
	"math/bits"
	"slices"
)

type twiddleSet struct {
	// For each stage s (m = 2<<s), fwd[s] (and inv[s]) has length m/2
	// holding w^j where w = psi^(n/m) for forward, and w = psiInv^(n/m) for inverse.
	fwd [][]uint64
	inv [][]uint64
	// fwdShoup[s][j] and invShoup[s][j] are the Shoup companions of fwd[s][j] and
	// inv[s][j], and nInvShoup is nInv's; see [PrimeField.mulShoup].
	fwdShoup  [][]uint64
	invShoup  [][]uint64
	nInv      uint64 // inverse of n (for inverse NTT scaling)
	nInvShoup uint64
}

func (r *PolyRing) getTwiddles(n int) (*twiddleSet, error) {
	r.mu.RLock()
	if ts, ok := r.twiddleCache[n]; ok {
		r.mu.RUnlock()
		return ts, nil
	}
	r.mu.RUnlock()

	// Build outside lock
	if n <= 1 {
		nInv := r.f.Inverse(uint64(n))
		ts := &twiddleSet{
			fwd:       [][]uint64{},
			inv:       [][]uint64{},
			fwdShoup:  [][]uint64{},
			invShoup:  [][]uint64{},
			nInv:      nInv,
			nInvShoup: r.f.shoupFactor(nInv),
		}

		r.mu.Lock()
		r.twiddleCache[n] = ts
		r.mu.Unlock()

		return ts, nil
	}
	psi, err := RootOfUnity(r.f, uint64(n))
	if err != nil {
		return nil, err
	}
	psiInv := r.f.Inverse(psi)

	var fwd [][]uint64
	var inv [][]uint64
	var fwdShoup [][]uint64
	var invShoup [][]uint64

	f := r.f
	// stages: m = 2,4,8,...,n  => stage index s = 0..(log2(n)-1)
	for m := 2; m <= n; m = m << 1 {
		half := m >> 1
		wmF := r.f.Pow(psi, uint64(n/m))    // forward stage root
		wmI := r.f.Pow(psiInv, uint64(n/m)) // inverse stage root

		rowF := make([]uint64, half)
		rowI := make([]uint64, half)
		rowFS := make([]uint64, half)
		rowIS := make([]uint64, half)

		wF := uint64(1)
		wI := uint64(1)
		for j := 0; j < half; j++ {
			rowF[j] = wF
			rowI[j] = wI
			rowFS[j] = f.shoupFactor(wF)
			rowIS[j] = f.shoupFactor(wI)
			wF = f.Mul(wF, wmF)
			wI = f.Mul(wI, wmI)
		}

		fwd = append(fwd, rowF)
		inv = append(inv, rowI)
		fwdShoup = append(fwdShoup, rowFS)
		invShoup = append(invShoup, rowIS)
	}

	nInv := r.f.Inverse(uint64(n))
	ts := &twiddleSet{
		fwd:       fwd,
		inv:       inv,
		fwdShoup:  fwdShoup,
		invShoup:  invShoup,
		nInv:      nInv,
		nInvShoup: f.shoupFactor(nInv),
	}

	r.mu.Lock()
	defer r.mu.Unlock()
	// Another goroutine may have won the race; keep the first one but return ours if we’re first.
	if existing, ok := r.twiddleCache[n]; ok {
		return existing, nil
	}

	r.twiddleCache[n] = ts

	return ts, nil
}

// NttForward transforms a into the NTT domain in place. If polynomial is already in NTT form, no operation is performed.
func (r *PolyRing) NttForward(a *Polynomial) error {
	if a == nil || len(a.inner) == 0 {
		return nil
	}
	if a.isNTT {
		return nil
	}
	n := len(a.inner)
	if !IsPowerOfTwo(uint64(n)) {
		return errors.New("NTTForward: length must be a power of two")
	}

	// Put the coefficients in the order nttRecursive's leaves would be in, so that the
	// stages below can work upward from there. See nttRecursive for why that order is
	// the bit-reversal of the natural one.
	bitReverseInPlace(a.inner)

	// Twiddles per stage
	ts, err := r.getTwiddles(n)
	if err != nil {
		return err
	}

	f := r.f
	inner := a.inner

	// One level of nttRecursive's tree per pass, deepest level first: m is the size of
	// the subproblem this stage finishes, so stage 0 merges single coefficients into
	// pairs and the last stage merges the two halves of the whole array. Where the
	// recursion reaches a level by returning into it, this sweeps levels bottom-up:
	// a merge needs only its own two children, never its cousins.
	for s, m := 0, 2; m <= n; s, m = s+1, m<<1 {
		half := m >> 1
		// s used to index the "level" of twiddles.
		// omegas for this stage.
		ws := ts.fwd[s][:half]
		wsp := ts.fwdShoup[s][:len(ws)]

		// k helps us choose the two halves to merge.
		// The permutation above is what leaves each subproblem contiguous, so [k, k+m) is
		// exactly one node's territory.
		for k := 0; k < n; k += m {
			// The halves are the two children, each already transformed by the
			// previous stage: lo is nttRecursive's yeven, hi its yodd.
			lo, hi := inner[k:k+half], inner[k+half:k+half+half]

			// The merge, written over the children in place -- lo[j] becomes y[j] and
			// hi[j] becomes y[j+half]. Both operands are read into u and t before
			// either write, because output and input share this memory: assigning
			// lo[j] first would clobber the value the next line still needs.
			for j, w := range ws {
				u := lo[j]
				t := f.mulShoup(w, wsp[j], hi[j])
				lo[j] = f.Add(u, t)
				hi[j] = f.Sub(u, t)
			}
		}
	}

	a.isNTT = true

	return nil
}

// NttBackward transforms a from the NTT domain back to coefficients, in place, and
// trims trailing zeros. It returns an error if a is not in the NTT domain or its
// length is not a power of two.
func (r *PolyRing) NttBackward(a *Polynomial) error {
	if err := r.nttBackwardNoTrim(a); err != nil {
		return err
	}
	a.trimTrailingZeros()

	return nil
}

func (r *PolyRing) nttBackwardNoTrim(a *Polynomial) error {
	if a == nil || len(a.inner) == 0 {
		return nil
	}
	if !a.isNTT {
		return errors.New("newMethod: polynomial is not in NTT form")
	}

	n := len(a.inner)
	if !IsPowerOfTwo(uint64(n)) {
		return errors.New("NTTBackward: length must be a power of two")
	}

	// Bit-reversal permutation
	bitReverseInPlace(a.inner)

	// Twiddles per stage
	ts, err := r.getTwiddles(n)
	if err != nil {
		return err
	}

	f := r.f
	inner := a.inner

	// The same three loops as NttForward -- level, node, merge; see there for what each
	// one walks -- run over the inverse stage roots. That computes n times the inverse
	// transform, which the scaling below divides out.
	for s, m := 0, 2; m <= n; s, m = s+1, m<<1 {
		half := m >> 1
		ws := ts.inv[s][:half]
		wsp := ts.invShoup[s][:len(ws)]
		for k := 0; k < n; k += m {
			lo, hi := inner[k:k+half], inner[k+half:k+half+half]
			for j, w := range ws {
				u := lo[j]
				t := f.mulShoup(w, wsp[j], hi[j])
				lo[j] = f.Add(u, t)
				hi[j] = f.Sub(u, t)
			}
		}
	}

	// scale by n^{-1}
	for i, v := range inner {
		inner[i] = f.mulShoup(ts.nInv, ts.nInvShoup, v)
	}

	a.isNTT = false
	return nil
}

// bitReverseInPlace rearranges the elements of xs in place according to the "bit-reversal permutation".
// i.e., similar to FFT's element ordering [0,1,2,3,4,5,6,7] -> [0,4,2,6,1,5,3,7].
// this function does so by finding for each index i its new position j by reversing the bits of i,
// and swapping xs[i] with xs[j] if i < j. It does this in O(n) time and O(1) space.
//
// For example, if xs has length 8, the indices 0..7 in binary are:
// 000, 001, 010, 011, 100, 101, 110, 111
// and their bit-reversals are:
// 000, 100, 010, 110, 001, 101, 011, 111
// so the elements of xs are rearranged to match these new indices.
//
// this is needed to run NTT from the bottom up, as the Cooley-Tukey algorithm requires.
func bitReverseInPlace(xs []uint64) {
	n := len(xs)
	if n <= 1 {
		return
	}

	shift := uint(64-bits.TrailingZeros(uint(n))) & 63
	for i := 1; i < n-1; i++ {
		j := int(bits.Reverse64(uint64(i))>>shift) & (n - 1)
		if i < j {
			xs[i], xs[j] = xs[j], xs[i]
		}
	}
}

// nttRecursive is the textbook recursive NTT, kept as the reference the iterative
// [PolyRing.NttForward] and [PolyRing.NttBackward] are derived from.
//
// len(a) must be a power of two that the field admits a transform of; anything else is
// a programming error and panics.
func nttRecursive(f Field, a []uint64) []uint64 {
	n := len(a)
	if n <= 1 {
		return append([]uint64(nil), a...)
	}

	w, err := RootOfUnity(f, uint64(n))
	if err != nil {
		panic("nttRecursive: " + err.Error())
	}

	half := n >> 1

	// split coeffs into even and odd.
	even := make([]uint64, half)
	odd := make([]uint64, half)
	for i := range half {
		even[i] = a[2*i]
		odd[i] = a[2*i+1]
	}

	yeven := nttRecursive(f, even)
	yodd := nttRecursive(f, odd)

	y := make([]uint64, n)
	for k := range half {
		// w^(k+n/2) = -w^k, so the upper half is the same product subtracted: that
		// sign is why one pass over half the indices fills all n, and it is what the
		// stage loops above write as lo[j], hi[j] = u+t, u-t.
		wk := f.Pow(w, uint64(k))
		y[k] = f.Add(yeven[k], f.Mul(wk, yodd[k]))
		y[k+half] = f.Sub(yeven[k], f.Mul(wk, yodd[k]))
	}

	return y
}

func nttRecursiveInverse(f Field, a []uint64) []uint64 {
	if len(a) <= 1 {
		return append([]uint64(nil), a...)
	}

	y := nttRecursive(f, a)
	slices.Reverse(y[1:])

	// normalize.
	nInv := f.Inverse(uint64(len(a)))
	for i, v := range y {
		y[i] = f.Mul(v, nInv)
	}

	return y
}
