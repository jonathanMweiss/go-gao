// improved from recursive variant(+cache of twiddles) using gpt.
package field

import "errors"

type twiddleSet struct {
	// For each stage s (m = 2<<s), fwd[s] (and inv[s]) has length m/2
	// holding w^j where w = psi^(n/m) for forward, and w = psiInv^(n/m) for inverse.
	fwd  [][]uint64
	inv  [][]uint64
	nInv uint64 // inverse of n (for inverse NTT scaling)
}

func (pr *DensePolyRing) getTwiddles(n int) (*twiddleSet, error) {
	pr.mu.RLock()
	if ts, ok := pr.twiddleCache[n]; ok {
		pr.mu.RUnlock()
		return ts, nil
	}
	pr.mu.RUnlock()

	// Build outside lock
	if n <= 1 {
		ts := &twiddleSet{
			fwd:  [][]uint64{},
			inv:  [][]uint64{},
			nInv: pr.Inverse(uint64(n)),
		}

		pr.mu.Lock()
		pr.twiddleCache[n] = ts
		pr.mu.Unlock()

		return ts, nil
	}
	psi, err := pr.GetRootOfUnity(uint64(n))
	if err != nil {
		return nil, err
	}
	psiInv := pr.Inverse(psi)

	var fwd [][]uint64
	var inv [][]uint64

	// stages: m = 2,4,8,...,n  => stage index s = 0..(log2(n)-1)
	for m := 2; m <= n; m = m << 1 {
		half := m >> 1
		wmF := pr.Pow(psi, uint64(n/m))    // forward stage root
		wmI := pr.Pow(psiInv, uint64(n/m)) // inverse stage root

		rowF := make([]uint64, half)
		rowI := make([]uint64, half)

		wF := uint64(1)
		wI := uint64(1)
		for j := 0; j < half; j++ {
			rowF[j] = wF
			rowI[j] = wI
			wF = pr.Mul(wF, wmF)
			wI = pr.Mul(wI, wmI)
		}

		fwd = append(fwd, rowF)
		inv = append(inv, rowI)
	}

	ts := &twiddleSet{
		fwd:  fwd,
		inv:  inv,
		nInv: pr.Inverse(uint64(n)),
	}

	pr.mu.Lock()
	defer pr.mu.Unlock()
	// Another goroutine may have won the race; keep the first one but return ours if we’re first.
	if existing, ok := pr.twiddleCache[n]; ok {
		return existing, nil
	}

	pr.twiddleCache[n] = ts

	return ts, nil
}
func (pr *DensePolyRing) NttForward(a *Polynomial) error {
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

	// Bit-reversal permutation (in place; allocation-free)
	bitReverseInPlace(a.inner)

	// Twiddles per stage
	ts, err := pr.getTwiddles(n)
	if err != nil {
		return err
	}

	// Stages: m = 2,4,8,...,n  with precomputed ws per stage.
	for s, m := 0, 2; m <= n; s, m = s+1, m<<1 {
		half := m >> 1
		ws := ts.fwd[s] // length = half
		for k := 0; k < n; k += m {
			// breadth-first butterflies
			for j := 0; j < half; j++ {
				u := a.inner[k+j]
				t := pr.Mul(ws[j], a.inner[k+j+half])
				a.inner[k+j] = pr.Add(u, t)
				a.inner[k+j+half] = pr.Sub(u, t)
			}
		}
	}

	a.isNTT = true

	return nil
}

func (pr *DensePolyRing) NttBackward(a *Polynomial) error {
	if err := pr.nttBackwardNoTrim(a); err != nil {
		return err
	}
	pr.trimTrailingZeros(a)

	return nil
}

func (pr *DensePolyRing) nttBackwardNoTrim(a *Polynomial) error {
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

	// Bit-reversal permutation (in place)
	bitReverseInPlace(a.inner)

	// Twiddles per stage
	ts, err := pr.getTwiddles(n)
	if err != nil {
		return err
	}

	// Inverse butterflies use inverse stage twiddles
	for s, m := 0, 2; m <= n; s, m = s+1, m<<1 {
		half := m >> 1
		ws := ts.inv[s]
		for k := 0; k < n; k += m {
			for j := 0; j < half; j++ {
				u := a.inner[k+j]
				t := pr.Mul(ws[j], a.inner[k+j+half])
				a.inner[k+j] = pr.Add(u, t)
				a.inner[k+j+half] = pr.Sub(u, t)
			}
		}
	}

	// scale by n^{-1}
	for i := 0; i < n; i++ {
		a.inner[i] = pr.Mul(a.inner[i], ts.nInv)
	}

	a.isNTT = false
	return nil
}

func bitReverseInPlace(xs []uint64) {
	n := len(xs)
	if n <= 1 {
		return
	}
	// Compute number of bits needed
	bits := 0
	for (1 << bits) < n {
		bits++
	}
	j := 0
	for i := 1; i < n-1; i++ {
		bit := n >> 1
		for j&bit != 0 {
			j &= ^bit
			bit >>= 1
		}
		j |= bit
		if i < j {
			xs[i], xs[j] = xs[j], xs[i]
		}
	}
}
