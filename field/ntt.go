package field

import (
	"fmt"
)

// n must be a power of two. Returns primitive n-th root w and its inverse.
func (r *DensePolyRing) primitiveRoot(n int) (w, winv uint64) {
	ru, err := r.GetRootOfUnity(uint64(n))
	if err != nil {
		panic(err) // simplify; surface differently if you prefer
	}
	w = ru
	winv = r.Inverse(w)
	return
}

// in-place bit-reversal permutation
func bitReverseInPlace(a []uint64) {
	n := len(a)
	j := 0
	for i := 1; i < n-1; i++ {
		bit := n >> 1
		for j&bit != 0 {
			j ^= bit
			bit >>= 1
		}
		j |= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
}

// convert p (coeff domain) -> NTT of size n; n must be pow2 >= len(p)
func (r *DensePolyRing) NTT(p *Polynomial) *Polynomial {
	if p.isNTT {
		panic("NTT called on already NTT polynomial")
	}

	if !IsPowerOfTwo(uint64(len(p.inner))) {
		panic("NTT called on non-power-of-two polynomial")
	}

	n := len(p.inner)

	out := p.Copy()
	out.isNTT = true

	// forward NTT
	w, _ := r.primitiveRoot(n)
	r.nttInPlace(out.inner, w)

	return out
}

func (r *DensePolyRing) nttInPlace(a []uint64, w uint64) {
	n := len(a)
	bitReverseInPlace(a)
	for len_ := 2; len_ <= n; len_ <<= 1 {
		// wn = w^(n/len_)
		step := n / len_
		wn := r.Pow(w, uint64(step))
		for i := 0; i < n; i += len_ {
			omega := uint64(1)
			half := len_ >> 1
			for j := 0; j < half; j++ {
				u := a[i+j]
				v := r.Mul(a[i+j+half], omega)
				a[i+j] = r.Add(u, v)
				a[i+j+half] = r.Sub(u, v)
				omega = r.Mul(omega, wn)
			}
		}
	}
}

func copyIntoReduced(f Field, src *Polynomial, c *Polynomial) {
	ensureLen(c, len(src.inner))
	for i, v := range src.inner {
		if src.isNTT { // isNTT == true
			c.inner[i] = v // already reduced in NTT
		} else {
			c.inner[i] = f.Reduce(v)
		}
	}
	c.f = f
	c.isNTT = src.isNTT
}

// convert p (NTT domain) -> coefficient domain, shrinking trailing zeros
func (r *DensePolyRing) INTT(p *Polynomial) *Polynomial {
	if !p.isNTT {
		panic("INTT called on non-NTT polynomial")
	}
	out := &Polynomial{f: r.Field, isNTT: false}
	ensureLen(out, len(p.inner))
	copy(out.inner, p.inner)
	_, winv := r.primitiveRoot(len(out.inner))
	r.inttInPlace(out.inner, winv)
	r.trimTrailingZeros(out)
	return out
}

func (r *DensePolyRing) inttInPlace(a []uint64, wInv uint64) {
	n := len(a)
	bitReverseInPlace(a)
	for len_ := 2; len_ <= n; len_ <<= 1 {
		step := n / len_
		wn := r.Pow(wInv, uint64(step))
		for i := 0; i < n; i += len_ {
			omega := uint64(1)
			half := len_ >> 1
			for j := 0; j < half; j++ {
				u := a[i+j]
				v := a[i+j+half]
				a[i+j] = r.Add(u, v)
				t := r.Sub(u, v)
				a[i+j+half] = r.Mul(t, omega)
				omega = r.Mul(omega, wn)
			}
		}
	}
	// scale by n^{-1}
	nInv := r.Inverse(uint64(n))
	for i := range a {
		a[i] = r.Mul(a[i], nInv)
	}
}

func NextPow2(n uint64) uint64 {
	if n == 0 {
		return 1
	}

	n |= n >> 1 // Divide by 2^k for consecutive doublings of k up to 32,
	n |= n >> 2 // and then or the results.
	n |= n >> 4
	n |= n >> 8
	n |= n >> 16
	n |= n >> 32

	// The result is a number of 1 bits equal to the number
	// of bits in the original number, plus 1. That's the
	// next highest power of 2.
	return n + 1
}

// PadForNTTWithFutureMults pads p.inner once so it can hold the result of a chain
// of multiplications by polynomials whose lengths are given in others.
//
// If p has length n0 and you plan to multiply by polys of lengths n1, n2, ... nk,
// the final (linear convolution) length is: L = n0 + n1 + ... + nk - k.
// We pad to the next power of two >= L.
//
// Example: starting len 1024, then multiply by 1024 and 1024:
// L = 1024 + 1024 + 1024 - 2 = 3070 → padded to 4096.
func PadForNTTWithFutureMults(p *Polynomial, others ...int) error {
	if p == nil {
		return fmt.Errorf("nil polynomial")
	}
	if p.isNTT {
		return fmt.Errorf("cannot pad: polynomial is already in NTT domain")
	}

	n0 := len(p.inner)
	if n0 == 0 {
		// Decide your policy; here we allow planning only if future mults exist.
		if len(others) == 0 {
			return nil
		}
	}

	// Compute target linear-convolution length after chaining all multiplies.
	targetLen := n0
	for _, n := range others {
		if n <= 0 {
			return fmt.Errorf("invalid operand length %d (must be > 0)", n)
		}
		targetLen += n
	}
	targetLen -= len(others) // subtract number of multiplications

	if targetLen <= 0 {
		// Shouldn't happen with valid inputs, but keep it safe.
		targetLen = 1
	}

	// Round up to the next power of two.
	padTo := NextPow2(uint64(targetLen))

	// Grow with zeros if needed.
	ensureLen(p, int(padTo))
	return nil
}
