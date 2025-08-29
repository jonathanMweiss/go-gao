// improved from recursive variant using gpt.
package field

import "errors"

// NTTForward converts a coefficient vector to NTT (evaluation-at-powers) form.
// Length must be a power of two and the field must provide an N-th primitive root.
// Transforms in-place and returns the same pointer for convenience.
func (pr *DensePolyRing) NttForward(a *Polynomial) (*Polynomial, error) {
	if a == nil || len(a.inner) == 0 {
		return a, nil
	}
	if a.isNTT {
		return a, nil // already in NTT domain
	}
	n := len(a.inner)
	if !IsPowerOfTwo(uint64(n)) {
		return nil, errors.New("NTTForward: length must be a power of two")
	}

	f := pr.GetField()
	psi, err := f.GetRootOfUnity(uint64(n))
	if err != nil {
		return nil, err
	}

	// Bit-reversal permutation (in-place)
	bitReverseInPlace(a.inner)

	// Iterative breadth-first butterflies (DIT)
	for m := 2; m <= n; m <<= 1 {
		// wm = psi^(N/m)  (principal m-th root)
		wm := f.Pow(psi, uint64(n/m))
		for k := 0; k < n; k += m {
			w := uint64(1)
			for j := 0; j < m/2; j++ {
				u := a.inner[k+j]
				t := f.Mul(w, a.inner[k+j+m/2])
				a.inner[k+j] = f.Add(u, t)
				a.inner[k+j+m/2] = f.Sub(u, t)
				w = f.Mul(w, wm)
			}
		}
	}

	a.isNTT = true
	return a, nil
}

// NTTBackward converts an NTT (evaluation) vector back to coefficient form.
// Uses the inverse root and multiplies by n^{-1} at the end.
// Transforms in-place and returns the same pointer for convenience.
func (pr *DensePolyRing) NttBackward(a *Polynomial) (*Polynomial, error) {
	if a == nil || len(a.inner) == 0 {
		return a, nil
	}
	if !a.isNTT {
		return a, nil // already in coefficient domain
	}
	n := len(a.inner)
	if !IsPowerOfTwo(uint64(n)) {
		return nil, errors.New("NTTBackward: length must be a power of two")
	}

	f := pr.GetField()
	psi, err := f.GetRootOfUnity(uint64(n))
	if err != nil {
		return nil, err
	}
	psiInv := f.Inverse(psi)

	// Bit-reversal permutation (in-place)
	bitReverseInPlace(a.inner)

	// Iterative butterflies with inverse roots
	for m := 2; m <= n; m <<= 1 {
		// wm = (psi^{-1})^(N/m)
		wm := f.Pow(psiInv, uint64(n/m))
		for k := 0; k < n; k += m {
			w := uint64(1)
			for j := 0; j < m/2; j++ {
				u := a.inner[k+j]
				t := f.Mul(w, a.inner[k+j+m/2])
				a.inner[k+j] = f.Add(u, t)
				a.inner[k+j+m/2] = f.Sub(u, t)
				w = f.Mul(w, wm)
			}
		}
	}

	// Multiply every entry by n^{-1} to finish the inverse
	nInv := f.Inverse(uint64(n))
	for i := 0; i < n; i++ {
		a.inner[i] = f.Mul(a.inner[i], nInv)
	}

	a.isNTT = false

	pr.trimTrailingZeros(a)

	return a, nil
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
