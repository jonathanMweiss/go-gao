// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"errors"
	"strconv"
	"strings"
)

// Polynomial is a polynomial over a fixed field, held either as coefficients or, after
// an NTT, as evaluations. Build one with PolyRing.NewPolynomial.
type Polynomial struct {
	f     Field
	inner []uint64
	isNTT bool
}

var (
	errModulusMismatch = errors.New("polynomials must be over the same field")
	errNTTMismatch     = errors.New("polynomials must be both in NTT or both in coefficient representation")
	errDegreeMismatch  = errors.New("polynomials must be of the same degree")
)

func preOpVerification(p, q *Polynomial) error {
	if p.f.Modulus() != q.f.Modulus() {
		return errModulusMismatch
	}

	if p.isNTT != q.isNTT {
		return errNTTMismatch
	}

	if p.isNTT && len(p.inner) != len(q.inner) {
		return errDegreeMismatch
	}

	return nil
}

// IsZero reports whether p is the zero polynomial: every coefficient zero.
func (p *Polynomial) IsZero() bool {
	// leadingCoeffPos is negative exactly when no coefficient is non-zero, which also
	// covers an empty inner.
	return p.leadingCoeffPos() < 0
}

// Equals reports whether p and q are the same polynomial, regardless of how many
// high-order zero coefficients either one is padded with. It returns false for
// polynomials over different fields or in different domains.
func (p *Polynomial) Equals(q *Polynomial) bool {
	if err := preOpVerification(p, q); err != nil {
		return false
	}

	deg := p.Degree()
	if deg != q.Degree() {
		return false
	}

	fld := p.f
	for i := 0; i <= deg; i++ {
		if !fld.Equals(p.inner[i], q.inner[i]) {
			return false
		}
	}

	return true
}

// Degree returns the degree of p, or -1 if p is the zero polynomial.
//
// Callers must treat a negative result as "no degree" before doing arithmetic on it;
// the convention throughout this package is to test Degree() < 0.
func (p *Polynomial) Degree() int {
	return p.leadingCoeffPos()
}

// LeadCoeff returns the highest-degree non-zero coefficient, or 0 if p is zero.
func (p *Polynomial) LeadCoeff() uint64 {
	if pos := p.leadingCoeffPos(); pos >= 0 {
		return p.inner[pos]
	}

	return 0
}

// leadingCoeffPos returns the index of the highest non-zero coefficient, or -1 when
// every coefficient is zero. -1 rather than a large negative sentinel so that callers
// computing degree differences cannot overflow.
func (p *Polynomial) leadingCoeffPos() int {
	for i := len(p.inner) - 1; i >= 0; i-- {
		if p.inner[i] != 0 {
			return i
		}
	}

	return -1
}

func (p *Polynomial) removeLeadingZeroes() {
	if p.isNTT {
		return
	}

	lead := p.leadingCoeffPos()
	if lead < 0 {
		if cap(p.inner) == 0 {
			p.inner = []uint64{0} // must allocate.

			return
		}

		// reusing the backing array
		p.inner = p.inner[:1]
		p.inner[0] = 0

		return
	}

	p.inner = p.inner[:lead+1]
}

// Copy returns a deep copy of p, sharing no memory with it.
func (p *Polynomial) Copy() *Polynomial {
	innercopy := make([]uint64, len(p.inner))
	copy(innercopy, p.inner)

	return &Polynomial{f: p.f, inner: innercopy, isNTT: p.isNTT}
}

// String renders p in descending degree order, as in "5*x^2 + 3*x^1 + 7". Zero
// coefficients are omitted, and the zero polynomial renders as "0".
//
// It is a debugging and test aid, not a parseable format: p is left untouched, at the
// cost of copying it.
func (p *Polynomial) String() string {
	q := p.Copy()
	q.removeLeadingZeroes()

	if len(q.inner) == 1 {
		return strconv.FormatUint(q.inner[0], 10)
	}

	bldr := strings.Builder{}

	for i := len(q.inner) - 1; i >= 0; i-- {
		if q.inner[i] == 0 {
			continue
		}

		if bldr.Len() > 0 {
			bldr.WriteString(" + ")
		}

		bldr.WriteString(strconv.FormatUint(q.inner[i], 10))

		if i != 0 {
			bldr.WriteString("*x^")
			bldr.WriteString(strconv.FormatInt(int64(i), 10))
		}
	}

	if bldr.Len() == 0 {
		return "0"
	}

	return bldr.String()
}

// ToSlice returns a copy of the coefficients, lowest degree first.
func (p *Polynomial) ToSlice() []uint64 {
	list := make([]uint64, len(p.inner))
	copy(list, p.inner)

	return list
}

// NoCopySlice returns the backing coefficient array without copying.
//
// Mutating it mutates p, and in the NTT domain resizing it breaks the transform's
// length invariant. Prefer ToSlice unless the copy is genuinely too costly.
func (p *Polynomial) NoCopySlice() []uint64 {
	return p.inner
}

// IsCoeffMode reports whether p holds coefficients rather than evaluations.
func (p *Polynomial) IsCoeffMode() bool {
	return !p.isNTT
}
