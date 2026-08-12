// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"errors"
	"strconv"
	"strings"
)

// Polynomials are built through PolyRing.NewPolynomial.
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

func (p *Polynomial) IsZero() bool {
	if len(p.inner) == 0 {
		return true
	}

	if len(p.inner) == 1 && p.inner[0] == 0 {
		return true
	}

	pos := p.leadingCoeffPos()
	for i := 0; i < pos; i++ {
		if p.inner[i] != 0 {
			return false
		}
	}

	return true
}

// Polynomial must be trim from leading zeros.
func (p *Polynomial) Equals(q *Polynomial) bool {
	if err := preOpVerification(p, q); err != nil {
		return false
	}

	if len(p.inner) != len(q.inner) {
		return false
	}

	fld := p.f
	for i := range p.inner {
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
		p.inner = []uint64{0}

		return
	}

	p.inner = p.inner[:lead+1]
}

func (p *Polynomial) Copy() *Polynomial {
	innercopy := make([]uint64, len(p.inner))
	copy(innercopy, p.inner)

	return &Polynomial{f: p.f, inner: innercopy, isNTT: p.isNTT}
}

// todo: fix
// Used mainly for testing., copies the polynomial.
func (p_ *Polynomial) String() string {
	p := p_.Copy()
	p.removeLeadingZeroes()

	if len(p.inner) == 1 {
		return strconv.FormatUint(p.inner[0], 10)
	}

	bldr := strings.Builder{}

	for i := len(p.inner) - 1; i >= 0; i-- {
		if p.inner[i] == 0 {
			continue
		}

		strI := strconv.FormatInt(int64(i), 10)

		strElem := strconv.FormatUint(p.inner[i], 10)
		bldr.WriteString(strElem)

		if i != 0 {
			bldr.WriteString("*x^")
			bldr.WriteString(strI)
			bldr.WriteString(" + ")
		}
	}

	return bldr.String()
}

func (p *Polynomial) ToSlice() []uint64 {
	list := make([]uint64, len(p.inner))
	copy(list, p.inner)

	return list
}

func (p *Polynomial) NoCopySlice() []uint64 {
	return p.inner
}

func (p *Polynomial) IsCoeffMode() bool {
	return !p.isNTT
}
