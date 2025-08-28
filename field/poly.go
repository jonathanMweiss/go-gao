package field

import (
	"math"
	"strconv"
	"strings"
)

type Polynomial struct {
	f     Field
	inner []uint64
	isNTT bool
}

/*
Polynomial expects the coefficients to be in the same field
and ordered from lowest to highest degree. (e.g. [1, 2, 3] is 1 + 2x + 3x^2)

Can be point representation, generated from numerous evaluation points.
*/
func NewPolynomial(f Field, inner []uint64, isPointRepresentation bool) *Polynomial {
	// validate inner are all in the same field
	// validate inner
	if len(inner) == 0 {
		panic("empty polynomial")
	}

	return &Polynomial{
		inner: inner,
		isNTT: isPointRepresentation,
		f:     f,
	}
}

func preOpVerification(p, q *Polynomial) bool {
	if p.f.Modulus() != q.f.Modulus() {
		return false
	}

	if p.isNTT != q.isNTT {
		return false
	}

	if p.isNTT {
		return len(p.inner) == len(q.inner)
	}

	return true
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

func (p *Polynomial) Equals(q *Polynomial) bool {
	if !preOpVerification(p, q) {
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

func (p *Polynomial) Degree() int {
	return p.leadingCoeffPos()
}

func (p *Polynomial) LeadCoeff() uint64 {
	if pos := p.leadingCoeffPos(); pos >= 0 {
		return p.inner[pos]
	}

	return 0
}

func (p *Polynomial) leadingCoeffPos() int {
	for i := len(p.inner) - 1; i >= 0; i-- {
		if p.inner[i] != 0 {
			return i
		}
	}

	return math.MinInt
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

	return NewPolynomial(p.f, innercopy, p.isNTT)
}

// todo: fix
func (p *Polynomial) String() string {
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

func (p *Polynomial) IsCoeffMode() bool {
	return p.isNTT
}
