package field

import (
	"math"
	"strconv"
	"strings"
)

type Polynomial struct {
	f                *PrimeField
	inner            []Elem
	isCoefficientMod bool
}

/*
Polynomial expects the coefficients to be in the same field
and ordered from lowest to highest degree. (e.g. [1, 2, 3] is 1 + 2x + 3x^2)

Can be point representation, generated from numerous evaluation points.
*/
func NewPolynomial(inner []Elem, isPointRepresentation bool) *Polynomial {
	// validate inner are all in the same field
	// validate inner
	if len(inner) == 0 {
		panic("empty polynomial")
	}

	return &Polynomial{
		inner:            inner,
		isCoefficientMod: isPointRepresentation,
		f:                inner[0].field,
	}
}

func preOpVerifcation(p, q *Polynomial) bool {
	if p.f.prime != q.f.prime {
		return false
	}

	if p.isCoefficientMod != q.isCoefficientMod {
		return false
	}

	if p.isCoefficientMod {
		return len(p.inner) == len(q.inner)
	}

	p.removeLeadingZeroes()
	q.removeLeadingZeroes()

	return true
}

func (p *Polynomial) IsZero() bool {
	return len(p.inner) == 0 || p.inner[0].Value() == 0
}

func (p *Polynomial) Add(q *Polynomial) *Polynomial {
	if !preOpVerifcation(p, q) {
		return nil // this is an error.
	}

	size := len(p.inner)
	if len(q.inner) > size {
		size = len(q.inner)
	}

	inner := make([]Elem, size)

	copy(inner, p.inner)

	for i := range q.inner {
		inner[i] = inner[i].Add(q.inner[i])
	}

	return NewPolynomial(inner, p.isCoefficientMod)
}
func (p *Polynomial) Sub(q *Polynomial) *Polynomial {
	if !preOpVerifcation(p, q) {
		return nil // this is an error.
	}

	negative := make([]Elem, len(q.inner))
	for i, v := range q.inner {
		negative[i] = v.Neg()
	}

	sub := p.Add(NewPolynomial(negative, q.isCoefficientMod))
	sub.removeLeadingZeroes()

	return sub
}

func (p *Polynomial) Mul(q *Polynomial) *Polynomial {
	if !preOpVerifcation(p, q) {
		return nil
	}

	var inner []Elem
	if p.isCoefficientMod {
		inner = make([]Elem, len(p.inner))
		for i := range p.inner {
			inner[i] = p.inner[i].Mul(q.inner[i])
		}
	} else {
		// regular polynomial multiplication O(n^2)
		inner = make([]Elem, len(p.inner)+len(q.inner)-1)
		for i := range inner {
			inner[i] = p.f.ElemFromUint64(0)
		}

		for i := range p.inner {
			for j := range q.inner {
				tmp := p.inner[i].Mul(q.inner[j])
				inner[i+j] = inner[i+j].Add(tmp)
			}
		}
	}

	prod := NewPolynomial(inner, p.isCoefficientMod)
	prod.removeLeadingZeroes()

	return prod
}

func (p *Polynomial) Eval(x uint64) Elem {
	xElem := p.f.ElemFromUint64(x)
	result := p.f.ElemFromUint64(0)

	// horner's rule:
	for i := len(p.inner) - 1; i >= 0; i-- {
		result = p.inner[i].Add(result.Mul(xElem))
	}

	return result
}
func (p *Polynomial) Equals(q *Polynomial) bool {
	if !preOpVerifcation(p, q) {
		return false
	}

	if len(p.inner) != len(q.inner) {
		return false
	}

	for i := range p.inner {
		if !p.inner[i].Equals(q.inner[i]) {
			return false
		}
	}

	return true
}

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (p *Polynomial) LongDiv(v *Polynomial) (q *Polynomial, r *Polynomial) {
	if !preOpVerifcation(p, v) {
		return nil, nil
	}

	if v.isCoefficientMod {
		return nil, nil
	}

	n, m := p.Degree(), v.Degree()

	b := v.Copy()
	tmp := v.LeadCoeff()
	u := tmp.Inverse()

	r = p.Copy()
	qInner := make([]Elem, n-m+1)

	for i := n - m; i >= 0; i-- {
		if r.Degree() == m+i {
			qInner[i] = r.LeadCoeff().Mul(u)
			r = r.Sub(monomialMultPoly(qInner[i], i, b))
		} else {
			qInner[i] = p.f.ElemFromUint64(0)
		}
	}

	r.removeLeadingZeroes()

	if len(qInner) == 0 {
		qInner = []Elem{p.f.ElemFromUint64(0)}
	}

	q = NewPolynomial(qInner, false)
	q.removeLeadingZeroes()

	return q, r
}

func monomialMultPoly(ai Elem, deg int, p *Polynomial) *Polynomial {
	newDegree := len(p.inner) + deg

	prod := make([]Elem, newDegree)
	for i := range p.inner {
		prod[i+deg] = ai.Mul(p.inner[i])
	}

	for i := range deg {
		prod[i] = p.f.ElemFromUint64(0)
	}

	return NewPolynomial(prod, p.isCoefficientMod)
}

func (p *Polynomial) Degree() int {
	return p.leadingCoeffPos()
}

func (p *Polynomial) LeadCoeff() Elem {
	if pos := p.leadingCoeffPos(); pos >= 0 {
		return p.inner[pos]
	}

	return p.f.ElemFromUint64(0)
}

func (p *Polynomial) leadingCoeffPos() int {
	for i := len(p.inner) - 1; i >= 0; i-- {
		if p.inner[i].Value() != 0 {
			return i
		}
	}

	return math.MinInt
}

func (p *Polynomial) removeLeadingZeroes() {
	if p.isCoefficientMod {
		return
	}

	lead := p.leadingCoeffPos()
	if lead < 0 {
		p.inner = []Elem{p.f.ElemFromUint64(0)}

		return
	}

	p.inner = p.inner[:lead+1]
}

func (p *Polynomial) Copy() *Polynomial {
	innercopy := make([]Elem, len(p.inner))
	for i := range p.inner {
		innercopy[i] = p.inner[i].Copy()
	}

	return NewPolynomial(innercopy, p.isCoefficientMod)
}

// todo: fix
func (p *Polynomial) String() string {
	p.removeLeadingZeroes()

	if len(p.inner) == 1 {
		return strconv.FormatUint(p.inner[0].Value(), 10)
	}

	bldr := strings.Builder{}

	for i := len(p.inner) - 1; i >= 0; i-- {
		if p.inner[i].Value() == 0 {
			continue
		}

		strI := strconv.FormatInt(int64(i), 10)

		strElem := strconv.FormatUint(p.inner[i].Value(), 10)
		bldr.WriteString(strElem)

		if i != 0 {
			bldr.WriteString("*x^")
			bldr.WriteString(strI)
			bldr.WriteString(" + ")
		}
	}

	return bldr.String()
}

// returns r= gcd(a,b), x, y such that ax + by = r.
// where r.Degree() < stopDegree.
func PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	if a.Degree() < stopDegree {
		gcd = a.Copy()
		x = a.f.constantPolynomial(1)
		y = a.f.constantPolynomial(0)

		return
	}

	quotient, r := a.LongDiv(b)
	gcd, x1, y1 := PartialExtendedEuclidean(b, r, stopDegree)
	x = y1
	y = x1.Sub((quotient).Mul(y1))

	return gcd, x, y
}

func (p *Polynomial) ToSlice() []uint64 {
	list := make([]uint64, len(p.inner))
	for i, e := range p.inner {
		list[i] = e.Value()
	}

	return list
}

// returns self for chaining/ fluent interface.
func (p *Polynomial) MulScalarInPlace(s Elem) *Polynomial {
	for i := range p.inner {
		p.inner[i] = p.inner[i].Mul(s)
	}

	return p
}

func (p *Polynomial) IsCoeffMode() bool {
	return p.isCoefficientMod
}
