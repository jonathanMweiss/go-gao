package field

import (
	"math"
	"strconv"
	"strings"
)

type Polynomial struct {
	f                Field
	inner            []uint64
	isCoefficientMod bool
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
		inner:            inner,
		isCoefficientMod: isPointRepresentation,
		f:                f,
	}
}

func preOpVerification(p, q *Polynomial) bool {
	if p.f.Modulus() != q.f.Modulus() {
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
	return len(p.inner) == 0 || p.f.Equals(p.inner[0], 0)
}

func (p *Polynomial) Add(q *Polynomial) *Polynomial {
	if !preOpVerification(p, q) {
		return nil // this is an error.
	}

	size := len(p.inner)
	if len(q.inner) > size {
		size = len(q.inner)
	}

	inner := make([]uint64, size)

	copy(inner, p.inner)

	for i := range q.inner {
		inner[i] = p.f.Add(inner[i], q.inner[i])
	}

	return NewPolynomial(p.f, inner, p.isCoefficientMod)
}
func (p *Polynomial) Sub(q *Polynomial) *Polynomial {
	if !preOpVerification(p, q) {
		return nil // this is an error.
	}

	negative := make([]uint64, len(q.inner))
	for i, v := range q.inner {
		negative[i] = p.f.Neg(v)
	}

	sub := p.Add(NewPolynomial(p.f, negative, q.isCoefficientMod))
	sub.removeLeadingZeroes()

	return sub
}

func (p *Polynomial) Mul(q *Polynomial) *Polynomial {
	if !preOpVerification(p, q) {
		return nil
	}
	fld := p.f

	var inner []uint64
	if p.isCoefficientMod {
		inner = make([]uint64, len(p.inner))
		for i := range p.inner {
			inner[i] = fld.Mul(p.inner[i], q.inner[i])
		}
	} else {
		// regular polynomial multiplication O(n^2)
		inner = make([]uint64, len(p.inner)+len(q.inner)-1)

		for i := range p.inner {
			for j := range q.inner {
				inner[i+j] = fld.Add(inner[i+j], fld.Mul(p.inner[i], q.inner[j]))
			}
		}
	}

	prod := NewPolynomial(fld, inner, p.isCoefficientMod)
	if !p.isCoefficientMod {
		prod.removeLeadingZeroes()
	}

	return prod
}

func (p *Polynomial) Eval(x uint64) uint64 {
	result := uint64(0)
	fld := p.f
	// horner's rule:
	for i := len(p.inner) - 1; i >= 0; i-- {
		result = fld.Mul(x, result)
		result = fld.Add(p.inner[i], result)
	}

	return result
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

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (p *Polynomial) LongDiv(v *Polynomial) (q *Polynomial, r *Polynomial) {
	if !preOpVerification(p, v) {
		return nil, nil
	}
	fld := p.f

	if v.isCoefficientMod {
		return nil, nil
	}

	n, m := p.Degree(), v.Degree()

	b := v.Copy()
	u := fld.Inverse(v.LeadCoeff())

	r = p.Copy()
	qInner := make([]uint64, n-m+1)

	for i := n - m; i >= 0; i-- {
		if r.Degree() == m+i {
			qInner[i] = fld.Mul(r.LeadCoeff(), u)
			r = r.Sub(monomialMultPoly(qInner[i], i, b))
		} else {
			qInner[i] = 0
		}
	}

	r.removeLeadingZeroes()

	if len(qInner) == 0 {
		qInner = []uint64{0}
	}

	q = NewPolynomial(fld, qInner, false)
	q.removeLeadingZeroes()

	return q, r
}

func monomialMultPoly(ai uint64, deg int, p *Polynomial) *Polynomial {
	newDegree := len(p.inner) + deg
	fld := p.f
	prod := make([]uint64, newDegree)

	for i := range p.inner {
		prod[i+deg] = fld.Mul(ai, p.inner[i])
	}

	for i := range deg {
		prod[i] = 0
	}

	return NewPolynomial(fld, prod, p.isCoefficientMod)
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
	if p.isCoefficientMod {
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
	for i := range p.inner {
		innercopy[i] = p.inner[i]
	}

	return NewPolynomial(p.f, innercopy, p.isCoefficientMod)
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

func makeConstantPoly(f Field, u uint64) *Polynomial {
	return NewPolynomial(f, []uint64{u}, false)
}

// returns r= gcd(a,b), x, y such that ax + by = r.
// where r.Degree() < stopDegree.
func PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	if a.Degree() < stopDegree {
		gcd = a.Copy()
		x = makeConstantPoly(a.f, 1)
		y = makeConstantPoly(a.f, 0)

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
		list[i] = e
	}

	return list
}

// returns self for chaining/ fluent interface.
func (p *Polynomial) MulScalarInPlace(s uint64) *Polynomial {
	fld := p.f
	for i := range p.inner {
		p.inner[i] = fld.Mul(p.inner[i], s)
	}

	return p
}

func (p *Polynomial) IsCoeffMode() bool {
	return p.isCoefficientMod
}

// PolyProductMonicNegRoots computes \prod (x - r_i).
func PolyProductMonicNegRoots(f Field, roots []uint64) *Polynomial {
	n := len(roots)
	if n == 0 {
		return makeConstantPoly(f, 1)
	}

	coeffs := make([]uint64, n+1)
	coeffs[0] = 1

	deg := 0
	for _, r := range roots {
		neg := f.Neg(f.Reduce(r)) // -r mod p
		coeffs[deg+1] = 0
		for j := deg; j >= 0; j-- {
			// new[j+1] += old[j] * 1
			coeffs[j+1] = f.Add(coeffs[j+1], coeffs[j])
			// new[j]   += old[j] * (-r)
			coeffs[j] = f.Mul(coeffs[j], neg)
		}
		deg++
	}

	out := make([]uint64, deg+1)
	for i := 0; i <= deg; i++ {
		out[i] = coeffs[i]
	}

	return &Polynomial{f: f, inner: out, isCoefficientMod: false}
}
