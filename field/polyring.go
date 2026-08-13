// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/bits"
	"sync"
)

// A PolyRing performs polynomial arithmetic over a fixed coefficient field.
type PolyRing interface {
	GetField() Field

	// NewPolynomial builds a polynomial over this ring's field. Coefficients run from
	// lowest to highest degree; an empty slice yields the zero polynomial.
	NewPolynomial(inner []uint64, isPointRepresentation bool) *Polynomial

	Evaluate(a *Polynomial, x uint64) uint64

	// Assumes polynomial of valid degree.
	NttForward(a *Polynomial) error
	NttBackward(a *Polynomial) error

	// compute c = a + b
	Add(a, b, c *Polynomial)
	// compute c = a - b
	Sub(a, b, c *Polynomial)

	// compute c = a * scalar
	MulScalar(a *Polynomial, scalar uint64, c *Polynomial)

	// compute c = a * b
	// performs smart dispatch between schoolbook and NTT based
	// on size and NTT support of the inner field.
	Mul(a, b, c *Polynomial)

	// Product multiplies a whole slice via a divide-and-conquer product tree,
	// in O(n log^2 n). An empty slice yields the constant polynomial p(x) = 1.
	Product(polys []*Polynomial) *Polynomial

	// Creates quotient q and remainder r.
	// chooses the algorithm based on size and NTT support of the inner field.
	Div(a, b *Polynomial) (q *Polynomial, r *Polynomial)

	// Extended Euclidean algorithm.
	PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)

	// Uses the half-GCD algorithm, less suitable for small inputs but asymptotically faster than PartialExtendedEuclidean.
	FastPartialGCD(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)
}

const nttMulThreshold = 16 // ~coeff count where NTT starts winning

// DensePolyRing implements PolyRing with optional NTT domain for polynomials.
type DensePolyRing struct {
	Field
	mu           sync.RWMutex
	twiddleCache map[int]*twiddleSet // key: n
}

type nttMulScratch struct {
	a []uint64
	b []uint64
}

var nttMulScratchPool = sync.Pool{
	New: func() any {
		return &nttMulScratch{}
	},
}

// NewDensePolyRing constructs a ring over the provided coefficient field.
//
// It panics if f is nil.
func NewDensePolyRing(f Field) PolyRing {
	if f == nil {
		panic("NewDensePolyRing: nil field")
	}

	return &DensePolyRing{
		Field:        f,
		twiddleCache: map[int]*twiddleSet{},
	}
}

// GetField returns the coefficient field this ring operates over.
func (r *DensePolyRing) GetField() Field { return r.Field }

// NewPolynomial builds a polynomial over this ring's field from coefficients ordered
// from lowest to highest degree. An empty inner yields the zero polynomial.
//
// The polynomial takes ownership of inner rather than copying it.
func (r *DensePolyRing) NewPolynomial(inner []uint64, isPointRepresentation bool) *Polynomial {
	// An empty sum of terms is zero.
	if len(inner) == 0 {
		inner = []uint64{0}
	}

	return &Polynomial{
		inner: inner,
		isNTT: isPointRepresentation,
		f:     r.Field,
	}
}

// Product multiplies every polynomial in polys using a divide-and-conquer product tree,
// in O(n log^2 n) with NTT-based multiplication rather than the O(n^2) of a running
// product.
//
// An empty slice yields the constant polynomial p(x) = 1
//
// polys is only read; the result never aliases any of its elements.
func (r *DensePolyRing) Product(polys []*Polynomial) *Polynomial {
	switch len(polys) {
	case 0:
		return makeConstantPoly(r.Field, 1)
	case 1:
		// Copy, so the caller cannot mutate the result through the input slice.
		return polys[0].Copy()
	default:
		return r.productTree(polys)
	}
}

// productTree recurses over polys. At a leaf it hands back the caller's polynomial
// uncopied, which is safe because Mul only ever reads its operands: the NTT path copies
// both into scratch buffers before transforming, and the schoolbook path writes solely
// to its output. Every internal node therefore returns a freshly allocated polynomial.
func (r *DensePolyRing) productTree(polys []*Polynomial) *Polynomial {
	if len(polys) == 1 {
		return polys[0]
	}

	mid := len(polys) / 2
	left := r.productTree(polys[:mid])
	right := r.productTree(polys[mid:])

	res := &Polynomial{}
	r.Mul(left, right, res)

	return res
}

// ---------- utilities ----------

func ensureLen(c *Polynomial, n int) {
	if len(c.inner) < n {
		tmp := make([]uint64, n)
		copy(tmp, c.inner)
		c.inner = tmp
	} else {
		c.inner = c.inner[:n]
	}
}

func ensureLenCheap(c *Polynomial, n int) {
	if len(c.inner) < n {
		c.inner = make([]uint64, n)
	} else {
		c.inner = c.inner[:n]
	}
}

func (r *DensePolyRing) trimTrailingZeros(p *Polynomial) {
	if len(p.inner) == 0 || p.isNTT {
		// In NTT domain we keep the fixed size.
		return
	}

	i := len(p.inner) - 1
	for i >= 0 && r.Equals(p.inner[i], 0) {
		i--
	}
	p.inner = p.inner[:i+1]
}

// ---------- Poly ops ----------

// Evaluate returns a(x) by Horner's rule.
// should receive a polynomial in coefficient form, not NTT form.
func (r *DensePolyRing) Evaluate(a *Polynomial, x uint64) uint64 {
	result := uint64(0)
	fld := r.Field

	// horner's rule:
	for i := len(a.inner) - 1; i >= 0; i-- {
		result = fld.Add(a.inner[i], fld.Mul(x, result))
	}

	return result
}

// MulScalar computes c = a * scalar, preserving a's domain. c may alias a.
func (r *DensePolyRing) MulScalar(a *Polynomial, scalar uint64, c *Polynomial) {
	s := r.Reduce(scalar)
	f := r.GetField()

	ensureLen(c, len(a.inner))
	for i := range a.inner {
		c.inner[i] = f.Mul(a.inner[i], s)
	}

	c.f = r.Field
	c.isNTT = a.isNTT // scalar mult preserves domain

	r.trimTrailingZeros(c)
}

// Add computes c = a + b. c may alias a or b.
//
// It panics if a and b are over different fields or in different domains.
func (r *DensePolyRing) Add(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)

	f := r.Field

	var av, bv uint64
	for i := 0; i < n; i++ {
		if i < alen {
			av = r.Reduce(a.inner[i])
		} else {
			av = 0
		}

		if i < blen {
			bv = r.Reduce(b.inner[i])
		} else {
			bv = 0
		}

		c.inner[i] = f.Add(av, bv)
	}

	c.f = r.Field
	c.isNTT = a.isNTT
	r.trimTrailingZeros(c)
}

// Sub computes c = a - b. c may alias a or b.
//
// It panics if a and b are over different fields or in different domains.
func (r *DensePolyRing) Sub(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	c.f = r.Field
	c.isNTT = a.isNTT

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)
	minLen := min(alen, blen)

	f := r.Field

	// Subtract overlapping part
	for i := 0; i < minLen; i++ {
		c.inner[i] = f.Sub(a.inner[i], b.inner[i])
	}

	// Handle remaining coefficients
	if alen > blen {
		// If a is longer, copy its remaining part
		copy(c.inner[minLen:], a.inner[minLen:])
	} else if blen > alen {
		// If b is longer, copy the negation of its remaining part
		for i := minLen; i < blen; i++ {
			c.inner[i] = f.Neg(b.inner[i])
		}
	}

	r.trimTrailingZeros(c)
}

// Mul computes c = a * b, choosing between schoolbook and NTT-based multiplication by size.
//
// It panics if a and b are over different fields or in different domains.
func (r *DensePolyRing) Mul(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	la, lb := len(a.inner), len(b.inner)
	if la == 0 || lb == 0 {
		c.f, c.inner, c.isNTT = r.Field, c.inner[:0], a.isNTT
		return
	}

	// In NTT domain, multiplication is pointwise and preserves NTT representation.
	// (preOpVerification ensures a.isNTT == b.isNTT)
	if a.isNTT {
		n := len(a.inner)
		ensureLen(c, n)
		r.pointwiseMult(a, b, c)

		c.f = r.Field
		c.isNTT = true
		return
	}

	// Coefficient-domain smart dispatch: schoolbook for small sizes, NTT otherwise.
	if min(la, lb) <= nttMulThreshold || !r.canUseNTTConvolutionLen(la+lb-1) {
		r.mulSchoolbook(a, b, c)
		return
	}

	r.mulViaNTT(a, b, c)
}

func (r *DensePolyRing) mulSchoolbook(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	if a.isNTT || b.isNTT {
		panic("mulSchoolbook cannot handle NTT polynomials")
	}

	f := r.Field

	newLen := len(a.inner) + len(b.inner) - 1

	// Decide where to write: use c.inner if capacity is enough; else allocate.
	var out []uint64
	if cap(c.inner) >= newLen {
		out = c.inner[:newLen]
		clear(out)

	} else {
		out = make([]uint64, newLen)
	}

	// Perform schoolbook convolution: O(n*m).
	// out[i+j] += a[i] * b[j]
	for i := range a.inner {
		ai := a.inner[i]
		if ai == 0 {
			continue
		}

		for j := range b.inner {
			out[i+j] = f.Add(out[i+j], f.Mul(ai, b.inner[j]))
		}
	}

	// Write result into c (safe even if c==a or c==b because we used `out`).
	c.f = a.f
	c.inner = out
	c.isNTT = false

	r.trimTrailingZeros(c)
}

func (r *DensePolyRing) monomialMultPoly(ai uint64, deg int, p *Polynomial) *Polynomial {
	newLen := len(p.inner) + deg
	fld := r.GetField()
	prod := make([]uint64, newLen)

	for i := range p.inner {
		prod[i+deg] = fld.Mul(ai, p.inner[i])
	}

	return r.NewPolynomial(prod, p.isNTT)
}

// Div returns the quotient and remainder of a divided by b, choosing between schoolbook
// and NTT-based division by size.
//
// similar to standard division, this function panics when its input doesn't make sense:
// a or b must not be nil, or in the NTT domain; they must have the same field, and
// b must not be the zero polynomial.
func (r *DensePolyRing) Div(a, b *Polynomial) (q *Polynomial, rem *Polynomial) {
	if a == nil || b == nil {
		panic("Div: nil polynomial")
	}

	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	if a.isNTT || b.isNTT {
		panic("Div expects coefficient-domain polynomials")
	}

	bDeg := b.Degree()
	aDeg := a.Degree()

	if bDeg < 0 {
		panic("division by zero polynomial")
	}

	if aDeg < bDeg {
		return polyZero(r.Field), a.Copy()
	}

	quotLen := aDeg - bDeg + 1
	// Since divViaNTT uses NTT-based convolution under the hood,
	// we check if the convolution size is supported before deciding which division algorithm to use.
	// 2*quotLen-1 is the convolution size needed
	maxConvLen := max(2*quotLen-1, aDeg+1)

	if quotLen > nttMulThreshold && r.canUseNTTConvolutionLen(maxConvLen) {
		return r.divViaNTT(a, b)
	}

	return r.divSchoolbook(a, b)
}

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (r *DensePolyRing) divSchoolbook(a, b *Polynomial) (q *Polynomial, rem *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}
	fld := r.Field

	if b.isNTT {
		panic("divSchoolbook expects coefficient-domain polynomials")
	}

	n, m := a.Degree(), b.Degree()

	u := fld.Inverse(b.LeadCoeff()) // Assumes inverse exists.

	rem = a.Copy()
	qInner := make([]uint64, n-m+1)

	for i := n - m; i >= 0; i-- {
		// TODO: keeping the degree in a variable might save time.
		if rem.Degree() == m+i {
			qInner[i] = fld.Mul(rem.LeadCoeff(), u)
			r.Sub(rem, r.monomialMultPoly(qInner[i], i, b), rem)
		} else {
			qInner[i] = 0
		}
	}

	r.trimTrailingZeros(rem)

	q = r.NewPolynomial(qInner, false)
	q.removeLeadingZeroes()

	return q, rem
}

func makeConstantPoly(f Field, u uint64) *Polynomial {
	return &Polynomial{f: f, inner: []uint64{u}, isNTT: false}
}

func polyZero(f Field) *Polynomial {
	return makeConstantPoly(f, 0)
}

func polyOne(f Field) *Polynomial {
	return makeConstantPoly(f, 1)
}

type polyMatrix2x2 struct {
	a00, a01 *Polynomial
	a10, a11 *Polynomial
}

func (m polyMatrix2x2) Mul(r *DensePolyRing, other polyMatrix2x2) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: polyAdd(r, polyMul(r, m.a00, other.a00), polyMul(r, m.a01, other.a10)),
		a01: polyAdd(r, polyMul(r, m.a00, other.a01), polyMul(r, m.a01, other.a11)),
		a10: polyAdd(r, polyMul(r, m.a10, other.a00), polyMul(r, m.a11, other.a10)),
		a11: polyAdd(r, polyMul(r, m.a10, other.a01), polyMul(r, m.a11, other.a11)),
	}
}

func (m polyMatrix2x2) MulVec(r *DensePolyRing, a, b *Polynomial) (*Polynomial, *Polynomial) {
	aOut := polyAdd(r, polyMul(r, m.a00, a), polyMul(r, m.a01, b))
	bOut := polyAdd(r, polyMul(r, m.a10, a), polyMul(r, m.a11, b))
	return aOut, bOut
}

func polyIdentity2x2(f Field) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: makeConstantPoly(f, 1), // x0
		a01: makeConstantPoly(f, 0), // y0
		a10: makeConstantPoly(f, 0), // x1
		a11: makeConstantPoly(f, 1), // y1
	}
}

// Creates:
// |1 0 |
// |1 -q|
func stepMatrix(r *DensePolyRing, q *Polynomial) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: polyZero(r.Field), a01: polyOne(r.Field),
		a10: polyOne(r.Field), a11: polySub(r, polyZero(r.Field), q),
	}
}

// PartialExtendedEuclidean runs the extended Euclidean algorithm, stopping early.
//
// returns r= gcd(a,b), x, y such that ax + by = r.
// where r.Degree() < stopDegree. For full GCD, use stopDegree=0.
func (r *DensePolyRing) PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// Work on local copies ensuring inputs aren't mutated.
	A := a.Copy()
	B := b.Copy()
	degA := A.Degree()
	degB := B.Degree()

	// M is the matrix that begins as identity, and each iteration it is updated by left-multiplying the Bézout matrix of the current division step:
	// M =  a00 a01 = 1 0
	//      a10 a11   0 1
	M := polyIdentity2x2(r.Field)

	// Reusable temporaries (avoid allocations).
	tmp1 := &Polynomial{f: r.Field} // holds q*M10 or q*M11
	tmp2 := &Polynomial{f: r.Field} // holds M00 - q*M10 or M01 - q*M11

	// Bezout identity is GCD(A,B)= ax +by.
	// letsdive into the iterative algorithm.
	// Initially, we have A=a, B=b, M = I, so the invariant holds:
	// | A |   | 1  0 |   | a |
	// | B | = | 0  1 | * | b |
	// In each step of the Euclidean algorithm, we perform the division $a = q \cdot b + r$, which implies $r = a - q \cdot b$.
	// Thus,
	// | B |   | 0  1 |   | a |   | b |
	// | r | = |1  -q| *  | b | = | a - q*b |
	// Thus,
	// each step performs left-multiplication by the Bézout matrix of the current division step:
	// | 0   1 |   | M00  M01 |   | M10            M11           |
	// |1   -q | * | M10  M11 | = | M00 - q*M10    M01 - q*M11   |
	for degA >= stopDegree {
		// If B == 0, can't divide further.
		if degB < 0 {
			break
		}

		// A = q*B + r
		q, rrem := r.Div(A, B)
		A, B = B, rrem // GCD recursive step: gcd(A, B) = gcd(B,rrem)
		degA, degB = degB, B.Degree()

		// Performing the matrix update one part at a time, to reuse the same temporary
		// for both x and y updates and avoid extra allocations.

		// left matrix update: (M00, M10) = (M10, M00 - q*M10)
		r.Mul(q, M.a10, tmp1)    // tmp1 = q * M10
		r.Sub(M.a00, tmp1, tmp2) // tmp2 = M00 - q*M10
		M.a00, M.a10, tmp2 = M.a10, tmp2, M.a00

		// right update: (M01, M11) = (M11, M01 - q*M11)
		r.Mul(q, M.a11, tmp1)    // tmp1 = q * M11
		r.Sub(M.a01, tmp1, tmp2) // tmp2 = M01 - q*M11
		M.a01, M.a11, tmp2 = M.a11, tmp2, M.a01
	}

	// Notice that in each step we updated A,B using the GCD step, so now B=0, A=GCD.
	// we can extract the coeffiecnts x,y from the matrix M,
	// since we maintained the invariant that A = M.a00*a + M.a01*b.
	// Namely, A is the top-left position in our vector, which is the GCD,
	// and M00, M01 are what we multiplied against a,b to get the result A.
	return A, M.a00, M.a01
}

// Reverse the top L coefficients: rev_L(f) = x^{L-1} * f(1/x) truncated to L.
// Uses the *true* degree (last non-zero) rather than len(inner)-1.
func (r *DensePolyRing) revTop(f *Polynomial, L int) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	if L <= 0 {
		return out
	}
	out.inner = make([]uint64, L)

	// Find true degree (ignore trailing zeros)
	n := len(f.inner) - 1
	for n >= 0 && r.Equals(f.inner[n], 0) {
		n--
	}
	if n < 0 {
		// zero polynomial
		return out
	}

	// b[i] = a[n - i] if n-i >= 0
	for i := 0; i < L; i++ {
		j := n - i
		if j >= 0 {
			out.inner[i] = r.Reduce(f.inner[j])
		} else {
			out.inner[i] = 0
		}
	}
	return out
}

func nextPow2(n int) int {
	return 1 << (bits.Len(uint(n - 1)))
}

func (r *DensePolyRing) mulViaNTT(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	// else, use mulTrunc with total length (coeff-domain out)
	la, lb := len(a.inner), len(b.inner)
	if la == 0 || lb == 0 {
		c.f, c.inner, c.isNTT = r.Field, c.inner[:0], false
		return
	}

	// if both are NTT, pointwise multiply
	if a.isNTT && b.isNTT {
		n := len(a.inner)
		ensureLen(c, n)
		r.pointwiseMult(a, b, c)

		c.f = r.Field
		c.isNTT = true

		return
	}

	total := la + lb - 1
	prod := r.mulTrunc(a, b, total) // NTT under the hood, coeff-domain out
	c.inner = prod.inner
	c.f, c.isNTT = r.Field, false
}

func (r *DensePolyRing) pointwiseMult(a, b, c *Polynomial) {
	f := r.Field
	for i := range c.inner {
		c.inner[i] = f.Mul(a.inner[i], b.inner[i])
	}
}

func (r *DensePolyRing) getNttMulScratch(n int) *nttMulScratch {
	s := nttMulScratchPool.Get().(*nttMulScratch)

	if cap(s.a) < n {
		s.a = make([]uint64, n)
	} else {
		s.a = s.a[:n]
		clear(s.a)
	}

	if cap(s.b) < n {
		s.b = make([]uint64, n)
	} else {
		s.b = s.b[:n]
		clear(s.b)
	}

	return s
}

func (r *DensePolyRing) putNttMulScratch(s *nttMulScratch) {
	nttMulScratchPool.Put(s)
}

// Multiply polynomials and then truncate to the lowest L terms.
// Use NTT under the hood (size = nextPow2(L + L - 1)), then slice [:L].
func (r *DensePolyRing) mulTrunc(a, b *Polynomial, L int) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	r.mulTruncInto(out, a, b, L)
	return out
}

func (r *DensePolyRing) mulTruncInto(dst *Polynomial, a, b *Polynomial, L int) {
	dst.f = r.Field
	dst.isNTT = false

	if L <= 0 {
		dst.inner = dst.inner[:0]
		return
	}

	if cap(dst.inner) < L {
		dst.inner = make([]uint64, L)
	} else {
		dst.inner = dst.inner[:L]
		clear(dst.inner)
	}

	if a == nil || b == nil {
		return
	}

	la := min(len(a.inner), L)
	lb := min(len(b.inner), L)
	if la == 0 || lb == 0 {
		return
	}

	total := la + lb - 1
	convLen := min(L, total)
	n := nextPow2(total)

	scratch := r.getNttMulScratch(n)
	defer r.putNttMulScratch(scratch)

	copy(scratch.a[:la], a.inner[:la])
	copy(scratch.b[:lb], b.inner[:lb])

	aNTT := &Polynomial{f: r.Field, inner: scratch.a, isNTT: false}
	bNTT := &Polynomial{f: r.Field, inner: scratch.b, isNTT: false}

	if err := r.NttForward(aNTT); err != nil {
		panic(err)
	}
	if err := r.NttForward(bNTT); err != nil {
		panic(err)
	}

	r.pointwiseMult(aNTT, bNTT, aNTT)

	if err := r.nttBackwardNoTrim(aNTT); err != nil {
		panic(err)
	}

	copy(dst.inner, aNTT.inner[:convLen])
}

// Series inverse modulo x^k using Newton iteration.
// Assumes b[0] != 0; returns t such that (b * t) ≡ 1 mod x^k.
// seriesInverse looks for the root of F(T)= T(X)^1 - B(X) = 0 mod x^k.
// That is, using Newton iteration:
// T_{n+1} = T_n - F(T_n)/F'(T_n) = T_n - (T_n - B*T_n^2) = 2*T_n - B*T_n^2 mod x^k
// where T_1 = b[0]^{-1}.
//
// Preconditions:
//   - b is in coefficient domain (isNTT == false)
//   - k >= 1
//   - b.inner[0] != 0 (invertible constant term)
func (r *DensePolyRing) seriesInverse(b *Polynomial, k int) *Polynomial {
	if k <= 0 {
		return &Polynomial{f: r.Field, isNTT: false}
	}
	if len(b.inner) == 0 || r.Equals(b.inner[0], 0) {
		panic("seriesInverse: constant term is zero")
	}

	b0 := r.Reduce(b.inner[0])
	t := &Polynomial{f: r.Field, isNTT: false, inner: []uint64{r.Inverse(b0)}}
	tmp := &Polynomial{f: r.Field, isNTT: false}
	next := &Polynomial{f: r.Field, isNTT: false}
	two := r.Reduce(2)

	f := r.Field
	for l := 1; l < k; {
		m := l << 1
		if m > k {
			m = k
		}

		// tmp = b*t mod x^m
		r.mulTruncInto(tmp, b, t, m)

		// tmp = 2 - tmp (mod x^m)
		tmp.inner[0] = f.Sub(two, tmp.inner[0])
		for i := 1; i < m; i++ {
			tmp.inner[i] = f.Neg(tmp.inner[i])
		}

		// t = t * tmp mod x^m
		r.mulTruncInto(next, t, tmp, m)
		t, next = next, t
		l = m
	}
	return t
}

// DivNTT follows `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard, section 9.1.
//
// The algorithm capitalize the notion that Newton Iteration can be used to compute the inverse of
// a polynomial in O(n log n) time.
func (r *DensePolyRing) divViaNTT(a, b *Polynomial) (q, rem *Polynomial) {
	if a == nil || b == nil || a.isNTT || b.isNTT {
		panic("LongDivNTT expects non-nil coefficient-domain polynomials")
	}
	n := len(a.inner) - 1
	m := len(b.inner) - 1
	if m < 0 {
		panic("division by zero polynomial")
	}
	if n < m {
		// q = 0, r = a
		return &Polynomial{f: r.Field, isNTT: false, inner: []uint64{0}}, a.Copy()
	}

	k := n - m + 1 // quotient length

	// 1) Reverse tops
	Astar := r.revTop(a, k)   // length k
	Bstar := r.revTop(b, m+1) // length m+1

	// lead(b) maps to Bstar[0]; must be invertible
	if len(Bstar.inner) == 0 || r.Equals(Bstar.inner[0], 0) {
		panic("division by polynomial with zero leading coefficient")
	}

	// 2) T = (Bstar)^{-1} mod x^k (Newton series inverse)
	T := r.seriesInverse(Bstar, k) // length k

	// 3) Q* = A* * T mod x^k
	Qstar := &Polynomial{f: r.Field, isNTT: false}
	r.mulTruncInto(Qstar, Astar, T, k)

	// 4) q = rev_k(Q*)
	q = r.revTop(Qstar, k) // coefficient domain

	// 5) rem = a − q*b
	prod := &Polynomial{f: r.Field, isNTT: false}
	r.mulTruncInto(prod, q, b, n+1) // full product length (deg = n)
	rem = &Polynomial{f: r.Field, isNTT: false}
	r.Sub(a, prod, rem)      // coeff-domain subtraction
	r.trimTrailingZeros(rem) // ensure deg(rem) < deg(b)

	return q, rem
}

func (r *DensePolyRing) canUseNTTConvolutionLen(convLen int) bool {
	if convLen <= 0 {
		return false
	}

	n := nextPow2(convLen)
	modMinusOne := r.Modulus() - 1
	if uint64(n) > modMinusOne {
		return false
	}

	return modMinusOne%uint64(n) == 0
}

func polyMul(r *DensePolyRing, a, b *Polynomial) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	r.Mul(a, b, out)
	return out
}

func polyAdd(r *DensePolyRing, a, b *Polynomial) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	r.Add(a, b, out)
	return out
}

func polySub(r *DensePolyRing, a, b *Polynomial) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	r.Sub(a, b, out)
	return out
}

// mulSubInto computes dst = a - q*b in one fused pass when q is small,
// avoiding a temporary allocation for the intermediate product q*b.
// When q is large, it falls back to mulFull + Sub.
// All polynomials must be in coefficient domain.
func (r *DensePolyRing) mulSubInto(dst, a, q, b *Polynomial) {
	lq := len(q.inner)
	lb := len(b.inner)

	// If q or b is empty, dst = a.
	if lq == 0 || lb == 0 {
		dst.f = r.Field
		dst.isNTT = false
		ensureLen(dst, len(a.inner))
		copy(dst.inner, a.inner)
		return
	}

	// For large q, fall back to mulFull + Sub with a temporary.
	if min(lq, lb) > nttMulThreshold {
		tmp := &Polynomial{f: r.Field}
		r.Mul(q, b, tmp)
		r.Sub(a, tmp, dst)
		return
	}

	// Fused schoolbook: dst[k] = a[k] - sum_{i+j=k} q[i]*b[j]
	// Product q*b has length lq + lb - 1.
	prodLen := lq + lb - 1
	n := max(len(a.inner), prodLen)

	f := r.Field
	dst.f = f
	dst.isNTT = false
	ensureLenCheap(dst, n)

	// Start with a copy of `a`, zero-extended.
	la := len(a.inner)
	if la > 0 {
		copy(dst.inner[:la], a.inner)
	}
	for i := la; i < n; i++ {
		dst.inner[i] = 0
	}

	// Subtract q*b from dst in-place: dst[i+j] -= q[i] * b[j]
	// Iterate over q (the small operand) in the outer loop.
	for i := 0; i < lq; i++ {
		qi := q.inner[i]
		if qi == 0 {
			continue
		}
		for j := 0; j < lb; j++ {
			dst.inner[i+j] = f.Sub(dst.inner[i+j], f.Mul(qi, b.inner[j]))
		}
	}

	r.trimTrailingZeros(dst)
}

// NttPartialExtendedEuclidean is PartialExtendedEuclidean using NTT-based division steps.
//
// for full explanation on the iterative algorithm go to PartialExtendedEuclidean.
// This is the same algorithm but uses schoolbook/DivNTT steps.
func (r *DensePolyRing) NttPartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// Work on local copies ensuring inputs aren't mutated (coeff domain expected).
	A := a.Copy()
	B := b.Copy()
	A.isNTT, B.isNTT = false, false

	M := polyIdentity2x2(r.Field)

	// Reusable temporaries (avoid allocations in the single-step path).
	tmp := &Polynomial{f: r.Field}

	degA := A.Degree()
	degB := B.Degree()

	for degA >= stopDegree {
		if degB < 0 || len(B.inner) == 0 {
			break
		}

		// A = q*B + r
		q, rrem := r.Div(A, B)
		A, B = B, rrem // gcd(A,B) = gcd(B,rrem)
		degA, degB = degB, B.Degree()

		// left side
		r.mulSubInto(tmp, M.a00, q, M.a01)
		M.a00, M.a01, tmp = M.a01, tmp, M.a00

		// right side
		r.mulSubInto(tmp, M.a10, q, M.a11)
		M.a10, M.a11, tmp = M.a11, tmp, M.a10
	}

	return A, M.a00, M.a10
}

/*
Half-GCD (HGCD) Implementation and Fast Partial GCD

The following functions implement the O(n log^2 n) version of the Euclidean Algorithm.
Traditional Euclidean Algorithm computes remainders one by one, taking O(n^2) time.
HGCD uses a divide-and-conquer approach by only looking at the most significant
coefficients (the "high parts") to compute transition matrices.

Terminology:
- Euclidean sequence: The sequence of remainders r_0, r_1, ... where r_0=a, r_1=b.
- Transition Matrix M: A 2x2 matrix such that [r_i, r_{i+1}]^T = M * [a, b]^T.
*/

// hgcdThreshold is the degree at which the half-GCD recursion bottoms out into the
// iterative routine, below which the recursion overhead exceeds its benefit.
//
// 256 chosen empirically: measured best or tied-best against 128, 512, 1024 and 4096
// at n = 2048, 8192 and 32768, worth 2-8% over the previous 128.
const hgcdThreshold = 256

/*
FastPartialGCD finds the first remainder in the Euclidean sequence of (a, b)
whose degree is strictly less than stopDegree.

Works on Copies of A and B.

Complexity: O(n log^2 n).
*/
func (r *DensePolyRing) FastPartialGCD(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// We work on copies to preserve the original inputs.
	A := a.Copy()
	B := b.Copy()

	r.ensureNotNttForm(A)
	r.ensureNotNttForm(B)

	// fastGCDRec is the recursive driver that uses HGCD to 'jump' through the sequence.
	M := r.fastGCDRec(A, B, stopDegree)

	// Apply the final transition matrix to the original inputs to get the desired remainder.
	AOut, _ := M.MulVec(r, A, B)
	r.trimTrailingZeros(AOut)

	// Return the remainder and its corresponding Bézout coefficients for 'a' and 'b':
	// AOut = M.a00 * a + M.a01 * b
	return AOut, M.a00, M.a01
}

func (r *DensePolyRing) ensureNotNttForm(A *Polynomial) {
	if A.isNTT {
		// Ignoring this error would clear isNTT below on a polynomial still holding
		// evaluations, silently mislabelling it as coefficients.
		if err := r.NttBackward(A); err != nil {
			panic("ensureNotNttForm: " + err.Error())
		}
	}

	A.isNTT = false
}

/*
fastGCDRec is the recursive engine for Fast GCD. It bridges the gap between
the starting degrees and the target stopDegree using HGCD for large steps.

Parameters:
- a, b: Current polynomials in the sequence.
- stopDegree: The degree boundary we are aiming to cross.
*/
func (r *DensePolyRing) fastGCDRec(a, b *Polynomial, stopDegree int) polyMatrix2x2 {
	// Base Case 1: Target reached.
	aDeg := a.Degree()
	if aDeg < stopDegree || b.Degree() < 0 {
		return polyIdentity2x2(r.Field)
	}

	// Base Case 2: Small polynomials, use iterative O(n^2) logic.
	if aDeg < hgcdThreshold {
		return r.iterativePartialExtendedEuclideanMatrix(a, b, stopDegree)
	}

	// Calculate required degree reduction 'm'.
	// We want to reduce deg(a) until it is < stopDegree.
	// Distance to target = deg(a) - (stopDegree - 1) = deg(a) - stopDegree + 1.
	n := aDeg
	m := n - stopDegree + 1

	// 1. Half-GCD Recursive Step:
	// Use HGCD to compute a matrix M that reduces degrees significantly.
	M := r.hgcd(a, b, m)
	aCur, bCur := M.MulVec(r, a, b)
	r.trimTrailingZeros(aCur)
	r.trimTrailingZeros(bCur)

	// Check if the HGCD step was enough to reach the target stopDegree.
	if aCur.Degree() < stopDegree || bCur.Degree() < 0 {
		return M
	}

	// 2. Standard Euclidean Step:
	// Perform exactly ONE division step: aCur = q*bCur + rem.
	// This step is mandatory to ensure progress. Without it, the algorithm
	// might call HGCD with the same parameters again, leading to an infinite loop.
	q, rem := r.Div(aCur, bCur)
	Mstep := stepMatrix(r, q)
	M = Mstep.Mul(r, M)

	// After one division step, the pair is (bCur, rem).
	// If bCur.Degree() < stopDegree, then bCur is the first remainder with degree < stopDegree.
	if bCur.Degree() < stopDegree || rem.Degree() < 0 {
		return M
	}

	// 3. Second Recursive Step:
	// Continue the process on the remainder.
	S := r.fastGCDRec(bCur, rem, stopDegree)
	return S.Mul(r, M)
}

/*
The half-GCD (HGCD) algorithm is a divide-and-conquer method that computes a GCD transition matrix.
HGCD returns a matrix M such that for (a', b') = M * (a, b),
the degree of b' is reduced by at least 'm' relative to the original degree of a.
That is: deg(b') < deg(a) - m.

The HGCD algorithm relies on an insight (proven in `Modern Computer Algebra` by
Joachim von zur Gathen and Jürgen Gerhard) that the high-order coefficients of
the polynomials share the same transition matrix as the full polynomials, allowing us
to use only the top coefficients of a and b to compute a transition matrix that will also reduce
the full polynomials by a significant amount.
Namely, it works by recursively applying itself to the "high parts" of the polynomials,
which are obtained by shifting the inputs to focus on the top coefficients. and thus performing
smaller multiplications and divisions as much as possible.

Procedure:
  - First Call: performs HGCD on the high parts of a and b, reducing the degree by about m/2.
  - Bridge: It performs one division to ensure progress.
  - Second Call: It calls itself recursively on the new pair (b,remainder) to advance m further.

The result
This implements the Schönhage strategy of high-part recursion.
*/
func (r *DensePolyRing) hgcd(a, b *Polynomial, m int) polyMatrix2x2 {
	n := a.Degree()
	// Base Case: target reduction reached or degree too small.
	if b.Degree() < n-m || n < hgcdThreshold {
		return r.iterativePartialExtendedEuclideanMatrix(a, b, n-m)
	}

	/*
		Divide and Conquer:
		To reduce by distance 'm', we first reduce by distance m/2.
		We 'shift' the polynomials by k to only look at the top bits.
		Shifting ensures that the recursive calls work on smaller polynomials (O(m) degrees)
		while the results remain valid for the high-order coefficients of the full inputs.
	*/
	m1 := m / 2
	k := n - 2*m1
	if k < 0 {
		k = 0
	}

	// 1. First Recursive Call (on high parts):
	// Reduction distance is m1.
	R := r.hgcd(r.shiftRight(a, k), r.shiftRight(b, k), m1)

	// 2. Apply the transition matrix R to the full inputs.
	aCur, bCur := R.MulVec(r, a, b)
	r.trimTrailingZeros(aCur)
	r.trimTrailingZeros(bCur)

	bCurDeg := bCur.Degree()
	// Check if the reduction target m has already been met.
	if bCurDeg < n-m || bCurDeg < 0 {
		return R
	}

	// 3. Standard Euclidean Step (Mandatory Progress):
	// Perform one division step to ensure the next recursive call makes progress.
	q, rem := r.Div(aCur, bCur)
	Mstep := stepMatrix(r, q)
	M := Mstep.Mul(r, R)

	remDeg := rem.Degree()
	// Check if the reduction target m is met after division.
	if remDeg < n-m || remDeg < 0 {
		return M
	}

	// 4. Second Recursive Call:
	// Calculate the remaining distance m2 to reach the total target m.
	n_new := bCur.Degree()
	m2 := m - (n - n_new)
	if m2 <= 0 {
		return M
	}
	// Recalculate shift relative to the new degree.
	k2 := n_new - 2*m2
	if k2 < 0 {
		k2 = 0
	}

	S := r.hgcd(r.shiftRight(bCur, k2), r.shiftRight(rem, k2), m2)
	return S.Mul(r, M)
}

/*
iterativePartialExtendedEuclideanMatrix is the O(n^2) fallback.
It computes the transition matrix M until deg(a) < stopDegree.
*/
func (r *DensePolyRing) iterativePartialExtendedEuclideanMatrix(a, b *Polynomial, stopDegree int) polyMatrix2x2 {
	A := a.Copy()
	B := b.Copy()
	M := polyIdentity2x2(r.Field)

	for A.Degree() >= stopDegree && B.Degree() >= 0 {
		q, rem := r.Div(A, B)
		A, B = B, rem

		Mnew := stepMatrix(r, q)
		M = Mnew.Mul(r, M)
	}
	return M
}

// shiftRight extracts the high coefficients of a polynomial by shifting.
func (r *DensePolyRing) shiftRight(p *Polynomial, m int) *Polynomial {
	if m <= 0 {
		return p.Copy()
	}
	if m >= len(p.inner) {
		return polyZero(r.Field)
	}
	return r.NewPolynomial(p.inner[m:], false)
}
