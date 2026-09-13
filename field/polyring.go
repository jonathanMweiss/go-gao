// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/bits"
	"sync"
)

const nttMulThreshold = 16 // ~coeff count where NTT starts winning

// A PolyRing performs polynomial arithmetic over a fixed coefficient field.
//
// Most operations write into a caller-supplied destination rather than allocating, so
// that a loop can reuse one polynomial: Add, Sub, Mul and MulScalar all take c as their
// last argument, and c may alias an input. Several dispatch on size; Mul and Div pick
// between schoolbook and NTT-based algorithms according to the operand lengths and what
// the field supports. Thus, a caller does not choose an algorithm, only a ring.
//
// A ring caches NTT twiddle factors, so reusing one across many operations over the same
// field is cheaper than building a ring per call.
//
// A *PolyRing is safe for concurrent use.
type PolyRing struct {
	f            Field
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

// NewPolyRing returns a ring performing polynomial arithmetic over f.
//
// It panics if f is nil.
func NewPolyRing(f Field) *PolyRing {
	if f == nil {
		panic("NewPolyRing: nil field")
	}

	return &PolyRing{
		f:            f,
		twiddleCache: map[int]*twiddleSet{},
	}
}

// GetField returns the coefficient field this ring operates over.
func (r *PolyRing) GetField() Field { return r.f }

// NewPolynomial builds a polynomial over this ring's field from coefficients ordered
// from lowest to highest degree. An empty inner yields the zero polynomial.
//
// The polynomial takes ownership of inner rather than copying it.
func (r *PolyRing) NewPolynomial(inner []uint64, isPointRepresentation bool) *Polynomial {
	// An empty sum of terms is zero.
	if len(inner) == 0 {
		inner = []uint64{0}
	}

	return &Polynomial{
		inner: inner,
		isNTT: isPointRepresentation,
		f:     r.f,
	}
}

// Product multiplies every polynomial in polys using a divide-and-conquer product tree,
// in O(n log^2 n) with NTT-based multiplication rather than the O(n^2) of a running
// product.
//
// An empty slice yields the constant polynomial p(x) = 1
//
// polys is only read; the result never aliases any of its elements.
func (r *PolyRing) Product(polys []*Polynomial) *Polynomial {
	switch len(polys) {
	case 0:
		return polyOne(r.f)
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
func (r *PolyRing) productTree(polys []*Polynomial) *Polynomial {
	if len(polys) == 1 {
		return polys[0]
	}

	mid := len(polys) / 2
	left := r.productTree(polys[:mid])
	right := r.productTree(polys[mid:])

	res := r.newDst()
	r.Mul(left, right, res)

	return res
}

// ---------- utilities ----------

// ensureLen resizes c.inner to exactly n, preserving the coefficients already there and
// zeroing whatever the resize exposes.
func ensureLen(c *Polynomial, n int) {
	if cap(c.inner) < n {
		tmp := make([]uint64, n)
		copy(tmp, c.inner)
		c.inner = tmp

		return
	}

	old := len(c.inner)
	c.inner = c.inner[:n]

	if old < n {
		clear(c.inner[old:])
	}
}

// ensureLenCheap resizes c.inner to exactly n without preserving or zeroing anything.
// The caller must write every element before reading it.
func ensureLenCheap(c *Polynomial, n int) {
	if cap(c.inner) < n {
		c.inner = make([]uint64, n)

		return
	}

	c.inner = c.inner[:n]
}

// ensureLenZeroed resizes c.inner to exactly n and zeroes every element.
func ensureLenZeroed(c *Polynomial, n int) {
	c.inner = resizeZeroed(c.inner, n)
}

// resizeZeroed returns buf resized to exactly n elements, all zero, reusing buf's
// array when it is already large enough.
func resizeZeroed(buf []uint64, n int) []uint64 {
	if cap(buf) < n {
		return make([]uint64, n)
	}

	buf = buf[:n]
	clear(buf)

	return buf
}

// newDst returns an empty coefficient-domain polynomial over this ring's field, for
// use as an operation's destination.
func (r *PolyRing) newDst() *Polynomial {
	return &Polynomial{f: r.f, isNTT: false}
}

func (r *PolyRing) trimTrailingZeros(p *Polynomial) {
	if len(p.inner) == 0 || p.isNTT {
		// In NTT domain we keep the fixed size.
		return
	}

	i := len(p.inner) - 1
	for i >= 0 && r.f.Equals(p.inner[i], 0) {
		i--
	}
	p.inner = p.inner[:i+1]
}

// ---------- Poly ops ----------

// Evaluate returns a(x) by Horner's rule.
//
// It panics if a is in the NTT domain: Horner over point values is not the value of the
// polynomial anywhere, and returning that number quietly would be worse than stopping.
func (r *PolyRing) Evaluate(a *Polynomial, x uint64) uint64 {
	if a.isNTT {
		panic("Evaluate expects a coefficient-domain polynomial")
	}

	result := uint64(0)
	fld := r.f

	// horner's rule:
	for i := len(a.inner) - 1; i >= 0; i-- {
		result = fld.Add(a.inner[i], fld.Mul(x, result))
	}

	return result
}

// MulScalar computes c = a * scalar, preserving a's domain. c may alias a.
func (r *PolyRing) MulScalar(a *Polynomial, scalar uint64, c *Polynomial) {
	s := r.f.Reduce(scalar)
	f := r.GetField()

	ensureLen(c, len(a.inner))
	for i := range a.inner {
		c.inner[i] = f.Mul(a.inner[i], s)
	}

	c.f = r.f
	c.isNTT = a.isNTT // scalar mult preserves domain

	r.trimTrailingZeros(c)
}

// Add computes c = a + b. c may alias a or b.
//
// It panics if a and b are over different fields or in different domains.
func (r *PolyRing) Add(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)

	f := r.f

	var av, bv uint64
	for i := 0; i < n; i++ {
		if i < alen {
			av = r.f.Reduce(a.inner[i])
		} else {
			av = 0
		}

		if i < blen {
			bv = r.f.Reduce(b.inner[i])
		} else {
			bv = 0
		}

		c.inner[i] = f.Add(av, bv)
	}

	c.f = r.f
	c.isNTT = a.isNTT
	r.trimTrailingZeros(c)
}

// Sub computes c = a - b. c may alias a or b.
//
// It panics if a and b are over different fields or in different domains.
func (r *PolyRing) Sub(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	c.f = r.f
	c.isNTT = a.isNTT

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)
	minLen := min(alen, blen)

	f := r.f

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
func (r *PolyRing) Mul(a, b, c *Polynomial) {
	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	la, lb := len(a.inner), len(b.inner)
	if la == 0 || lb == 0 {
		c.f, c.inner, c.isNTT = r.f, c.inner[:0], a.isNTT
		return
	}

	// In NTT domain, multiplication is pointwise and preserves NTT representation.
	// (preOpVerification ensures a.isNTT == b.isNTT)
	if a.isNTT {
		r.mulPointwiseInto(a, b, c)

		return
	}

	// Coefficient-domain smart dispatch: schoolbook for small sizes, NTT otherwise.
	if min(la, lb) <= nttMulThreshold || !r.canUseNTTConvolutionLen(la+lb-1) {
		r.mulSchoolbook(a, b, c)
		return
	}

	r.mulViaNTT(a, b, c)
}

func (r *PolyRing) mulSchoolbook(a, b, c *Polynomial) {
	f := r.f

	newLen := len(a.inner) + len(b.inner) - 1

	// Decide where to write: use c.inner if capacity is enough; else allocate.
	var out []uint64
	if c != a && c != b {
		out = resizeZeroed(c.inner, newLen)
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
	c.f = r.f
	c.inner = out
	c.isNTT = false

	r.trimTrailingZeros(c)
}

// subMonomialMul subtracts ai * x^deg * b from rem, in place.
//
// This is one long-division step. Using a specific function rather than Mul+Sub saves
// allocations and a full-length multiplication per step.
//
// The caller guarantees deg(rem) == deg(b)+deg, so rem is long enough to hold every term
// written here. Coefficients of b above its degree are zero and are skipped: they would
// contribute nothing, and indexing them is what would run off the end of rem when b
// carries trailing zeros.
func (r *PolyRing) subMonomialMul(rem *Polynomial, ai uint64, deg int, b *Polynomial) {
	fld := r.f

	for j, bDeg := 0, b.Degree(); j <= bDeg; j++ {
		rem.inner[j+deg] = fld.Sub(rem.inner[j+deg], fld.Mul(ai, b.inner[j]))
	}
}

// Div returns q and rem with a = q*b + rem and deg(rem) < deg(b), choosing between
// schoolbook and NTT-based division by size.
//
// Like Add, Sub and Mul it panics rather than reporting malformed input: a nil
// polynomial, operands over different fields, or an operand in the NTT domain (call
// NttBackward first). Dividing by the zero polynomial panics for the same reason
// integer division by zero does, and for the same reason math/big panics; the
// quotient does not exist, and a caller that can reach a zero divisor is expected to
// say what it means before dividing, not afterwards.
func (r *PolyRing) Div(a, b *Polynomial) (q *Polynomial, rem *Polynomial) {
	if a == nil || b == nil {
		panic("Div: nil polynomial")
	}

	if err := preOpVerification(a, b); err != nil {
		panic(err)
	}

	if a.isNTT || b.isNTT {
		panic("Div expects coefficient-domain polynomials")
	}

	// True degrees is needed. (using len() caused a bug before).
	bDeg := b.Degree()
	aDeg := a.Degree()

	if bDeg < 0 {
		panic("division by the zero polynomial")
	}

	if aDeg < bDeg {
		return polyZero(r.f), a.Copy()
	}

	quotLen := aDeg - bDeg + 1
	// Since divViaNTT uses NTT-based convolution under the hood,
	// we check if the convolution size is supported before deciding which division algorithm to use.
	// 2*quotLen-1 is the convolution size needed
	maxConvLen := max(2*quotLen-1, aDeg+1)

	if quotLen > nttMulThreshold && r.canUseNTTConvolutionLen(maxConvLen) {
		return r.divViaNTT(a, b, aDeg, bDeg)
	}

	return r.divSchoolbook(a, b, aDeg, bDeg)
}

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (r *PolyRing) divSchoolbook(a, b *Polynomial, n, m int) (q *Polynomial, rem *Polynomial) {
	fld := r.f

	u := fld.Inverse(b.LeadCoeff()) // Assumes inverse exists.

	rem = a.Copy()
	qInner := make([]uint64, n-m+1)

	for i := n - m; i >= 0; i-- {
		// TODO: keeping the degree in a variable might save time.
		if rem.Degree() == m+i {
			qInner[i] = fld.Mul(rem.LeadCoeff(), u)
			// the step rem = rem - qInner[i] * x^i * b is done in place,
			// so that we don't allocate a new polynomial for each step.
			r.subMonomialMul(rem, qInner[i], i, b)
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

func (m polyMatrix2x2) Mul(r *PolyRing, other polyMatrix2x2) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: polyAdd(r, polyMul(r, m.a00, other.a00), polyMul(r, m.a01, other.a10)),
		a01: polyAdd(r, polyMul(r, m.a00, other.a01), polyMul(r, m.a01, other.a11)),
		a10: polyAdd(r, polyMul(r, m.a10, other.a00), polyMul(r, m.a11, other.a10)),
		a11: polyAdd(r, polyMul(r, m.a10, other.a01), polyMul(r, m.a11, other.a11)),
	}
}

// mulVecFirst is MulVec's first component on its own: m.a00*a + m.a01*b.
// so it's just one part of the matrix-vector product.
//
// PartialGCD wants only the remainder, not the pair.
func (m polyMatrix2x2) mulVecFirst(r *PolyRing, a, b *Polynomial) *Polynomial {
	return polyAdd(r, polyMul(r, m.a00, a), polyMul(r, m.a01, b))
}

// full matrix vector product: (m.a00*a + m.a01*b, m.a10*a + m.a11*b)
func (m polyMatrix2x2) MulVec(r *PolyRing, a, b *Polynomial) (*Polynomial, *Polynomial) {
	aOut := m.mulVecFirst(r, a, b)
	bOut := polyAdd(r, polyMul(r, m.a10, a), polyMul(r, m.a11, b))
	return aOut, bOut
}

// applyMatrix advances the pair (a, b) by M, trimming both components of the result.
func (r *PolyRing) applyMatrix(M polyMatrix2x2, a, b *Polynomial) (*Polynomial, *Polynomial) {
	aOut, bOut := M.MulVec(r, a, b)

	r.trimTrailingZeros(aOut)
	r.trimTrailingZeros(bOut)

	return aOut, bOut
}

func polyIdentity2x2(f Field) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: polyOne(f),  // x0
		a01: polyZero(f), // y0
		a10: polyZero(f), // x1
		a11: polyOne(f),  // y1
	}
}

// applyStep left-multiplies M by the Euclidean step matrix for quotient q,
//
//	| 0   1 |
//	| 1  -q |
//
// and returns the product of the two:
//
//	(M_{00}, M_{10}) = (M_{10}, M_{00} - q*M_{10})
//	(M_{01}, M_{11}) = (M_{11}, M_{01} - q*M_{11})
//
// writing a specialized function rather than Mul+Sub saves allocations and a full-length multiplication per step.
func (r *PolyRing) applyStep(M polyMatrix2x2, q *Polynomial) polyMatrix2x2 {
	// mulSubInto writes a - q*b into a destination that must alias neither.
	lo := r.newDst()
	hi := r.newDst()

	r.mulSubInto(lo, M.a00, q, M.a10) // lo = M.a00 - q*M.a10
	r.mulSubInto(hi, M.a01, q, M.a11) // hi = M.a01 - q*M.a11

	// building the new matrix
	return polyMatrix2x2{a00: M.a10, a01: M.a11, a10: lo, a11: hi}
}

// rev reverses p into a window of exactly L coefficients: out[i] = p[L-1-i], with any
// index outside p read as zero.
//
// L is the caller's, never derived from p's degree, and that is the point. Reversing a
// whole polynomial means passing deg(p)+1; reversing a truncated series back into a
// polynomial means passing the length it was truncated to. The two differ exactly when
// p's top coefficient is zero -- routine for a quotient series -- and an earlier version
// that found its own anchor by scanning for the degree silently shifted the quotient
// down by one degree for every low-order zero of the dividend.
func (r *PolyRing) rev(p *Polynomial, L int) *Polynomial {
	out := r.newDst()
	if L <= 0 {
		return out
	}

	// specifically not using Copy(): revInPlace below requires that out.inner be exactly L
	// and Copy() would preserve p's length, which may be shorter.
	out.inner = make([]uint64, L)
	for i, v := range p.inner[:min(L, len(p.inner))] {
		out.inner[i] = r.f.Reduce(v)
	}

	revInPlace(out, L)

	return out
}

// revInPlace reverses p's first L coefficients in place, where rev would allocate a
// second buffer to copy them into. Only valid when p already holds exactly L
// coefficients, so that the reversal is a permutation of what is there.
func revInPlace(p *Polynomial, L int) {
	for i := 0; i < L/2; i++ {
		j := L - 1 - i
		p.inner[i], p.inner[j] = p.inner[j], p.inner[i]
	}
}

func nextPow2(n int) int {
	return 1 << (bits.Len(uint(n - 1)))
}

func (r *PolyRing) mulViaNTT(a, b, c *Polynomial) {
	// if both are NTT, pointwise multiply
	if a.isNTT && b.isNTT {
		r.mulPointwiseInto(a, b, c)

		return
	}

	// else, use mulTrunc with total length (coeff-domain out)
	la, lb := len(a.inner), len(b.inner)
	total := la + lb - 1

	// mulTruncInto clears its destination before reading the operands, so it cannot be
	// pointed at one of them. Where c is distinct it writes straight into c's existing
	// array; where c aliases, the result is built to the side first.
	if c != a && c != b {
		r.mulTruncInto(c, a, b, total)

		return
	}

	prod := r.mulTrunc(a, b, total) // NTT under the hood, coeff-domain out
	c.inner = prod.inner
	c.f, c.isNTT = r.f, false
}

// mulPointwiseInto writes the product of a and b into c, all three in the NTT domain,
// where multiplication is pointwise. a and b must carry the same number of evaluations.
func (r *PolyRing) mulPointwiseInto(a, b, c *Polynomial) {
	ensureLen(c, len(a.inner))
	r.pointwiseMult(a, b, c)

	c.f = r.f
	c.isNTT = true
}

func (r *PolyRing) pointwiseMult(a, b, c *Polynomial) {
	f := r.f
	for i := range c.inner {
		c.inner[i] = f.Mul(a.inner[i], b.inner[i])
	}
}

func (r *PolyRing) getNttMulScratch(n int) *nttMulScratch {
	s := nttMulScratchPool.Get().(*nttMulScratch)

	s.a = resizeZeroed(s.a, n)
	s.b = resizeZeroed(s.b, n)

	return s
}

func (r *PolyRing) putNttMulScratch(s *nttMulScratch) {
	nttMulScratchPool.Put(s)
}

// Multiply polynomials and then truncate to the lowest L terms.
// Use NTT under the hood (size = nextPow2(L + L - 1)), then slice [:L].
func (r *PolyRing) mulTrunc(a, b *Polynomial, L int) *Polynomial {
	out := r.newDst()
	r.mulTruncInto(out, a, b, L)

	return out
}

func (r *PolyRing) mulTruncInto(dst *Polynomial, a, b *Polynomial, L int) {
	dst.f = r.f
	dst.isNTT = false

	if L <= 0 {
		dst.inner = dst.inner[:0]
		return
	}

	ensureLenZeroed(dst, L)

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

	aNTT := &Polynomial{f: r.f, inner: scratch.a, isNTT: false}
	bNTT := &Polynomial{f: r.f, inner: scratch.b, isNTT: false}

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
func (r *PolyRing) seriesInverse(b *Polynomial, k int) *Polynomial {
	if k <= 0 {
		return r.newDst()
	}

	if len(b.inner) == 0 || r.f.Equals(b.inner[0], 0) {
		panic("seriesInverse: constant term is zero")
	}

	// Newton doubles m up to k, and mulTruncInto reallocates a destination whose capacity
	// is short. Starting all three at capacity k costs one allocation each instead of one
	// per doubling, and the loop below never grows them again.
	//
	// t and next trade places every iteration, so which of the two is returned depends on
	// the number of steps. tmp never escapes.
	b0 := r.f.Reduce(b.inner[0])
	t := &Polynomial{f: r.f, isNTT: false, inner: make([]uint64, 1, k)}
	t.inner[0] = r.f.Inverse(b0)
	next := &Polynomial{f: r.f, isNTT: false, inner: make([]uint64, 0, k)}
	tmp := &Polynomial{f: r.f, inner: make([]uint64, k)}
	two := r.f.Reduce(2)

	f := r.f
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
func (r *PolyRing) divViaNTT(a, b *Polynomial, n, m int) (q, rem *Polynomial) {
	k := n - m + 1 // quotient length

	// 1) Reverse both inputs whole.
	Astar := r.rev(a, n+1) // n+1 means full reversal of a.
	Bstar := r.rev(b, m+1) // m+1 means full reversal of b.

	// lead(b) maps to Bstar[0]; must be invertible
	if len(Bstar.inner) == 0 || r.f.Equals(Bstar.inner[0], 0) {
		panic("division by polynomial with zero leading coefficient")
	}

	// 2) T = (Bstar)^{-1} mod x^k (Newton series inverse)
	T := r.seriesInverse(Bstar, k) // length k

	// 3) Q* = A* * T mod x^k
	Qstar := r.mulTrunc(Astar, T, k)

	// 4) Reverse Q* back into q, over the quotient length k rather than Q*'s own degree.
	//
	// In the book's notation this step is rev_{k-1}, not rev_k: its subscript is a degree
	// bound, whereas k here is a coefficient count, and deg(q) = n-m = k-1.
	revInPlace(Qstar, k)
	q = Qstar

	// 5) rem = a − q*b
	tmp := r.mulTrunc(q, b, n+1) // full product length (deg = n)
	rem = r.newDst()
	r.Sub(a, tmp, rem)       // coeff-domain subtraction
	r.trimTrailingZeros(rem) // ensure deg(rem) < deg(b)

	return q, rem
}

func (r *PolyRing) canUseNTTConvolutionLen(convLen int) bool {
	if convLen <= 0 {
		return false
	}

	n := nextPow2(convLen)
	modMinusOne := r.f.Modulus() - 1
	if uint64(n) > modMinusOne {
		return false
	}

	return modMinusOne%uint64(n) == 0
}

func polyMul(r *PolyRing, a, b *Polynomial) *Polynomial {
	out := r.newDst()
	r.Mul(a, b, out)
	return out
}

func polyAdd(r *PolyRing, a, b *Polynomial) *Polynomial {
	out := r.newDst()
	r.Add(a, b, out)
	return out
}

// mulSubInto computes dst = a - q*b in one fused pass when q is small,
// avoiding a temporary allocation for the intermediate product q*b.
// When q is large, it falls back to mulFull + Sub.
// All polynomials must be in coefficient domain.
func (r *PolyRing) mulSubInto(dst, a, q, b *Polynomial) {
	lq := len(q.inner)
	lb := len(b.inner)

	// If q or b is empty, dst = a.
	if lq == 0 || lb == 0 {
		dst.f = r.f
		dst.isNTT = false
		ensureLen(dst, len(a.inner))
		copy(dst.inner, a.inner)

		return
	}

	// For large q, fall back to mulFull + Sub with a temporary.
	if min(lq, lb) > nttMulThreshold {
		tmp := r.newDst()
		r.Mul(q, b, tmp)
		r.Sub(a, tmp, dst)

		return
	}

	// Fused schoolbook: dst[k] = a[k] - sum_{i+j=k} q[i]*b[j]
	// Product q*b has length lq + lb - 1.
	prodLen := lq + lb - 1
	n := max(len(a.inner), prodLen)

	f := r.f
	dst.f = f
	dst.isNTT = false
	ensureLenCheap(dst, n)

	// Start with a copy of `a`, zero-extended.
	la := len(a.inner)
	if la > 0 {
		copy(dst.inner[:la], a.inner)
	}
	clear(dst.inner[la:])

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

/*
Half-GCD (HGCD) Implementation and Fast Partial GCD

The following functions implement the O(n log^2 n) version of the Euclidean Algorithm.
Traditional Euclidean Algorithm computes remainders one by one, taking O(n^2) time.
HGCD uses a divide-and-conquer approach by only looking at the most significant
coefficients (the "high parts") to compute transition matrices.

Terminology:
- Euclidean sequence: The sequence of remainders r_0, r_1, ... where r_0=a, r_1=b.
- Transition Matrix M: A 2x2 matrix such that [r_i, r_{i+1}]^T = M * [a, b]^T.
- The pair: the two adjacent terms (r_i, r_{i+1}) the algorithm currently holds, which
  is the whole of its state -- one step advances it to (r_{i+1}, r_{i+2}), and the
  transition matrix records how far it has come. Degrees fall strictly along the
  sequence, so the first entry (the "leading" one) always outranks the second.
*/

// hgcdThreshold is the degree at which the half-GCD recursion bottoms out into the
// iterative routine, below which the recursion overhead exceeds its benefit.
//
// 256 chosen empirically: measured best or tied-best against 128, 512, 1024 and 4096
// at n = 2048, 8192 and 32768, worth 2-8% over the previous 128.
const hgcdThreshold = 256

// PartialGCD runs the extended Euclidean algorithm on (a, b) and stops early, at the
// first remainder in the Euclidean sequence whose degree is strictly less than
// stopDegree. It returns that remainder along with the Bezout coefficients producing
// it: gcd = x*a + y*b. Pass stopDegree = 0 for an ordinary GCD.
//
// Stopping early is what a Reed-Solomon decoder wants. Gao's algorithm halts the
// sequence at (n+k)/2 and reads the message straight out of the pair it stops on,
// rather than running to the true GCD and working backwards.
//
// a and b are only read; both are copied, and an operand in the NTT domain is
// transformed back before use.
//
// Which algorithm runs is not the caller's choice. Below hgcdThreshold the classical
// quadratic sequence is faster outright, and the half-GCD recursion -- O(n log^2 n),
// and an order of magnitude ahead by n = 32768 -- only pays above it. PartialGCD picks
// by degree, so the crossover stays a number this package can re-measure rather than a
// decision frozen into the API.
func (r *PolyRing) PartialGCD(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// We work on copies to preserve the original inputs.
	A := a.Copy()
	B := b.Copy()

	r.ensureNotNttForm(A)
	r.ensureNotNttForm(B)

	// fastGCDMatrix is the driver that uses HGCD to 'jump' through the sequence.
	M := r.fastGCDMatrix(A, B, stopDegree)

	// Apply the final transition matrix to the original inputs to get the desired remainder.
	AOut := M.mulVecFirst(r, A, B)
	r.trimTrailingZeros(AOut)

	// Return the remainder and its corresponding Bézout coefficients for 'a' and 'b':
	// AOut = M.a00 * a + M.a01 * b
	return AOut, M.a00, M.a01
}

func (r *PolyRing) ensureNotNttForm(A *Polynomial) {
	if A.isNTT {
		// Ignoring this error would clear isNTT below on a polynomial still holding
		// evaluations, silently mislabelling it as coefficients.
		if err := r.NttBackward(A); err != nil {
			panic("ensureNotNttForm: " + err.Error())
		}
	}

	A.isNTT = false
}

// fastGCDMatrix is the driver for Fast GCD. It bridges the gap between the starting
// degrees and the target stopDegree using HGCD for the large step.

// Parameters:
// - a, b: Current polynomials in the sequence.
// - stopDegree: The degree boundary we are aiming to cross.
func (r *PolyRing) fastGCDMatrix(a, b *Polynomial, stopDegree int) polyMatrix2x2 {
	// Base Case 1: Target reached.
	aDeg := a.Degree()
	if aDeg < stopDegree || b.Degree() < 0 {
		return polyIdentity2x2(r.f)
	}

	// Base Case 2: Small polynomials, use iterative O(n^2) logic.
	if aDeg < hgcdThreshold {
		// iterative GCD, stops when the degree of the first remainder is below stopDegree.
		return r.iterativePartialExtendedEuclideanMatrix(a, b, stopDegree, stopAtFirstBelow)
	}

	// the distance hgcd must cover to bring the sequence down to stopDegree.
	reduceBy := aDeg - stopDegree

	// 1. Half-GCD Step:
	// Use HGCD to compute a matrix M that reduces degrees significantly.
	M := r.hgcd(a, b, reduceBy)
	aCur, bCur := r.applyMatrix(M, a, b)

	// hgcd stopped because it couldn't guarantee the next step would not cross stopDegree,
	// so we revert to the iterative routine to complete the task.
	return r.iterativePartialExtendedEuclideanMatrix(aCur, bCur, stopDegree, stopAtFirstBelow).Mul(r, M)
}

/*
The half-GCD (HGCD) algorithm is a divide-and-conquer method that computes a GCD transition matrix.
HGCD returns a matrix M such that for (a', b') = M * (a, b),
the degree of b' is reduced by at least 'reduceBy' relative to the original degree of a.
That is: deg(b') < deg(a) - reduceBy.

The HGCD algorithm relies on an insight (proven in `Modern Computer Algebra` by
Joachim von zur Gathen and Jürgen Gerhard) that the high-order coefficients of
the polynomials share the same transition matrix as the full polynomials, allowing us
to use only the top coefficients of a and b to compute a transition matrix that will also reduce
the full polynomials by a significant amount.
Namely, it works by recursively applying itself to the "high parts" of the polynomials,
which are obtained by shifting the inputs to focus on the top coefficients. and thus performing
smaller multiplications and divisions as much as possible.

Procedure:
  - First Call: performs HGCD on the high parts of a and b, covering about half of reduceBy.
  - Bridge: It performs one division to ensure progress.
  - Second Call: It calls itself recursively on the new pair (b,remainder) to advance m further.

The result
This implements the Schönhage strategy of high-part recursion.
*/
func (r *PolyRing) hgcd(a, b *Polynomial, reduceBy int) polyMatrix2x2 {
	degA := a.Degree()
	targetDeg := degA - reduceBy

	// Base Case: target reduction reached or degree too small.
	if b.Degree() < targetDeg || degA < hgcdThreshold {
		// hgcd needs a matrix that gets us as close as possible to the target reduction,
		// so we want it to return before it crosses the target (otherwise, the high-part recursion
		return r.iterativePartialExtendedEuclideanMatrix(a, b, targetDeg, stopBeforeCrossing)
	}

	/*
		Divide and Conquer:
		To cover reduceBy, we first cover half of it.
		We 'shift' the polynomials by firstShift to only look at the top bits.
		Shifting ensures that the recursive calls work on smaller polynomials
		(O(reduceBy) degrees) while the results remain valid for the high-order
		coefficients of the full inputs.
	*/
	halfReduce := reduceBy / 2

	// firstShift discards the low coefficients, retaining 2*halfReduce+1 of them -- the
	// window the recursion is allowed to spend covering halfReduce.
	firstShift := max(degA-2*halfReduce, 0)

	// 1. First Recursive Call (on high parts):
	// Covers halfReduce of the distance.
	R := r.hgcd(r.shiftRight(a, firstShift), r.shiftRight(b, firstShift), halfReduce)

	// 2. Apply the transition matrix R to the full inputs.
	aCur, bCur := r.applyMatrix(R, a, b)

	// reachedDeg is how far down the sequence has come, so degA-reachedDeg is the
	// distance covered so far.
	reachedDeg := bCur.Degree()
	// Check if the reduction target has already been met.
	if reachedDeg < targetDeg || reachedDeg < 0 {
		return R
	}

	// 3. Standard Euclidean Step (Mandatory Progress):
	// Perform one division step to ensure the next recursive call makes progress.
	q, rem := r.Div(aCur, bCur)
	M := r.applyStep(R, q)

	remDeg := rem.Degree()
	// Check if the reduction target is met after division.
	if remDeg < targetDeg || remDeg < 0 {
		return M
	}

	// 4. Second Recursive Call:
	// Whatever distance the first call and the division left uncovered.
	restReduce := reduceBy - (degA - reachedDeg)
	if restReduce <= 0 {
		return M
	}

	// Recalculate the shift against the degree the pair now leads with.
	secondShift := max(reachedDeg-2*restReduce, 0)

	S := r.hgcd(r.shiftRight(bCur, secondShift), r.shiftRight(rem, secondShift), restReduce)

	return S.Mul(r, M)
}

// A gcdStop says where the Euclidean loop leaves the pair, relative to stopDegree.
type gcdStop int

const (
	// stopAtFirstBelow runs until the leading entry has dropped below stopDegree. It is
	// what a caller asking for the answer itself wants.
	stopAtFirstBelow gcdStop = iota

	// stopBeforeCrossing stops one step earlier, leaving the leading entry at or above
	// stopDegree.
	//
	// hgcd needs this. It computes its matrices from polynomials shifted down by k and
	// then applies them to the unshifted pair, and that transfer is licensed only while
	// the matrix has advanced the sequence no further than the retained coefficients
	// determine. Past that point the last quotient depends on coefficients the shift
	// discarded.
	//
	// Stopping at the first remainder below stopDegree always overruns that licence, so
	// hgcd stops short of it.
	stopBeforeCrossing
)

// euclidStep advances (A, B) by one division and folds the quotient into M.
func (r *PolyRing) euclidStep(
	A, B *Polynomial, M polyMatrix2x2,
) (*Polynomial, *Polynomial, polyMatrix2x2) {
	q, rem := r.Div(A, B)

	return B, rem, r.applyStep(M, q)
}

// iterativePartialExtendedEuclideanMatrix is the O(n^2) fallback. It computes the
// transition matrix M for the Euclidean sequence of (a, b), stopping on the side of
// stopDegree that stop selects.
//
// It is the classical Euclidean algorithm, differing only in what it keeps: the Bezout
// coefficients are left in the matrix instead of being read out of it, so the levels
// above can compose this one with theirs.
func (r *PolyRing) iterativePartialExtendedEuclideanMatrix(
	a, b *Polynomial, stopDegree int, stop gcdStop,
) polyMatrix2x2 {
	A := a.Copy()
	B := b.Copy()
	M := polyIdentity2x2(r.f)

	for B.Degree() >= stopDegree && B.Degree() >= 0 {
		A, B, M = r.euclidStep(A, B, M)
	}

	// check if we need one more step.
	if stop == stopAtFirstBelow && A.Degree() >= stopDegree && B.Degree() >= 0 {
		_, _, M = r.euclidStep(A, B, M)
	}

	return M
}

// shiftRight extracts the high coefficients of a polynomial by shifting.
func (r *PolyRing) shiftRight(p *Polynomial, m int) *Polynomial {
	if m <= 0 {
		return p.Copy()
	}
	if m >= len(p.inner) {
		return polyZero(r.f)
	}
	return r.NewPolynomial(p.inner[m:], false)
}
