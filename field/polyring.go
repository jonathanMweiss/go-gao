package field

import (
	"math/bits"
	"sync"
)

type PolyRing interface {
	GetField() Field

	Evaluate(a *Polynomial, x uint64) uint64

	// Assumes polynomial of valid degree.
	NttForward(a *Polynomial) error
	NttBackward(a *Polynomial) error

	// compute c = a * scalar
	MulScalar(a *Polynomial, scalar uint64, c *Polynomial)

	// compute c = a * b
	Mul(a, b, c *Polynomial)
	MulNTT(a, b, c *Polynomial) // multiply in NTT domain, pointwise

	// compute c = a + b
	Add(a, b, c *Polynomial)
	// compute c = a - b
	Sub(a, b, c *Polynomial)

	// Creates quotient and remainder
	Div(a, b *Polynomial) (q *Polynomial, r *Polynomial) // returns quotient, remainder
	DivNTT(a, b *Polynomial) (q, r *Polynomial)          // returns quotient, remainder

	// Extended Euclidean algorithm.
	PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)
	NttPartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)

	PartialGCD(a, b *Polynomial, stopDegree int) (gcd, y *Polynomial)
	NttPartialGCD(a, b *Polynomial, stopDegree int) (gcd, y *Polynomial)
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
func NewDensePolyRing(f Field) PolyRing {
	return &DensePolyRing{
		Field:        f,
		twiddleCache: map[int]*twiddleSet{},
	}
}

func (r *DensePolyRing) GetField() Field { return r.Field }

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
func (r *DensePolyRing) Evaluate(a *Polynomial, x uint64) uint64 {
	if a.isNTT {
		panic("Evaluate not supported in NTT domain")
	}

	result := uint64(0)
	fld := r.Field

	// horner's rule:
	for i := len(a.inner) - 1; i >= 0; i-- {
		result = fld.Add(a.inner[i], fld.Mul(x, result))
	}

	return result
}

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

func (r *DensePolyRing) Add(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
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

func (r *DensePolyRing) Sub(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
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

func (r *DensePolyRing) Mul(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
	}

	if a.isNTT || b.isNTT {
		panic("use MulNTT for pointwise multiplication")
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

	return NewPolynomial(fld, prod, p.isNTT)
}

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (r *DensePolyRing) Div(a, b *Polynomial) (q *Polynomial, rem *Polynomial) {
	if !preOpVerification(a, b) {
		return nil, nil
	}
	fld := r.Field

	if b.isNTT {
		return nil, nil
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

	if len(qInner) == 0 {
		qInner = []uint64{0}
	}

	q = NewPolynomial(fld, qInner, false)
	q.removeLeadingZeroes()

	return q, rem
}

func makeConstantPoly(f Field, u uint64) *Polynomial {
	return NewPolynomial(f, []uint64{u}, false)
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

func polyIdentity2x2(f Field) polyMatrix2x2 {
	return polyMatrix2x2{
		a00: makeConstantPoly(f, 1), // x0
		a01: makeConstantPoly(f, 0), // y0
		a10: makeConstantPoly(f, 0), // x1
		a11: makeConstantPoly(f, 1), // y1
	}
}

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

func (r *DensePolyRing) MulNTT(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
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
func (r *DensePolyRing) DivNTT(a, b *Polynomial) (q, rem *Polynomial) {
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

// mulFull computes c = a*b in coefficient domain, length len(a)+len(b)-1.
// It uses MulNTT when big enough; otherwise falls back to Mul (classic algorithm).
// The decision is based on min(la, lb): when one operand is small, schoolbook
// multiplication is O(min·max) which beats NTT's O(N log N).
func (r *DensePolyRing) mulFull(a, b, c *Polynomial) {
	la, lb := len(a.inner), len(b.inner)
	if la == 0 || lb == 0 {
		c.f, c.inner, c.isNTT = r.Field, c.inner[:0], false
		return
	}
	if min(la, lb) <= nttMulThreshold || !r.canUseNTTConvolutionLen(la+lb-1) {
		r.Mul(a, b, c)
	} else {
		r.MulNTT(a, b, c)
	}
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
	r.mulFull(a, b, out)
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

func (r *DensePolyRing) polyDivSmart(a, b *Polynomial) (q, rem *Polynomial) {
	if b == nil || b.Degree() < 0 {
		panic("division by zero polynomial")
	}

	if a == nil {
		return polyZero(r.Field), polyZero(r.Field)
	}

	if a.Degree() < b.Degree() {
		return polyZero(r.Field), a.Copy()
	}

	quotDeg := a.Degree() - b.Degree()
	// When the quotient degree is small (common in PEEA with errors),
	// schoolbook division is O(quotDeg·degB) which beats DivNTT's O(n log n).
	if quotDeg > nttMulThreshold {
		return r.DivNTT(a, b)
	}

	return r.Div(a, b)
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
		r.mulFull(q, b, tmp)
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
		q, rrem := r.polyDivSmart(A, B)
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

// NttPartialGCD is like NttPartialExtendedEuclidean but only returns (gcd, y),
// skipping computation of the unused x Bézout coefficient for ~2x speedup on the matrix updates.
func (r *DensePolyRing) NttPartialGCD(a, b *Polynomial, stopDegree int) (gcd, y *Polynomial) {
	A := a.Copy()
	B := b.Copy()
	A.isNTT, B.isNTT = false, false

	// Only track right column: (M10, M11).
	M10 := polyZero(r.Field)
	M11 := polyOne(r.Field)

	tmp := &Polynomial{f: r.Field}

	degA := A.Degree()
	degB := B.Degree()

	for degA >= stopDegree {
		if degB < 0 || len(B.inner) == 0 {
			break
		}

		q, rrem := r.polyDivSmart(A, B)
		A, B = B, rrem
		degA, degB = degB, B.Degree()

		// Update only right column: (M10, M11) = (M11, M10 - q*M11)
		r.mulSubInto(tmp, M10, q, M11)
		M10, M11, tmp = M11, tmp, M10
	}

	return A, M10
}

// PartialGCD is like PartialExtendedEuclidean but only returns (gcd, y),
// skipping computation of the unused x Bézout coefficient for ~2x speedup.
func (r *DensePolyRing) PartialGCD(a, b *Polynomial, stopDegree int) (gcd, y *Polynomial) {
	A := a.Copy()
	B := b.Copy()
	degA := A.Degree()
	degB := B.Degree()

	// Only track right column: (M01, M11).
	M01 := polyZero(r.Field)
	M11 := polyOne(r.Field)

	tmp1 := &Polynomial{f: r.Field}
	tmp2 := &Polynomial{f: r.Field}

	for degA >= stopDegree {
		if degB < 0 {
			break
		}

		q, rrem := r.Div(A, B)
		A, B = B, rrem
		degA, degB = degB, B.Degree()

		// (M01, M11) = (M11, M01 - q*M11)
		r.Mul(q, M11, tmp1)
		r.Sub(M01, tmp1, tmp2)
		M01, M11, tmp2 = M11, tmp2, M01
	}

	return A, M01
}
