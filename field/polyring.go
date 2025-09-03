package field

import "sync"

type PolyRing interface {
	Field
	GetField() Field

	Evaluate(a *Polynomial, x uint64) uint64
	// compute c = a * scalar
	MulScalar(a *Polynomial, scalar uint64, c *Polynomial)

	// compute c = a * b
	MulPoly(a, b, c *Polynomial)
	// compute c = a + b
	AddPoly(a, b, c *Polynomial)
	// compute c = a - b
	SubPoly(a, b, c *Polynomial)

	// Creates quotient and remainder
	LongDiv(a, b *Polynomial) (q *Polynomial, r *Polynomial) // returns quotient, remainder
	LongDivNTT(a, b *Polynomial) (q, r *Polynomial)          // returns quotient, remainder

	// Extended Euclidean algorithm.
	PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)
	NttPartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)

	// Assumes it is a polynomial of a valid degree.
	NttForward(a *Polynomial) error
	NttBackward(a *Polynomial) error
}

// DensePolyRing implements PolyRing with optional NTT domain for polynomials.
type DensePolyRing struct {
	Field
	mu           sync.RWMutex
	twiddleCache map[int]*twiddleSet // key: n
}

// NewDensePolyRing constructs a ring over the provided coefficient field.
func NewDensePolyRing(f Field) PolyRing {
	return &DensePolyRing{
		Field:        f,
		mu:           sync.RWMutex{},
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

func (r *DensePolyRing) AddPoly(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
	}

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)

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

		c.inner[i] = r.Add(av, bv)
	}

	c.f = r.Field
	c.isNTT = a.isNTT
	r.trimTrailingZeros(c)
}

func (r *DensePolyRing) SubPoly(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
	}

	alen := len(a.inner)
	blen := len(b.inner)
	n := max(alen, blen)
	ensureLen(c, n)

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

		c.inner[i] = r.Sub(av, bv)
	}

	c.f = r.Field
	c.isNTT = a.isNTT

	r.trimTrailingZeros(c)
}

func (r *DensePolyRing) MulPoly(a, b, c *Polynomial) {
	if !preOpVerification(a, b) {
		panic("preOpVerification failed")
	}

	// Case 1: both inputs are already NTT with same length -> pointwise
	if a.isNTT && b.isNTT {
		n := len(a.inner)
		ensureLen(c, n)
		for i := 0; i < n; i++ {
			c.inner[i] = r.Mul(a.inner[i], b.inner[i])
		}

		c.f = r.Field
		c.isNTT = true

		return
	}

	newLen := len(a.inner) + len(b.inner) - 1

	// Decide where to write: use c.inner if capacity is enough; else allocate.
	var out []uint64
	if cap(c.inner) >= newLen {
		out = c.inner[:newLen]

		for i := range out {
			out[i] = 0
		}

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
			out[i+j] = r.Add(out[i+j], r.Mul(ai, b.inner[j]))
		}
	}

	// Write result into c (safe even if c==a or c==b because we used `out`).
	c.f = a.f
	c.inner = out
	c.isNTT = false

	r.trimTrailingZeros(c)
}

func (r *DensePolyRing) monomialMultPoly(ai uint64, deg int, p *Polynomial) *Polynomial {
	newDegree := len(p.inner) + deg
	fld := r.GetField()
	prod := make([]uint64, newDegree)

	for i := range p.inner {
		prod[i+deg] = fld.Mul(ai, p.inner[i])
	}

	for i := range deg {
		prod[i] = 0
	}

	return NewPolynomial(fld, prod, p.isNTT)
}

// Following Algorithm 2.5 (Polynomial division with remainder) in
// `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard
//
// returns q, r such that p = q*v + r.
func (r *DensePolyRing) LongDiv(a, b *Polynomial) (q *Polynomial, rem *Polynomial) {
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
			r.SubPoly(rem, r.monomialMultPoly(qInner[i], i, b), rem)
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

// returns r= gcd(a,b), x, y such that ax + by = r.
// where r.Degree() < stopDegree.
//
// improved from recursive function using gpt:
func (r *DensePolyRing) PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// Work on local copies ensuring inputs aren't mutated.
	A := a.Copy()
	B := b.Copy()

	// Invariants:
	//   A = x0*a_orig + y0*b_orig
	//   B = x1*a_orig + y1*b_orig
	x0 := makeConstantPoly(r.Field, 1) // 1
	x1 := makeConstantPoly(r.Field, 0) // 0
	y0 := makeConstantPoly(r.Field, 0) // 0
	y1 := makeConstantPoly(r.Field, 1) // 1

	// Reusable temporaries (avoid allocations).
	tmp1 := &Polynomial{f: r.Field} // holds q*x1 or q*y1
	tmp2 := &Polynomial{f: r.Field} // holds x0 - q*x1 or y0 - q*y1

	for A.Degree() >= stopDegree {
		// If B == 0, can't divide further.
		if B.Degree() < 0 {
			break
		}

		// A = q*B + r
		q, rrem := r.LongDiv(A, B)
		A, B = B, rrem // GCD recursive step: gcd(A, B) = gcd(B,rrem)

		// following Bézout's identity:
		// x update: (x0, x1) = (x1, x0 - q*x1)
		r.MulPoly(q, x1, tmp1)    // tmp1 = q * x1
		r.SubPoly(x0, tmp1, tmp2) // tmp2 = x0 - q*x1
		x0, x1, tmp2 = x1, tmp2, x0

		// y update: (y0, y1) = (y1, y0 - q*y1)
		r.MulPoly(q, y1, tmp1)    // tmp1 = q * y1
		r.SubPoly(y0, tmp1, tmp2) // tmp2 = y0 - q*y1
		y0, y1, tmp2 = y1, tmp2, y0
	}

	// gcd = A, x = x0, y = y0
	return A, x0, y0
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

	return &Polynomial{f: f, inner: out, isNTT: false}
}

// NTTDIV: Used GPT instead of implementing.

// Reverse the top L coefficients: rev_L(f) = x^{L-1} * f(1/x) truncated to L.
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
	if n == 0 {
		return 1
	}

	p := 1
	for p < n {
		p <<= 1
	}
	return p
}

// Multiply polynomials and then truncate to the lowest L terms.
// Use NTT under the hood (size = nextPow2(L + L - 1)), then slice [:L].
func (r *DensePolyRing) mulTrunc(a, b *Polynomial, L int) *Polynomial {
	out := &Polynomial{f: r.Field, isNTT: false}
	if L <= 0 {
		return out
	}
	if a == nil || b == nil {
		out.inner = make([]uint64, 1) // 0
		return out
	}

	la := min(len(a.inner), L)
	lb := min(len(b.inner), L)
	if la == 0 || lb == 0 {
		// zero product
		return out
	}

	total := la + lb - 1
	convLen := min(L, total)
	n := nextPow2(total)

	// Prepare coeff-domain buffers of length n
	aNTT := &Polynomial{f: r.Field, inner: make([]uint64, n), isNTT: false}
	for i := 0; i < la; i++ {
		aNTT.inner[i] = r.Reduce(a.inner[i])
	}

	bNTT := &Polynomial{f: r.Field, inner: make([]uint64, n), isNTT: false}
	for i := 0; i < lb; i++ {
		bNTT.inner[i] = r.Reduce(b.inner[i])
	}

	// Forward NTT (these should toggle isNTT to true internally)
	if err := r.NttForward(aNTT); err != nil {
		panic(err)
	}
	if err := r.NttForward(bNTT); err != nil {
		panic(err)
	}

	// Pointwise multiply into aNTT
	for i := 0; i < n; i++ {
		aNTT.inner[i] = r.Mul(aNTT.inner[i], bNTT.inner[i])
	}

	// Inverse NTT back to coeff domain (should toggle isNTT back to false)
	if err := r.nttBackwardNoTrim(aNTT); err != nil {
		panic(err)
	}

	// Truncate to the lowest convLen terms and return in coeff domain
	out.inner = aNTT.inner[:convLen]
	return out
}

// Series inverse modulo x^k using Newton iteration.
// Assumes b[0] != 0; returns t such that (b * t) ≡ 1 mod x^k.
// seriesInverse computes t such that (b * t) ≡ 1 (mod x^k), using Newton iteration.
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
	two := r.Reduce(2)

	for l := 1; l < k; {
		m := l << 1
		if m > k {
			m = k
		}

		// tmp = b*t mod x^m
		tmp := r.mulTrunc(b, t, m)

		// tmp = 2 - tmp (mod x^m)
		if len(tmp.inner) < m {
			z := make([]uint64, m)
			copy(z, tmp.inner)
			tmp.inner = z
		}
		tmp.inner[0] = r.Sub(two, tmp.inner[0])
		for i := 1; i < m; i++ {
			tmp.inner[i] = r.Neg(tmp.inner[i])
		}

		// t = t * tmp mod x^m
		t = r.mulTrunc(t, tmp, m)
		l = m
	}
	return t
}

// LongDivNTT follows `Modern Computer Algebra` by Joachim von zur Gathen and Jürgen Gerhard, section 9.1.
//
// The algorithm applies q' = Rev(a.Copy(),len(a)) * Rev(b.Copy(),len(b)) ^{-1}.
// Then it computes Rev(q',len(a)-len(b)+1) to compute the quotient.
// To compute the remainder rem it follows the relation a = q*b+rem. Namely, rem = a-q*b.
// Rev(*,*) has O(n) complexity. there are finite number of multiplications, and each uses NTT
// thus O(nlogn) complexity.
// The inverse of Rev(b.Copy(),len(b)) is computed via Newton iteration in the method seriesInverse
// with total complexity of O(nlogn).
func (r *DensePolyRing) LongDivNTT(a, b *Polynomial) (q, rem *Polynomial) {
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
	Qstar := r.mulTrunc(Astar, T, k)

	// 4) q = rev_k(Q*)
	q = r.revTop(Qstar, k) // coefficient domain

	// 5) rem = a − q*b
	prod := r.mulTrunc(q, b, n+1) // full product length (deg = n)
	rem = &Polynomial{f: r.Field, isNTT: false}
	r.SubPoly(a, prod, rem)  // coeff-domain subtraction
	r.trimTrailingZeros(rem) // ensure deg(rem) < deg(b)

	return q, rem
}

const nttMulThreshold = 256 // ~coeff count where NTT starts winning

// mulFull computes c = a*b in coefficient domain, length len(a)+len(b)-1.
// It uses mulTrunc with L = total when big enough; otherwise falls back to Mul.
func (r *DensePolyRing) mulFull(a, b, c *Polynomial) {
	la, lb := len(a.inner), len(b.inner)
	if la == 0 || lb == 0 {
		c.f, c.inner, c.isNTT = r.Field, c.inner[:0], false
		return
	}
	total := la + lb - 1
	if total >= nttMulThreshold {
		prod := r.mulTrunc(a, b, total) // NTT under the hood, coeff-domain out
		// write into c without extra allocs when possible
		if cap(c.inner) < total {
			c.inner = make([]uint64, total)
		} else {
			c.inner = c.inner[:total]
		}
		copy(c.inner, prod.inner)
		c.f, c.isNTT = r.Field, false
	} else {
		// naive dense mul
		r.MulPoly(a, b, c)
	}
}

func (r *DensePolyRing) NttPartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// Work on local copies ensuring inputs aren't mutated (coeff domain expected).
	A := a.Copy()
	B := b.Copy()
	A.isNTT, B.isNTT = false, false

	// Invariants:
	//   A = x0*a_orig + y0*b_orig
	//   B = x1*a_orig + y1*b_orig
	x0 := makeConstantPoly(r.Field, 1) // 1
	x1 := makeConstantPoly(r.Field, 0) // 0
	y0 := makeConstantPoly(r.Field, 0) // 0
	y1 := makeConstantPoly(r.Field, 1) // 1

	// Reusable temporaries (avoid allocations).
	tmp1 := &Polynomial{f: r.Field} // holds q*x1 or q*y1
	tmp2 := &Polynomial{f: r.Field} // holds x0 - q*x1 or y0 - q*y1

	for A.Degree() >= stopDegree {
		// If B == 0, can't divide further.
		if B.Degree() < 0 || len(B.inner) == 0 {
			break
		}

		// A = q*B + r  (use NTT-accelerated division when large)
		var q, rrem *Polynomial
		if len(A.inner)+len(B.inner) >= nttMulThreshold { // simple heuristic
			q, rrem = r.LongDivNTT(A, B)
		} else {
			q, rrem = r.LongDiv(A, B)
		}
		A, B = B, rrem // gcd(A,B) = gcd(B,rrem)

		// x update: (x0, x1) = (x1, x0 - q*x1)
		r.mulFull(q, x1, tmp1)    // tmp1 = q * x1   (NTT if big)
		r.SubPoly(x0, tmp1, tmp2) // tmp2 = x0 - q*x1
		x0, x1, tmp2 = x1, tmp2, x0

		// y update: (y0, y1) = (y1, y0 - q*y1)
		r.mulFull(q, y1, tmp1)    // tmp1 = q * y1   (NTT if big)
		r.SubPoly(y0, tmp1, tmp2) // tmp2 = y0 - q*y1
		y0, y1, tmp2 = y1, tmp2, y0
	}

	// gcd = A, x = x0, y = y0
	return A, x0, y0
}
