package field

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

	// Extended Euclidean algorithm.
	PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial)

	// Assumes it is a polynomial of a valid degree.
	NTT(a *Polynomial) *Polynomial
	INTT(a *Polynomial) *Polynomial
}

// DensePolyRing implements PolyRing with optional NTT domain for polynomials.
type DensePolyRing struct {
	Field
}

// NewDensePolyRing constructs a ring over the provided coefficient field.
func NewDensePolyRing(f Field) PolyRing { return &DensePolyRing{Field: f} }

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

	u := fld.Inverse(b.LeadCoeff())

	rem = a.Copy()
	qInner := make([]uint64, n-m+1)

	for i := n - m; i >= 0; i-- {
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
func (r *DensePolyRing) PartialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	if a.Degree() < stopDegree {
		gcd = a.Copy()
		x = makeConstantPoly(r.Field, 1)
		y = makeConstantPoly(r.Field, 0)

		return
	}

	quotient, rem := r.LongDiv(a, b)
	gcd, x1, y1 := r.PartialExtendedEuclidean(b, rem, stopDegree)
	x = y1.Copy()
	y = y1

	r.MulPoly(quotient, y1, y)
	r.SubPoly(x1, y, y)

	return gcd, x, y
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
