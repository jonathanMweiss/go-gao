// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

// A DivisorCache holds the part of a polynomial division that depends only on the
// divisor: the series inverse of the reversed divisor, which [PolyRing.Div] otherwise
// rebuilds by Newton iteration for each dividend. Dividing many polynomials by one
// divisor is then a multiplication apiece.
//
// A cache fixes a largest quotient length. A shorter quotient reuses the same inverse:
// T*rev(b) = 1 mod x^k gives the same mod x^j for every j <= k, because x^j divides x^k.
// DivBy cuts it to the quotient length before multiplying, which is a saving rather than
// a requirement -- the product is truncated there anyway.
//
// A quotient longer than the cache covers falls back to an ordinary division, as does one
// small enough that Div would not have used the transform at all.
//
// The cache copies the divisor, so later changes to the caller's polynomial do not reach
// it. It is read-only once built and safe for concurrent use.
type DivisorCache struct {
	b *Polynomial
	// t is rev(b) inverted as a power series, modulo x^maxQuot.
	t       *Polynomial
	bDeg    int
	maxQuot int
}

// NewDivisorCache prepares b as a divisor for quotients of up to maxQuotLen coefficients.
//
// It panics on a nil, NTT-domain, or zero divisor, matching [PolyRing.Div].
func (r *PolyRing) NewDivisorCache(b *Polynomial, maxQuotLen int) *DivisorCache {
	if b == nil {
		panic("NewDivisorCache: nil polynomial")
	}

	if b.isNTT {
		panic("NewDivisorCache expects a coefficient-domain polynomial")
	}

	bDeg := b.Degree()
	if bDeg < 0 {
		panic("division by the zero polynomial")
	}

	maxQuotLen = max(maxQuotLen, 1)

	rev := r.newDst()
	r.revInto(rev, b, bDeg+1)

	// lead(b) maps to rev[0]; must be invertible.
	if len(rev.inner) == 0 || r.f.Equals(rev.inner[0], 0) {
		panic("division by polynomial with zero leading coefficient")
	}

	// seriesInverse hands back a pooled polynomial, and this one outlives the call.
	scratch := r.seriesInverse(rev, maxQuotLen)
	t := scratch.Copy()
	r.returnPoly(scratch)

	return &DivisorCache{
		b:       b.Copy(),
		t:       t,
		bDeg:    bDeg,
		maxQuot: maxQuotLen,
	}
}

// Divisor returns a copy of the polynomial the cache divides by.
func (dc *DivisorCache) Divisor() *Polynomial { return dc.b.Copy() }

// MaxQuotientLen is the longest quotient the cache serves without falling back.
func (dc *DivisorCache) MaxQuotientLen() int { return dc.maxQuot }

// DivBy returns q and rem with a = q*dc.Divisor() + rem, reusing the cached inverse
// where it applies. The results match [PolyRing.Div] exactly.
func (r *PolyRing) DivBy(a *Polynomial, dc *DivisorCache) (q, rem *Polynomial) {
	if a == nil || dc == nil {
		panic("DivBy: nil argument")
	}

	if err := preOpVerification(a, dc.b); err != nil {
		panic(err)
	}

	if a.isNTT {
		panic("DivBy expects coefficient-domain polynomials")
	}

	aDeg := a.Degree()
	if aDeg < dc.bDeg {
		return polyZero(r.f), a.Copy()
	}

	quotLen := aDeg - dc.bDeg + 1
	// 2*quotLen-1 is the convolution size the reversal-and-multiply below needs.
	maxConvLen := max(2*quotLen-1, aDeg+1)

	// Outside what the cache covers, or below the size where the transform pays for
	// itself, plain division is both correct and faster.
	if quotLen > dc.maxQuot || quotLen <= nttMulThreshold || !r.canUseNTTConvolutionLen(maxConvLen) {
		return r.Div(a, dc.b)
	}

	aStar := r.borrowPoly(0)
	defer r.returnPoly(aStar)

	r.revInto(aStar, a, aDeg+1)

	// Q* = A* * T mod x^quotLen, then reversed back over the quotient length.
	q = r.mulTrunc(aStar, dc.truncated(quotLen), quotLen)
	revInPlace(q, quotLen)

	rem = r.newDst()
	r.mulSubInto(rem, a, q, dc.b)

	return q, rem
}

// truncated views the first k coefficients of the cached inverse. The view shares the
// backing array, which nothing writes to, so concurrent divisions may hold one at once.
func (dc *DivisorCache) truncated(k int) *Polynomial {
	if k >= len(dc.t.inner) {
		return dc.t
	}

	return &Polynomial{f: dc.t.f, inner: dc.t.inner[:k], isNTT: false}
}
