// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"math/rand"
	"testing"
)

func randomPolyWithDegree(f Field, degree int, rng *rand.Rand) *Polynomial {
	coeffs := make([]uint64, degree+1)
	for i := range coeffs {
		coeffs[i] = f.Reduce(uint64(rng.Uint32()))
	}
	if coeffs[0] == 0 {
		coeffs[0] = 1
	}
	if coeffs[degree] == 0 {
		coeffs[degree] = 1
	}
	return newPolynomial(f, coeffs, false)
}

// Baseline copy of old mulTrunc implementation before optimization.
func oldMulTrunc(r *PolyRing, a, b *Polynomial, L int) *Polynomial {
	out := &Polynomial{f: r.f, isNTT: false}
	if L <= 0 {
		return out
	}
	if a == nil || b == nil {
		out.inner = make([]uint64, 1)
		return out
	}

	la := min(len(a.inner), L)
	lb := min(len(b.inner), L)
	if la == 0 || lb == 0 {
		return out
	}

	total := la + lb - 1
	convLen := min(L, total)
	n := nextPow2(total)

	aNTT := &Polynomial{f: r.f, inner: make([]uint64, n), isNTT: false}
	copy(aNTT.inner, a.inner[:la])

	bNTT := &Polynomial{f: r.f, inner: make([]uint64, n), isNTT: false}
	copy(bNTT.inner, b.inner[:lb])

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

	out.inner = aNTT.inner[:convLen]
	return out
}

func oldSeriesInverse(r *PolyRing, b *Polynomial, k int) *Polynomial {
	if k <= 0 {
		return &Polynomial{f: r.f, isNTT: false}
	}
	if len(b.inner) == 0 || r.f.Equals(b.inner[0], 0) {
		panic("seriesInverse: constant term is zero")
	}

	b0 := r.f.Reduce(b.inner[0])
	t := &Polynomial{f: r.f, isNTT: false, inner: []uint64{r.f.Inverse(b0)}}
	two := r.f.Reduce(2)

	f := r.f
	for l := 1; l < k; {
		m := l << 1
		if m > k {
			m = k
		}

		tmp := oldMulTrunc(r, b, t, m)

		if len(tmp.inner) < m {
			z := make([]uint64, m)
			copy(z, tmp.inner)
			tmp.inner = z
		}
		tmp.inner[0] = f.Sub(two, tmp.inner[0])
		for i := 1; i < m; i++ {
			tmp.inner[i] = f.Neg(tmp.inner[i])
		}

		t = oldMulTrunc(r, t, tmp, m)
		l = m
	}
	return t
}

func oldDivNTT(r *PolyRing, a, b *Polynomial) (q, rem *Polynomial) {
	if a == nil || b == nil || a.isNTT || b.isNTT {
		panic("LongDivNTT expects non-nil coefficient-domain polynomials")
	}
	n := len(a.inner) - 1
	m := len(b.inner) - 1
	if m < 0 {
		panic("division by zero polynomial")
	}
	if n < m {
		return &Polynomial{f: r.f, isNTT: false, inner: []uint64{0}}, a.Copy()
	}

	k := n - m + 1
	Astar := r.rev(a, n+1)
	Bstar := r.rev(b, m+1)

	if len(Bstar.inner) == 0 || r.f.Equals(Bstar.inner[0], 0) {
		panic("division by polynomial with zero leading coefficient")
	}

	T := oldSeriesInverse(r, Bstar, k)
	Qstar := oldMulTrunc(r, Astar, T, k)
	// Anchored like divViaNTT's, so this benchmarks the old algorithm rather than the
	// quotient-shift bug both copies of it used to share.
	q = r.rev(Qstar, k)

	prod := oldMulTrunc(r, q, b, n+1)
	rem = &Polynomial{f: r.f, isNTT: false}
	r.Sub(a, prod, rem)
	r.trimTrailingZeros(rem)

	return q, rem
}

func oldNttPartialExtendedEuclidean(r *PolyRing, a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	A := a.Copy()
	B := b.Copy()
	A.isNTT, B.isNTT = false, false

	x0 := makeConstantPoly(r.f, 1)
	x1 := makeConstantPoly(r.f, 0)
	y0 := makeConstantPoly(r.f, 0)
	y1 := makeConstantPoly(r.f, 1)

	tmp1 := &Polynomial{f: r.f}
	tmp2 := &Polynomial{f: r.f}

	for A.Degree() >= stopDegree {
		if B.Degree() < 0 || len(B.inner) == 0 {
			break
		}

		var q, rrem *Polynomial
		if len(A.inner)+len(B.inner) >= nttMulThreshold {
			q, rrem = oldDivNTT(r, A, B)
		} else {
			q, rrem = r.Div(A, B)
		}
		A, B = B, rrem

		r.Mul(q, x1, tmp1)
		r.Sub(x0, tmp1, tmp2)
		x0, x1, tmp2 = x1, tmp2, x0

		r.Mul(q, y1, tmp1)
		r.Sub(y0, tmp1, tmp2)
		y0, y1, tmp2 = y1, tmp2, y0
	}

	return A, x0, y0
}

func benchmarkPolys(tb testing.TB, dividendDegree, divisorDegree int) (*PolyRing, *Polynomial, *Polynomial) {
	tb.Helper()

	fld, err := NewPrimeField(65537)
	if err != nil {
		tb.Fatal(err)
	}
	pr := NewPolyRing(fld)

	rng := rand.New(rand.NewSource(1337))
	a := randomPolyWithDegree(fld, dividendDegree, rng)
	b := randomPolyWithDegree(fld, divisorDegree, rng)

	return pr, a, b
}

func BenchmarkDivNTTLarge_OldVsOptimized(b *testing.B) {
	cases := []struct {
		dividend int
		divisor  int
	}{
		{16384, 8192},
		{32768, 16384},
	}

	for _, tc := range cases {
		name := fmt.Sprintf("degA=%d/degB=%d", tc.dividend, tc.divisor)
		b.Run(name+"/old", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.dividend, tc.divisor)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _ = oldDivNTT(pr, a, d)
			}
		})

		b.Run(name+"/optimized", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.dividend, tc.divisor)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _ = pr.divViaNTT(a, d)
			}
		})
	}
}

func BenchmarkNttPartialEEALarge_OldVsOptimized(b *testing.B) {
	cases := []struct {
		dividend   int
		divisor    int
		stopDegree int
	}{
		{16384, 8192, 8192},
	}

	for _, tc := range cases {
		name := fmt.Sprintf("degA=%d/degB=%d/stop=%d", tc.dividend, tc.divisor, tc.stopDegree)
		b.Run(name+"/old", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.dividend, tc.divisor)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _, _ = oldNttPartialExtendedEuclidean(pr, a, d, tc.stopDegree)
			}
		})

		b.Run(name+"/optimized", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.dividend, tc.divisor)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _, _ = pr.nttPartialExtendedEuclidean(a, d, tc.stopDegree)
			}
		})
	}
}

func BenchmarkMulTruncLarge_OldVsOptimized(b *testing.B) {
	cases := []struct {
		degA int
		degB int
		L    int
	}{
		{16384, 8192, 8193},
		{16384, 16384, 16385},
	}

	for _, tc := range cases {
		name := fmt.Sprintf("degA=%d/degB=%d/L=%d", tc.degA, tc.degB, tc.L)
		b.Run(name+"/old", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.degA, tc.degB)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = oldMulTrunc(pr, a, d, tc.L)
			}
		})

		b.Run(name+"/optimized", func(b *testing.B) {
			pr, a, d := benchmarkPolys(b, tc.degA, tc.degB)
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = pr.mulTrunc(a, d, tc.L)
			}
		})
	}
}
