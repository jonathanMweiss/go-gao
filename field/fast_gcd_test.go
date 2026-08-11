// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestFastPartialGCDLarge(t *testing.T) {
	a := assert.New(t)
	f, err := NewPrimeField(65537)
	a.NoError(err)

	pr := NewDensePolyRing(f).(*DensePolyRing)

	n := 8192
	k := 4096

	// Mimic Gao decoder stop degree: (n+k)/2 = (8192+4096)/2 = 6144
	stopDegree := (n + k) / 2

	p := randomPolynomial(f, 1337, n)
	q := randomPolynomial(f, 7331, n-1)

	gcd, x, y := pr.FastPartialGCD(p, q, stopDegree)
	a.True(bezoutIdentityHolds(pr, p, q, gcd, x, y), "Bézout identity should hold for Fast GCD")
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

func BenchmarkGCDScaling(b *testing.B) {
	f, err := NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}
	pr := NewDensePolyRing(f).(*DensePolyRing)

	// Test across a range of degrees to see the crossover point and scaling.
	degrees := []int{128, 512, 2048, 8192, 16384}

	for _, n := range degrees {
		// Prepare inputs: a and b such that a target reduction is needed.
		// Gao decoder typical case: n points, k data, stopDegree = (n+k)/2.
		// Let's use k = n/2, so stopDegree = 3n/4.
		stopDegree := (3 * n) / 4
		p := randomPolynomial(f, 42, n)
		q := randomPolynomial(f, 24, n-1)

		b.Run(fmt.Sprintf("Iterative/n=%d", n), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _ = pr.PartialGCD(p, q, stopDegree)
			}
		})

		b.Run(fmt.Sprintf("Fast/n=%d", n), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_, _, _ = pr.FastPartialGCD(p, q, stopDegree)
			}
		})
	}
}

func TestHGCDSecondRecursionCoverage(t *testing.T) {
	a_assert := assert.New(t)
	f, err := NewPrimeField(65537)
	a_assert.NoError(err)
	pr := NewDensePolyRing(f).(*DensePolyRing)

	// Fibonacci polynomials: F_n = x*F_{n-1} + F_{n-2}
	n := 300
	fibs := make([]*Polynomial, n+1)
	fibs[0] = makeConstantPoly(f, 1)
	fibs[1] = NewPolynomial(f, []uint64{0, 1}, false)

	for i := 2; i <= n; i++ {
		xPoly := NewPolynomial(f, []uint64{0, 1}, false)
		prod := &Polynomial{}
		pr.Mul(xPoly, fibs[i-1], prod)
		fibs[i] = &Polynomial{f: f}
		pr.Add(prod, fibs[i-2], fibs[i])
	}

	p := fibs[n]
	q := fibs[n-1]

	stopDegree := 50

	gcdFast, xFast, yFast := pr.FastPartialGCD(p, q, stopDegree)

	// 1. Verify degree condition
	a_assert.Less(gcdFast.Degree(), stopDegree, "Fast GCD degree should be less than stopDegree")
	// 2. Verify Bézout identity: gcd = x*p + y*q
	a_assert.True(bezoutIdentityHolds(pr, p, q, gcdFast, xFast, yFast), "Bézout identity should hold for Fast GCD")

	// 3. Verify that it is a valid remainder in the Euclidean sequence.
	// In the Euclidean sequence of Fibonacci polynomials, every remainder is a Fibonacci polynomial.
	found := false
	for j := 0; j < n; j++ {
		if gcdFast.Equals(fibs[j]) {
			found = true
			break
		}
	}
	a_assert.True(found, "Fast GCD must be a polynomial in the original Fibonacci sequence")
}

func bezoutIdentityHolds(pr *DensePolyRing, a, b, gcd, x, y *Polynomial) bool {
	// ax + by should equal gcd
	ax := polyMul(pr, a, x)
	by := polyMul(pr, b, y)
	return polyAdd(pr, ax, by).Equals(gcd)
}
