// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMonomialQuickDiv(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewPolyRing(f)
	t.Run("simple", func(t *testing.T) {
		m1 := newPolynomial(f, []uint64{5, 1}, false)
		m2 := newPolynomial(f, []uint64{3, 1}, false)

		m := &Polynomial{}
		pr.Mul(m1, m2, m)

		q, r := pr.Div(m, m1)

		a.True(r.IsZero())
		a.Equal(m2.ToSlice(), q.ToSlice())

		intr := NewInterpolator(pr)

		q_ := intr.mDivMi(m, m1)
		a.Equal(q.ToSlice(), q_.ToSlice())

		q, r = pr.Div(m, m2)
		a.True(r.IsZero())
		a.Equal(m1.ToSlice(), q.ToSlice())

		q_ = intr.mDivMi(m, m2)
		a.Equal(q.ToSlice(), q_.ToSlice())
	})

	t.Run("complex", func(t *testing.T) {
		xs := []uint64{1, 2, 3, 5, 6, 7}

		intr := NewInterpolator(pr)

		miSlice := intr.createMiSlice(xs)
		m := pr.Product(miSlice)

		for _, mi := range miSlice {
			qQuickDiv := intr.mDivMi(m, mi)
			qLongdiv, _ := pr.Div(m, mi)
			a.Equal(qQuickDiv.ToSlice(), qLongdiv.ToSlice())
		}
	})
}

func TestInterpolation(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewPolyRing(f)

	coeffs := []uint64{0, 1, 2}
	p := newPolynomial(f, coeffs, false)

	intr := NewInterpolator(pr)

	xs, ys := evalPolyForTest(pr, p, 0, 3)

	interpolated, err := intr.Interpolate(xs, ys)
	a.NoError(err)

	a.Equal(p.ToSlice(), interpolated.ToSlice())
}

func FuzzInterpolation(f *testing.F) {
	testcases := []uint64{1, 5, 1 << 62, (1 << 63) - 1}
	for _, tc := range testcases {
		f.Add(tc) // Use f.Add to provide a seed corpus
	}

	fld, err := NewPrimeField(largePrime)
	if err != nil {
		f.FailNow()
	}

	pr := NewPolyRing(fld)

	f.Fuzz(func(t *testing.T, randomSeed uint64) {
		a := assert.New(t)
		const boundingDegree = 10

		p := randomPolynomial(fld, randomSeed, boundingDegree)

		// interpolate a random polynomial
		intr := NewInterpolator(pr)

		xs, ys := evalPolyForTest(pr, p, int(randomSeed), boundingDegree)
		q, err := intr.Interpolate(xs, ys)
		a.NoError(err)

		fmt.Println()
		a.Equal(p.ToSlice(), q.ToSlice())
	})

}

func evalPolyForTest(pr *PolyRing, p *Polynomial, randomSeed, numEvals int) ([]uint64, []uint64) {
	xs := make([]uint64, numEvals)
	for i := range xs {
		xs[i] = p.f.Reduce(uint64(randomSeed + i + 1))
	}

	ys := make([]uint64, len(xs))

	for i, x := range xs {
		ys[i] = pr.Evaluate(p, x)
	}

	return xs, ys
}

func BenchmarkMDivMi(b *testing.B) {
	a := assert.New(b)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewPolyRing(f)

	xs := []uint64{1, 2, 3, 5, 6, 7}

	intr := NewInterpolator(pr)

	miSlice := intr.createMiSlice(xs)
	m := pr.Product(miSlice)

	mi := miSlice[0]

	b.Run("mDivMi", func(b *testing.B) {
		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			intr.mDivMi(m, mi)
		}
	})

	b.Run("LongDiv", func(b *testing.B) {
		b.ResetTimer()

		for i := 0; i < b.N; i++ {
			pr.Div(m, mi)
		}
	})
}

// used to compare the performance of the O(n log^2 n) tree-based product vs the O(n^2) simple product.
func simplePolyProduct(pr *PolyRing, miSlice []*Polynomial) *Polynomial {
	m := makeConstantPoly(pr.GetField(), 1)
	for _, mi := range miSlice {
		pr.Mul(m, mi, m)
	}

	return m
}

func BenchmarkPolyProductComparison(b *testing.B) {
	f, err := NewPrimeField(65537)
	if err != nil {
		b.Fatalf("failed to create field: %v", err)
	}

	pr := NewPolyRing(f)
	intr := NewInterpolator(pr)

	cases := []int{8, 32, 128, 512, 2048}
	for _, n := range cases {
		xs := make([]uint64, n)
		for i := range xs {
			xs[i] = uint64(i + 1)
		}

		miSlice := intr.createMiSlice(xs)

		// Sanity check once per case so benchmarked runs only measure performance.
		tree := pr.Product(miSlice)
		linear := simplePolyProduct(pr, miSlice)
		if !tree.Equals(linear) {
			b.Fatalf("mismatch for n=%d", n)
		}

		b.Run(fmt.Sprintf("n=%d/tree", n), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = pr.Product(miSlice)
			}
		})

		b.Run(fmt.Sprintf("n=%d/simple", n), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = simplePolyProduct(pr, miSlice)
			}
		})
	}
}
