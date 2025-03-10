package field

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

const largePrime = 9191248642791733759

func TestCheck(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	slice := f.ElemSlice([]uint64{1, 2, 0, 3})

	fmt.Println(NewPolynomial(slice, false))
}

func TestPolyAdd(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	t.Run("sameSize", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 0, 3})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice, false)

		sum := p1.Add(p2)

		a.Equal([]uint64{2, 4, 0, 6}, sum.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 0, 3})
		slice2 := f.ElemSlice([]uint64{1, 2, 0})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice2, false)

		sum := p1.Add(p2)
		sum2 := p2.Add(p1)
		a.Equal([]uint64{2, 4, 0, 3}, sum.ToSlice())
		a.Equal([]uint64{2, 4, 0, 3}, sum2.ToSlice())
	})

	t.Run("WrapAroundElems", func(t *testing.T) {
		q := f.prime - 1

		slice := f.ElemSlice([]uint64{q, q, q, q})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(f.ElemSlice([]uint64{1, 1, 1, 1}), false)

		sum := p1.Add(p2)
		a.True(sum.IsZero())
	})
}

func TestPolySub(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	t.Run("sameSize", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 0, 3})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice, false)

		sum := p1.Sub(p2)

		a.Equal([]uint64{0}, sum.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 0, 3})
		slice2 := f.ElemSlice([]uint64{1, 2, 0})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice2, false)

		sum := p1.Sub(p2)
		a.Equal([]uint64{0, 0, 0, 3}, sum.ToSlice())

		sum2 := p2.Sub(p1)
		a.Equal([]uint64{0, 0, 0, 154}, sum2.ToSlice())
	})
}

func TestPolyMul(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	t.Run("sameSize", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 3})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice, false)

		prod := p1.Mul(p2)

		a.Equal([]uint64{1, 4, 0, 2, 4}, prod.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 0, 3})
		slice2 := f.ElemSlice([]uint64{1, 2, 0})

		p1 := NewPolynomial(slice, false)
		p2 := NewPolynomial(slice2, false)

		prod := p1.Mul(p2)
		prod2 := p2.Mul(p1)
		a.Equal([]uint64{1, 4, 4, 3, 1}, prod.ToSlice())
		a.True(prod.Equals(prod2))
	})

	t.Run("coeffmod", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 3})

		p1 := NewPolynomial(slice, true)
		p2 := NewPolynomial(slice, true)

		prod := p1.Mul(p2)

		a.Equal([]uint64{1, 4, 4}, prod.ToSlice())
	})
}

func TestPolyLongDiv(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	t.Run("simple", func(t *testing.T) {
		p1 := NewPolynomial(f.ElemSlice([]uint64{1, 2, 3}), false)
		p2 := NewPolynomial(f.ElemSlice([]uint64{1, 2, 3}), false)

		quotient, remainder := p1.LongDiv(p2)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.Equal([]uint64{0}, remainder.ToSlice())

		quotient, remainder = p2.LongDiv(p1)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.Equal([]uint64{0}, remainder.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		p1 := NewPolynomial(f.ElemSlice([]uint64{1, 2, 3}), false)
		p2 := NewPolynomial(f.ElemSlice([]uint64{1, 2}), false)

		quotient, remainder := p1.LongDiv(p2)

		a.Equal([]uint64{4, 4}, quotient.ToSlice())
		a.Equal([]uint64{2}, remainder.ToSlice())

		q, r := p2.LongDiv(p1)
		a.True(p2.Equals(r))
		a.True(q.IsZero())

		p1 = NewPolynomial(f.ElemSlice([]uint64{1, 2, 0, 0, 3}), false)
		p2 = NewPolynomial(f.ElemSlice([]uint64{1, 2}), false)

		quotient, remainder = p1.LongDiv(p2)
		a.Equal([]uint64{3, 1, 3, 4}, quotient.ToSlice())
		a.Equal([]uint64{3}, remainder.ToSlice())
	})

	t.Run("complex", func(t *testing.T) {
		p1 := NewPolynomial(f.ElemSlice([]uint64{1, 0, 0, 0, 2, 3}), false)
		p2 := NewPolynomial(f.ElemSlice([]uint64{1, 0, 1, 0, 2}), false)

		quotient, remainder := p1.LongDiv(p2)

		a.Equal([]uint64{1, 4}, quotient.ToSlice())
		a.Equal([]uint64{0, 1, 4, 1}, remainder.ToSlice())
	})
}

func TestPolyEvaluation(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	t.Run("simple", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{1, 2, 3})

		p := NewPolynomial(slice, false)

		// pairs of {x,p(x)}
		test := [][2]uint64{{1, 1}, {2, 2}, {3, 4}, {4, 2}}
		for _, tt := range test {
			a.Equal(tt[1], p.Eval(tt[0]).Value())
		}
	})

	t.Run("zero", func(t *testing.T) {
		slice := f.ElemSlice([]uint64{0, 0, 0})

		p := NewPolynomial(slice, false)

		// pairs of {x,p(x)}
		test := [][2]uint64{{1, 0}, {2, 0}, {3, 0}, {4, 0}}
		for _, tt := range test {
			a.Equal(tt[1], p.Eval(tt[0]).Value())
		}
	})
}

// Testing the correctness of the partial Extended Euclidean Algorithm
func FuzzPEEA(f *testing.F) {
	testcases := []uint64{1, 5, 1 << 62, (1 << 63) - 1}
	for _, tc := range testcases {
		f.Add(tc) // Use f.Add to provide a seed corpus
	}

	fld, err := NewPrimeField(largePrime)
	if err != nil {
		f.FailNow()
	}

	f.Fuzz(func(t *testing.T, randomSeed uint64) {
		// Create random polynomials.
		maxDegree := 10
		randomPolynomialDegree := randomSeed % (uint64(maxDegree) - 1)

		a := randomPolynomial(fld, randomSeed, maxDegree)
		b := randomPolynomial(fld, randomSeed, int(randomPolynomialDegree))

		for i := 1; i < maxDegree-1; i++ {
			partialDegree := i

			gcd, x, y := PartialExtendedEuclidean(a, b, partialDegree)

			ax := a.Mul(x)
			by := b.Mul(y)
			ax_plus_by := ax.Add(by)

			if !ax_plus_by.Equals(gcd) {
				t.Fatalf("expected %v, got %v", ax_plus_by, gcd)
			}
		}
	})
}

func randomPolynomial(f *PrimeField, seed uint64, maxDegree int) *Polynomial {
	coefficients := make([]Elem, maxDegree)
	for i := 0; i < maxDegree; i++ {
		coefficients[i] = f.ElemFromUint64(seed + uint64(i))
	}

	return NewPolynomial(coefficients, false)
}

func BenchmarkPolyDiv(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}

	p1 := randomPolynomial(f, largePrime/4, 8192)
	p2 := randomPolynomial(f, largePrime/4, 8192/2)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		p1.LongDiv(p2)
	}
}

// TODO: Optimise object creation. We spend a lot of time creating new objects.
func BenchmarkPEEA(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}

	polyMaxDegree := 8193
	p1 := randomPolynomial(f, largePrime/4, polyMaxDegree)   // Large Polynomial.
	p2 := randomPolynomial(f, largePrime/7, polyMaxDegree-1) // The degree of g0 in Gao's decoder, for a polynomial of degree p1.

	for i := 0; i <= 11; i++ {
		b.Run(fmt.Sprintf("partialGCD:remainderDeg<2^%d", i), func(b *testing.B) {
			partialDegree := 1 << i
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				PartialExtendedEuclidean(p1, p2, partialDegree)
			}
		})
	}

	for n := 12; n < 14; n++ {
		k := n / 2

		b.Run(fmt.Sprintf("Partial GCD for %d faults", (n-k)/2), func(b *testing.B) {
			// max errors is (n-k)/2=4096
			p1 := randomPolynomial(f, largePrime/4, n+1) // Large Polynomial.
			p2 := randomPolynomial(f, largePrime/7, n)   // The degree of g0 in Gao's decoder, for a polynomial of degree p1.

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				PartialExtendedEuclidean(p1, p2, (n+k)/2) // Gao's decoder partialEEA.
			}
		})
	}
}
