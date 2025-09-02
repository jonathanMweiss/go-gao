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

	slice := []uint64{1, 2, 0, 3}

	fmt.Println(NewPolynomial(f, slice, false))
}

func TestPolyAdd(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewDensePolyRing(f)
	t.Run("NotInPlace", func(t *testing.T) {

		t.Run("sameSize", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, slice, false)
			sum := &Polynomial{f: f}

			pr.AddPoly(p1, p2, sum)
			a.Equal([]uint64{2, 4, 0, 6}, sum.ToSlice())
		})

		t.Run("differentSizes", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}
			slice2 := []uint64{1, 2, 0}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, slice2, false)
			sum, sum2 := &Polynomial{f: f}, &Polynomial{f: f}

			pr.AddPoly(p1, p2, sum)
			pr.AddPoly(p2, p1, sum2)

			a.Equal([]uint64{2, 4, 0, 3}, sum.ToSlice())
			a.Equal([]uint64{2, 4, 0, 3}, sum2.ToSlice())
		})

		t.Run("WrapAroundElems", func(t *testing.T) {
			q := f.Modulus() - 1

			slice := []uint64{q, q, q, q}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, []uint64{1, 1, 1, 1}, false)
			sum := &Polynomial{f: f}
			pr.AddPoly(p1, p2, sum)
			a.True(sum.IsZero())
		})
	})

	t.Run("InPlace", func(t *testing.T) {

		t.Run("sameSize", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, slice, false)

			p1cpy := p1.Copy()
			p2cpy := p2.Copy()
			pr.AddPoly(p1, p2, p1)
			a.Equal([]uint64{2, 4, 0, 6}, p1.ToSlice())

			p1 = p1cpy.Copy()
			p2 = p2cpy.Copy()

			pr.AddPoly(p1, p2, p2)
			a.Equal([]uint64{2, 4, 0, 6}, p2.ToSlice())
		})

		t.Run("differentSizes", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}
			slice2 := []uint64{1, 2, 0}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, slice2, false)
			p1cpy, p2cpy := p1.Copy(), p2.Copy()
			reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

			pr.AddPoly(p1, p2, p1)
			a.Equal([]uint64{2, 4, 0, 3}, p1.ToSlice())
			reset()

			pr.AddPoly(p2, p1, p1)
			a.Equal([]uint64{2, 4, 0, 3}, p1.ToSlice())
			reset()

			// Now the other side:
			pr.AddPoly(p1, p2, p2)
			a.Equal([]uint64{2, 4, 0, 3}, p2.ToSlice())
			reset()

			pr.AddPoly(p2, p1, p2)
			a.Equal([]uint64{2, 4, 0, 3}, p2.ToSlice())
			reset()
		})

		t.Run("WrapAroundElems", func(t *testing.T) {
			q := f.Modulus() - 1

			slice := []uint64{q, q, q, q}

			p1 := NewPolynomial(f, slice, false)
			p2 := NewPolynomial(f, []uint64{1, 1, 1, 1}, false)
			cpy := p1.Copy()
			pr.AddPoly(cpy, p2, cpy)
			a.True(cpy.IsZero())

			cpy = p2.Copy()
			pr.AddPoly(cpy, p1, cpy)
			a.True(cpy.IsZero())
		})
	})

}

func TestPolySub(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewDensePolyRing(f)

	t.Run("sameSize", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}

		p1 := NewPolynomial(f, slice, false)
		p2 := NewPolynomial(f, slice, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()
		reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

		pr.SubPoly(p1, p2, p1)
		a.True(p1.IsZero())
		reset()

		pr.SubPoly(p1, p2, p2)
		a.True(p2.IsZero())
		reset()
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}
		slice2 := []uint64{1, 2, 0}

		p1 := NewPolynomial(f, slice, false)
		p2 := NewPolynomial(f, slice2, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()
		reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

		pr.SubPoly(p1, p2, p1)
		a.Equal([]uint64{0, 0, 0, 3}, p1.ToSlice())
		reset()

		pr.SubPoly(p1, p2, p2)
		a.Equal([]uint64{0, 0, 0, 3}, p2.ToSlice())
		reset()

		pr.SubPoly(p2, p1, p1)
		a.Equal([]uint64{0, 0, 0, 154}, p1.ToSlice())
		reset()

		pr.SubPoly(p2, p1, p2)
		a.Equal([]uint64{0, 0, 0, 154}, p2.ToSlice())
	})

}

func TestPolyMul(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewDensePolyRing(f)
	t.Run("sameSize", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p1 := NewPolynomial(f, slice, false)
		p2 := NewPolynomial(f, slice, false)

		pr.MulPoly(p1, p2, p1)

		a.Equal([]uint64{1, 4, 0, 2, 4}, p1.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}
		slice2 := []uint64{1, 2, 0}

		p1 := NewPolynomial(f, slice, false)
		p2 := NewPolynomial(f, slice2, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()

		pr.MulPoly(p1, p2, p1)
		pr.MulPoly(p2cpy, p1cpy, p2)
		a.Equal([]uint64{1, 4, 4, 3, 1}, p1.ToSlice())
		a.True(p1.Equals(p2))
	})

	t.Run("inNTT", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p1 := NewPolynomial(f, slice, true)
		p2 := NewPolynomial(f, slice, true)

		pr.MulPoly(p1, p2, p1)
		a.Equal([]uint64{1, 4, 4}, p1.ToSlice())
	})
}

func TestPolyLongDiv(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewDensePolyRing(f)
	t.Run("simple", func(t *testing.T) {
		p1 := NewPolynomial(f, []uint64{1, 2, 3}, false)
		p2 := NewPolynomial(f, []uint64{1, 2, 3}, false)

		quotient, remainder := pr.LongDiv(p1, p2)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.True(remainder.IsZero())

		quotient, remainder = pr.LongDiv(p2, p1)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.True(remainder.IsZero())
	})

	t.Run("differentSizes", func(t *testing.T) {
		p1 := NewPolynomial(f, []uint64{1, 2, 3}, false)
		p2 := NewPolynomial(f, []uint64{1, 2}, false)

		quotient, remainder := pr.LongDiv(p1, p2)

		a.Equal([]uint64{4, 4}, quotient.ToSlice())
		a.Equal([]uint64{2}, remainder.ToSlice())

		q, r := pr.LongDiv(p2, p1)
		a.True(p2.Equals(r))
		a.True(q.IsZero())

		p1 = NewPolynomial(f, []uint64{1, 2, 0, 0, 3}, false)
		p2 = NewPolynomial(f, []uint64{1, 2}, false)

		quotient, remainder = pr.LongDiv(p1, p2)
		a.Equal([]uint64{3, 1, 3, 4}, quotient.ToSlice())
		a.Equal([]uint64{3}, remainder.ToSlice())
	})

	t.Run("complex", func(t *testing.T) {
		p1 := NewPolynomial(f, []uint64{1, 0, 0, 0, 2, 3}, false)
		p2 := NewPolynomial(f, []uint64{1, 0, 1, 0, 2}, false)

		quotient, remainder := pr.LongDiv(p1, p2)

		a.Equal([]uint64{1, 4}, quotient.ToSlice())
		a.Equal([]uint64{0, 1, 4, 1}, remainder.ToSlice())
	})
}

func TestPolyEvaluation(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewDensePolyRing(f)
	t.Run("simple", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p := NewPolynomial(f, slice, false)

		// pairs of {x,p(x)}
		test := [][2]uint64{{1, 1}, {2, 2}, {3, 4}, {4, 2}}
		for _, tt := range test {
			a.Equal(tt[1], pr.Evaluate(p, tt[0]))
		}
	})

	t.Run("zero", func(t *testing.T) {
		slice := []uint64{0, 0, 0}

		p := NewPolynomial(f, slice, false)

		// pairs of {x,p(x)}
		test := [][2]uint64{{1, 0}, {2, 0}, {3, 0}, {4, 0}}
		for _, tt := range test {
			a.Equal(tt[1], pr.Evaluate(p, tt[0]))
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

	pr := NewDensePolyRing(fld)

	f.Fuzz(func(t *testing.T, randomSeed uint64) {
		// Create random polynomials.
		maxDegree := 10
		randomPolynomialDegree := randomSeed % (uint64(maxDegree) - 1)

		a := randomPolynomial(fld, randomSeed, maxDegree)
		b := randomPolynomial(fld, randomSeed, int(randomPolynomialDegree))

		for i := 1; i < maxDegree-1; i++ {
			partialDegree := i

			gcd, x, y := pr.PartialExtendedEuclidean(a, b, partialDegree)

			ax, by, ax_plus_by := &Polynomial{}, &Polynomial{}, &Polynomial{}
			pr.MulPoly(a, x, ax)
			pr.MulPoly(b, y, by)
			pr.AddPoly(ax, by, ax_plus_by)

			if !ax_plus_by.Equals(gcd) {
				t.Fatalf("expected %v, got %v", ax_plus_by, gcd)
			}
		}
	})
}

func randomPolynomial(f Field, seed uint64, maxDegree int) *Polynomial {
	coefficients := make([]uint64, maxDegree)
	for i := 0; i < maxDegree; i++ {
		coefficients[i] = f.Reduce(seed + uint64(i))
	}

	return NewPolynomial(f, coefficients, false)
}

func BenchmarkPolyDiv(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}
	pr := NewDensePolyRing(f)

	p1 := randomPolynomial(f, largePrime/4, 8192)
	p2 := randomPolynomial(f, largePrime/4, 8192/2)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		pr.LongDiv(p1, p2)
	}
}

// TODO: Optimise object creation. We spend a lot of time creating new objects.
func BenchmarkPEEA(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}

	pr := NewDensePolyRing(f)

	polyMaxDegree := 8193
	p1 := randomPolynomial(f, largePrime/4, polyMaxDegree)   // Large Polynomial.
	p2 := randomPolynomial(f, largePrime/7, polyMaxDegree-1) // The degree of g0 in Gao's decoder, for a polynomial of degree p1.

	for i := 0; i <= 11; i++ {
		b.Run(fmt.Sprintf("partialGCD:remainderDeg<2^%d", i), func(b *testing.B) {
			partialDegree := 1 << i
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				pr.PartialExtendedEuclidean(p1, p2, partialDegree)
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
				pr.PartialExtendedEuclidean(p1, p2, (n+k)/2) // Gao's decoder partialEEA.
			}
		})
	}
}

func makeRoots(n int) []uint64 {
	roots := make([]uint64, n)
	for i := 0; i < n; i++ {
		roots[i] = uint64(i*7 + 3)
	}
	return roots
}

func TestLocatorPolynomial(t *testing.T) {
	roots := makeRoots(15)

	f, err := NewPrimeField(largePrime)
	if err != nil {
		t.Fatal(err)
	}
	pr := NewDensePolyRing(f)

	p := PolyProductMonicNegRoots(f, roots)

	intr := NewInterpolator(pr)

	q := PolyProduct(pr, intr.createMiSlice(roots))

	if !p.Equals(q) {
		t.FailNow()
	}
}

var benchPolySink *Polynomial // avoid DCE

/*
pkg: github.com/jonathanmweiss/go-gao/field
BenchmarkPolyProductMonicNegRoots
BenchmarkPolyProductMonicNegRoots/n=15
BenchmarkPolyProductMonicNegRoots/n=15-10         	 1799329	       661.1 ns/op	     304 B/op	       3 allocs/op
BenchmarkPolyProductMonicNegRoots/n=32
BenchmarkPolyProductMonicNegRoots/n=32-10         	  432550	      2724 ns/op	     624 B/op	       3 allocs/op
BenchmarkPolyProductMonicNegRoots/n=64
BenchmarkPolyProductMonicNegRoots/n=64-10         	   95090	     11852 ns/op	    1200 B/op	       3 allocs/op
BenchmarkPolyProductMonicNegRoots/n=128
BenchmarkPolyProductMonicNegRoots/n=128-10        	   14961	     79504 ns/op	    2352 B/op	       3 allocs/op
BenchmarkPolyProductMonicNegRoots/n=256
BenchmarkPolyProductMonicNegRoots/n=256-10        	    2966	    406873 ns/op	    4656 B/op	       3 allocs/op
*/
func BenchmarkPolyProductMonicNegRoots(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.Fatal(err)
	}

	for _, n := range []int{15, 32, 64, 128, 256} {
		roots := makeRoots(n) // prepare inputs outside timed loop

		b.Run(fmt.Sprintf("n=%d", n), func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			var p *Polynomial
			for i := 0; i < b.N; i++ {
				p = PolyProductMonicNegRoots(f, roots)
			}
			b.StopTimer()
			benchPolySink = p
		})
	}
}

func TestDivNTT(t *testing.T) {
	a := assert.New(t)
	f, err := NewPrimeField(65537)
	a.NoError(err)

	for _, maxDegree := range []int{16, 64, 256, 1024} {
		p := randomPolynomial(f, 12345, maxDegree)
		q := randomPolynomial(f, 67890, maxDegree/2)

		pr := NewDensePolyRing(f)
		quo1, rem1 := pr.LongDiv(p.Copy(), q.Copy())
		quo2, rem2 := pr.LongDivNTT(p.Copy(), q.Copy())
		a.True(quo1.Equals(quo2))
		a.True(rem1.Equals(rem2))
	}
}
