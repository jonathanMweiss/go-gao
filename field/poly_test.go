// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
)

const largePrime = 9191248642791733759

func TestCheck(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	slice := []uint64{1, 2, 0, 3}

	fmt.Println(newPolynomial(f, slice, false))
}

func TestPolyAdd(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewPolyRing(f)
	t.Run("NotInPlace", func(t *testing.T) {

		t.Run("sameSize", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, slice, false)
			sum := &Polynomial{f: f}

			pr.Add(p1, p2, sum)
			a.Equal([]uint64{2, 4, 0, 6}, sum.ToSlice())
		})

		t.Run("differentSizes", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}
			slice2 := []uint64{1, 2, 0}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, slice2, false)
			sum, sum2 := &Polynomial{f: f}, &Polynomial{f: f}

			pr.Add(p1, p2, sum)
			pr.Add(p2, p1, sum2)

			a.Equal([]uint64{2, 4, 0, 3}, sum.ToSlice())
			a.Equal([]uint64{2, 4, 0, 3}, sum2.ToSlice())
		})

		t.Run("WrapAroundElems", func(t *testing.T) {
			q := f.Modulus() - 1

			slice := []uint64{q, q, q, q}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, []uint64{1, 1, 1, 1}, false)
			sum := &Polynomial{f: f}
			pr.Add(p1, p2, sum)
			a.True(sum.IsZero())
		})
	})

	t.Run("InPlace", func(t *testing.T) {

		t.Run("sameSize", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, slice, false)

			p1cpy := p1.Copy()
			p2cpy := p2.Copy()
			pr.Add(p1, p2, p1)
			a.Equal([]uint64{2, 4, 0, 6}, p1.ToSlice())

			p1 = p1cpy.Copy()
			p2 = p2cpy.Copy()

			pr.Add(p1, p2, p2)
			a.Equal([]uint64{2, 4, 0, 6}, p2.ToSlice())
		})

		t.Run("differentSizes", func(t *testing.T) {
			slice := []uint64{1, 2, 0, 3}
			slice2 := []uint64{1, 2, 0}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, slice2, false)
			p1cpy, p2cpy := p1.Copy(), p2.Copy()
			reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

			pr.Add(p1, p2, p1)
			a.Equal([]uint64{2, 4, 0, 3}, p1.ToSlice())
			reset()

			pr.Add(p2, p1, p1)
			a.Equal([]uint64{2, 4, 0, 3}, p1.ToSlice())
			reset()

			// Now the other side:
			pr.Add(p1, p2, p2)
			a.Equal([]uint64{2, 4, 0, 3}, p2.ToSlice())
			reset()

			pr.Add(p2, p1, p2)
			a.Equal([]uint64{2, 4, 0, 3}, p2.ToSlice())
			reset()
		})

		t.Run("WrapAroundElems", func(t *testing.T) {
			q := f.Modulus() - 1

			slice := []uint64{q, q, q, q}

			p1 := newPolynomial(f, slice, false)
			p2 := newPolynomial(f, []uint64{1, 1, 1, 1}, false)
			cpy := p1.Copy()
			pr.Add(cpy, p2, cpy)
			a.True(cpy.IsZero())

			cpy = p2.Copy()
			pr.Add(cpy, p1, cpy)
			a.True(cpy.IsZero())
		})
	})

}

func TestPolySub(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	pr := NewPolyRing(f)

	t.Run("sameSize", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}

		p1 := newPolynomial(f, slice, false)
		p2 := newPolynomial(f, slice, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()
		reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

		pr.Sub(p1, p2, p1)
		a.True(p1.IsZero())
		reset()

		pr.Sub(p1, p2, p2)
		a.True(p2.IsZero())
		reset()
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}
		slice2 := []uint64{1, 2, 0}

		p1 := newPolynomial(f, slice, false)
		p2 := newPolynomial(f, slice2, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()
		reset := func() { p1, p2 = p1cpy.Copy(), p2cpy.Copy() }

		pr.Sub(p1, p2, p1)
		a.Equal([]uint64{0, 0, 0, 3}, p1.ToSlice())
		reset()

		pr.Sub(p1, p2, p2)
		a.Equal([]uint64{0, 0, 0, 3}, p2.ToSlice())
		reset()

		pr.Sub(p2, p1, p1)
		a.Equal([]uint64{0, 0, 0, 154}, p1.ToSlice())
		reset()

		pr.Sub(p2, p1, p2)
		a.Equal([]uint64{0, 0, 0, 154}, p2.ToSlice())
	})

}

func TestPolyMul(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewPolyRing(f)
	t.Run("sameSize", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p1 := newPolynomial(f, slice, false)
		p2 := newPolynomial(f, slice, false)

		pr.Mul(p1, p2, p1)

		a.Equal([]uint64{1, 4, 0, 2, 4}, p1.ToSlice())
	})

	t.Run("differentSizes", func(t *testing.T) {
		slice := []uint64{1, 2, 0, 3}
		slice2 := []uint64{1, 2, 0}

		p1 := newPolynomial(f, slice, false)
		p2 := newPolynomial(f, slice2, false)
		p1cpy, p2cpy := p1.Copy(), p2.Copy()

		pr.Mul(p1, p2, p1)
		pr.Mul(p2cpy, p1cpy, p2)
		a.Equal([]uint64{1, 4, 4, 3, 1}, p1.ToSlice())
		a.True(p1.Equals(p2))
	})

	t.Run("inNTT", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p1 := newPolynomial(f, slice, true)
		p2 := newPolynomial(f, slice, true)

		pr.Mul(p1, p2, p1)
		a.Equal([]uint64{1, 4, 4}, p1.ToSlice())
	})
}

func TestPolyDiv(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewPolyRing(f)
	t.Run("simple", func(t *testing.T) {
		p1 := newPolynomial(f, []uint64{1, 2, 3}, false)
		p2 := newPolynomial(f, []uint64{1, 2, 3}, false)

		quotient, remainder := pr.Div(p1, p2)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.True(remainder.IsZero())

		quotient, remainder = pr.Div(p2, p1)
		a.Equal([]uint64{1}, quotient.ToSlice())
		a.True(remainder.IsZero())
	})

	t.Run("differentSizes", func(t *testing.T) {
		p1 := newPolynomial(f, []uint64{1, 2, 3}, false)
		p2 := newPolynomial(f, []uint64{1, 2}, false)

		quotient, remainder := pr.Div(p1, p2)

		a.Equal([]uint64{4, 4}, quotient.ToSlice())
		a.Equal([]uint64{2}, remainder.ToSlice())

		q, r := pr.Div(p2, p1)
		a.True(p2.Equals(r))
		a.True(q.IsZero())

		p1 = newPolynomial(f, []uint64{1, 2, 0, 0, 3}, false)
		p2 = newPolynomial(f, []uint64{1, 2}, false)

		quotient, remainder = pr.Div(p1, p2)
		a.Equal([]uint64{3, 1, 3, 4}, quotient.ToSlice())
		a.Equal([]uint64{3}, remainder.ToSlice())
	})

	t.Run("complex", func(t *testing.T) {
		p1 := newPolynomial(f, []uint64{1, 0, 0, 0, 2, 3}, false)
		p2 := newPolynomial(f, []uint64{1, 0, 1, 0, 2}, false)

		quotient, remainder := pr.Div(p1, p2)

		a.Equal([]uint64{1, 4}, quotient.ToSlice())
		a.Equal([]uint64{0, 1, 4, 1}, remainder.ToSlice())
	})
}

func TestPolyEvaluation(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(5)
	a.NoError(err)

	pr := NewPolyRing(f)
	t.Run("simple", func(t *testing.T) {
		slice := []uint64{1, 2, 3}

		p := newPolynomial(f, slice, false)

		// pairs of {x,p(x)}
		test := [][2]uint64{{1, 1}, {2, 2}, {3, 4}, {4, 2}}
		for _, tt := range test {
			a.Equal(tt[1], pr.Evaluate(p, tt[0]))
		}
	})

	t.Run("zero", func(t *testing.T) {
		slice := []uint64{0, 0, 0}

		p := newPolynomial(f, slice, false)

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

	pr := NewPolyRing(fld)

	f.Fuzz(func(t *testing.T, randomSeed uint64) {
		// Create random polynomials.
		maxDegree := 10
		randomPolynomialDegree := randomSeed % (uint64(maxDegree) - 1)

		a := randomPolynomial(fld, randomSeed, maxDegree)
		b := randomPolynomial(fld, randomSeed, int(randomPolynomialDegree))

		for i := 1; i < maxDegree-1; i++ {
			partialDegree := i

			gcd, x, y := pr.partialExtendedEuclidean(a, b, partialDegree)

			ax, by, ax_plus_by := &Polynomial{}, &Polynomial{}, &Polynomial{}
			pr.Mul(a, x, ax)
			pr.Mul(b, y, by)
			pr.Add(ax, by, ax_plus_by)

			if !ax_plus_by.Equals(gcd) {
				t.Fatalf("expected %v, got %v", ax_plus_by, gcd)
			}
		}
	})
}

// newPolynomial builds a polynomial directly over f, bypassing PolyRing.
func newPolynomial(f Field, inner []uint64, isPointRepresentation bool) *Polynomial {
	if len(inner) == 0 {
		inner = []uint64{0}
	}

	return &Polynomial{f: f, inner: inner, isNTT: isPointRepresentation}
}

// randomPolynomial builds a pseudo-random polynomial of the given degree, deterministic
// in seed.
//
// It must be genuinely random. It previously used coefficients seed, seed+1, seed+2,
// ... — an arithmetic progression, which (1-x)^2 nearly annihilates. Such polynomials
// have a Euclidean remainder sequence that collapses in 3 steps at any size, so every
// GCD test and benchmark built on them was exercising a degenerate case: the benchmarks
// made the half-GCD look 3x slower than the iterative version when it is in fact
// several times faster on realistic input.
func randomPolynomial(f Field, seed uint64, maxDegree int) *Polynomial {
	if maxDegree <= 0 {
		return newPolynomial(f, nil, false)
	}

	rng := rand.New(rand.NewSource(int64(seed)))

	coefficients := make([]uint64, maxDegree)
	for i := range coefficients {
		coefficients[i] = f.Reduce(rng.Uint64())
	}

	// Keep the degree deterministic: a zero leading coefficient would silently shorten
	// the polynomial and make degree-dependent assertions flaky.
	if coefficients[maxDegree-1] == 0 {
		coefficients[maxDegree-1] = 1
	}

	return newPolynomial(f, coefficients, false)
}

func BenchmarkPolyDiv(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}
	pr := NewPolyRing(f)

	p1 := randomPolynomial(f, largePrime/4, 8192)
	p2 := randomPolynomial(f, largePrime/4, 8192/2)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		pr.Div(p1, p2)
	}
}

// TODO: Optimise object creation. We spend a lot of time creating new objects.
func BenchmarkPEEA(b *testing.B) {
	f, err := NewPrimeField(largePrime)
	if err != nil {
		b.FailNow()
	}

	pr := NewPolyRing(f)

	polyMaxDegree := 8193
	p1 := randomPolynomial(f, largePrime/4, polyMaxDegree)   // Large Polynomial.
	p2 := randomPolynomial(f, largePrime/7, polyMaxDegree-1) // The degree of g0 in Gao's decoder, for a polynomial of degree p1.

	for i := 0; i <= 11; i++ {
		b.Run(fmt.Sprintf("partialGCD:remainderDeg<2^%d", i), func(b *testing.B) {
			partialDegree := 1 << i
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				pr.partialExtendedEuclidean(p1, p2, partialDegree)
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
				pr.partialExtendedEuclidean(p1, p2, (n+k)/2) // Gao's decoder partialEEA.
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

// polyProductMonicNegRoots computes \prod (x - r_i).
func polyProductMonicNegRoots(f Field, roots []uint64) *Polynomial {
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

func TestLocatorPolynomial(t *testing.T) {
	roots := makeRoots(15)

	f, err := NewPrimeField(largePrime)
	if err != nil {
		t.Fatal(err)
	}
	pr := NewPolyRing(f)

	p := polyProductMonicNegRoots(f, roots)

	intr := NewInterpolator(pr)

	q := pr.Product(intr.createMiSlice(roots))

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
				p = polyProductMonicNegRoots(f, roots)
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

		pr := NewPolyRing(f)
		// Ensuring both methods produce the same quotient and remainder.
		quo1, rem1 := pr.divSchoolbook(p.Copy(), q.Copy(), p.Degree(), q.Degree())
		quo2, rem2 := pr.divViaNTT(p.Copy(), q.Copy(), p.Degree(), q.Degree())
		a.True(quo1.Equals(quo2))
		a.True(rem1.Equals(rem2))
	}
}

/*
BenchmarkDivs/A=2048_B=1024/Div-10         	      88	  13087283 ns/op	13364703 B/op	    2056 allocs/op
BenchmarkDivs/A=2048_B=1024/DivNTT
BenchmarkDivs/A=2048_B=1024/DivNTT-10      	     356	   3374486 ns/op	  441755 B/
*/
func BenchmarkDivs(b *testing.B) {
	f, err := NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}
	pr := NewPolyRing(f)

	type cfg struct{ degA, degB int }
	cases := []cfg{
		// {16, 8},
		// {64, 32},
		// {256, 128},
		// {1024, 512},
		// {2048, 1024}, // add larger sizes as your NTT supports
		// {4096, 128},
		// {8192, 256},
		// {16384, 512},
		{1 << 15, 1 << 9},
	}

	// For stability across runs
	baseSeed := uint64(time.Now().UnixNano() / 1e6)

	for _, tc := range cases {
		name := fmt.Sprintf("A=%d_B=%d", tc.degA, tc.degB)

		// Pre-generate inputs once per size (outside timers).
		p := randomPolynomial(f, baseSeed+12345+uint64(tc.degA), tc.degA)
		q := randomPolynomial(f, baseSeed+67890+uint64(tc.degB), tc.degB)

		// Sanity check: both paths agree (outside the timer).
		quo1, rem1 := pr.divSchoolbook(p.Copy(), q.Copy(), p.Degree(), q.Degree())
		quo2, rem2 := pr.divViaNTT(p.Copy(), q.Copy(), p.Degree(), q.Degree())
		if !quo1.Equals(quo2) || !rem1.Equals(rem2) {
			b.Fatalf("mismatch for %s", name)
		}

		b.Run(name+"/Div", func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				// Copy so each iteration gets identical inputs
				_ = pr // avoid inline
				qq, rr := pr.Div(p.Copy(), q.Copy())
				// prevent compiler from optimizing away
				if qq == nil || rr == nil {
					b.Fatal("nil result")
				}
			}
		})

		b.Run(name+"/DivViaNTT", func(b *testing.B) {
			b.ReportAllocs()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				qq, rr := pr.divViaNTT(p.Copy(), q.Copy(), p.Degree(), q.Degree())
				if qq == nil || rr == nil {
					b.Fatal("nil result")
				}
			}
		})
	}
}

// Mul documents that c may alias an input, and both dispatch paths have to honour it.
// The schoolbook path used to clear c's array in place before reading the operands, so
// an aliased destination with enough capacity to hold the product silently produced the
// zero polynomial.
func TestMulAliasedDestination(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(65537)
	a.NoError(err)

	r := NewPolyRing(f)

	// degree 1 and 8 dispatch to schoolbook, 64 to the NTT path.
	for _, degree := range []int{1, 8, 64} {
		p := randomPolynomial(f, 12345+uint64(degree), degree)

		want := &Polynomial{}
		r.Mul(p, p, want)

		// Spare capacity is what made the in-place clear reachable.
		withRoom := func() *Polynomial {
			q := p.Copy()
			q.inner = append(make([]uint64, 0, 4*degree+4), q.inner...)

			return q
		}

		both := withRoom()
		r.Mul(both, both, both)
		a.True(want.Equals(both), "degree %d, c == a == b: got %s, want %s", degree, both, want)

		lhs := withRoom()
		r.Mul(lhs, p.Copy(), lhs)
		a.True(want.Equals(lhs), "degree %d, c == a: got %s, want %s", degree, lhs, want)

		rhs := withRoom()
		r.Mul(p.Copy(), rhs, rhs)
		a.True(want.Equals(rhs), "degree %d, c == b: got %s, want %s", degree, rhs, want)
	}
}
