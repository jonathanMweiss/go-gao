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

	t.Run("simple", func(t *testing.T) {
		m1 := NewPolynomial(f, []uint64{5, 1}, false)
		m2 := NewPolynomial(f, []uint64{3, 1}, false)

		m := m1.Mul(m2)

		q, r := m.LongDiv(m1)

		a.Equal(makeConstantPoly(f, 0).ToSlice(), r.ToSlice())
		a.Equal(m2.ToSlice(), q.ToSlice())

		intr := NewInterpolator(f)

		q_ := intr.mDivMi(m, m1)
		a.Equal(q.ToSlice(), q_.ToSlice())

		q, r = m.LongDiv(m2)
		a.Equal(makeConstantPoly(f, 0).ToSlice(), r.ToSlice())
		a.Equal(m1.ToSlice(), q.ToSlice())

		q_ = intr.mDivMi(m, m2)
		a.Equal(q.ToSlice(), q_.ToSlice())
	})

	t.Run("complex", func(t *testing.T) {
		xs := []uint64{1, 2, 3, 5, 6, 7}

		intr := NewInterpolator(f)

		miSlice := intr.createMiSlice(xs)
		m := PolyProduct(f, miSlice)

		for _, mi := range miSlice {
			qQuickDiv := intr.mDivMi(m, mi)
			qLongdiv, _ := m.LongDiv(mi)
			a.Equal(qQuickDiv.ToSlice(), qLongdiv.ToSlice())
		}
	})
}

func TestInterpolation(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(157)
	a.NoError(err)

	coeffs := []uint64{0, 1, 2}
	p := NewPolynomial(f, coeffs, false)

	intr := NewInterpolator(f)

	xs, ys := evalPolyForTest(p, 0, 3)

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

	f.Fuzz(func(t *testing.T, randomSeed uint64) {
		a := assert.New(t)
		const boundingDegree = 10

		p := randomPolynomial(fld, randomSeed, boundingDegree)

		// interpolate a random polynomial
		intr := NewInterpolator(fld)

		xs, ys := evalPolyForTest(p, int(randomSeed), boundingDegree)
		q, err := intr.Interpolate(xs, ys)
		a.NoError(err)

		fmt.Println()
		a.Equal(p.ToSlice(), q.ToSlice())
	})

}

func evalPolyForTest(p *Polynomial, randomSeed, numEvals int) ([]uint64, []uint64) {
	xs := make([]uint64, numEvals)
	for i := range xs {
		xs[i] = p.f.Reduce(uint64(randomSeed + i + 1))
	}

	ys := make([]uint64, len(xs))

	for i, x := range xs {
		ys[i] = p.Eval(x)
	}

	return xs, ys
}

func BenchmarkMDivMi(b *testing.B) {
	a := assert.New(b)

	f, err := NewPrimeField(157)
	a.NoError(err)

	xs := []uint64{1, 2, 3, 5, 6, 7}

	intr := NewInterpolator(f)

	miSlice := intr.createMiSlice(xs)
	m := PolyProduct(f, miSlice)

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
			m.LongDiv(mi)
		}
	})
}
