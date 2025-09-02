package field

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestNTTForward(t *testing.T) {
	a := assert.New(t)
	f, err := NewPrimeField(3329)
	a.NoError(err)

	p := NewPolynomial(f, []uint64{1, 2, 3, 4, 5, 6, 7, 8}, false)
	expected := []uint64{36, 3240, 3067, 427, 3325, 2894, 254, 81}

	pr := NewDensePolyRing(f)

	a.NoError(pr.NttForward(p))
	a.Equal(expected, p.ToSlice())
}

func TestNTTForwardBackward(t *testing.T) {
	// Test the forward and backward NTT transforms.
	a := assert.New(t)
	f, err := NewPrimeField(65537)
	a.NoError(err)

	pr := NewDensePolyRing(f)
	for i := range 8 {
		cappingDegree := 1 << (i + 1)

		p1 := randomPolynomial(f, 12345+uint64(i), cappingDegree)
		pcpy := p1.Copy()

		a.NoError(pr.NttForward(p1))

		a.NoError(pr.NttBackward(p1))

		a.True(pcpy.Equals(p1))
	}
}

func TestPolyMult(t *testing.T) {
	a := assert.New(t)
	f, err := NewPrimeField(65537)
	a.NoError(err)

	pr := NewDensePolyRing(f)

	for i := 0; i < 8; i++ {

		degree := 1 << (i + 1)

		p1 := randomPolynomial(f, 12345+uint64(i), degree)

		regMulRes := &Polynomial{}
		pr.MulPoly(p1, p1, regMulRes)

		// padding p1 with zeros:
		p1.inner = append(p1.inner, make([]uint64, degree)...)

		a.NoError(pr.NttForward(p1))

		nttRes := &Polynomial{}
		pr.MulPoly(p1, p1, nttRes)

		a.NoError(pr.NttBackward(nttRes))

		if !regMulRes.Equals(nttRes) {
			regMulRes.Equals(nttRes)
			t.Errorf("Mismatch between regular and NTT multiplication results")
		}
		a.True(regMulRes.Equals(nttRes))
	}
}
