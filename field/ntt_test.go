package field

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

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
		p1, err = pr.NttForward(p1)
		a.NoError(err)

		p1, err = pr.NttBackward(p1)
		a.NoError(err)

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

		p1, err = pr.NttForward(p1)
		a.NoError(err)

		nttRes := &Polynomial{}
		pr.MulPoly(p1, p1, nttRes)

		nttRes, err = pr.NttBackward(nttRes)
		a.NoError(err)

		if !regMulRes.Equals(nttRes) {
			regMulRes.Equals(nttRes)
			t.Errorf("Mismatch between regular and NTT multiplication results")
		}
		a.True(regMulRes.Equals(nttRes))
	}
}
