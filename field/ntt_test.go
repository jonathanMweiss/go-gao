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
