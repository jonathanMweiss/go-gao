package field

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// TestHGCDSecondRecursivePath forces the second recursive call in hgcd and fastGCDRec.
func TestHGCDSecondRecursivePath(t *testing.T) {
	// Use a standard prime field.
	f, err := NewPrimeField(65537)
	if err != nil {
		t.Fatal(err)
	}

	pr := NewDensePolyRing(f).(*DensePolyRing)

	t.Run("VeryLargeRandom", func(t *testing.T) {
		// Use a huge N to ensure we hit all recursive branches.
		n := 16384
		stopDegree := 8192

		p := randomPolynomial(f, 42, n)
		q := randomPolynomial(f, 24, n-1)

		gcd, x, y := pr.FastPartialGCD(p, q, stopDegree)
		assert.True(t, bezoutIdentityHolds(pr, p, q, gcd, x, y), "Bézout identity should hold for Fast GCD")
	})
}
