// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

// TestHGCDSecondRecursivePath forces the second recursive call in hgcd -- the one that
// advances the remaining distance after the mandatory division. At these sizes it is
// taken 63 times.
func TestHGCDSecondRecursivePath(t *testing.T) {
	// Use a standard prime field.
	f, err := NewPrimeField(65537)
	if err != nil {
		t.Fatal(err)
	}

	pr := NewPolyRing(f)

	t.Run("VeryLargeRandom", func(t *testing.T) {
		// Use a huge N to ensure we hit all recursive branches.
		n := 16384
		stopDegree := 8192

		p := randomPolynomial(f, 42, n)
		q := randomPolynomial(f, 24, n-1)

		gcd, x, y := pr.PartialGCD(p, q, stopDegree)
		assert.True(t, bezoutIdentityHolds(pr, p, q, gcd, x, y), "Bézout identity should hold for Fast GCD")
	})
}
