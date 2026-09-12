// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
)

// gcdPair builds a random (a, b) with deg(a) = degA and deg(b) = degA-1.
func gcdPair(pr PolyRing, rng *rand.Rand, degA int) (a, b *Polynomial) {
	coeffs := func(n int) []uint64 {
		c := make([]uint64, n)
		for i := range c {
			c[i] = rng.Uint64() % 65537
		}

		c[n-1] = rng.Uint64()%65536 + 1

		return c
	}

	return pr.NewPolynomial(coeffs(degA+1), false), pr.NewPolynomial(coeffs(degA), false)
}

// TestFastPartialGCDMatchesClassical is a regression test, and the differential check
// that should have existed from the start: the half-GCD is an optimization of the
// classical algorithm, so the two must return the same gcd and the same Bezout pair for
// the same stop degree. Nothing else pins that -- each was only ever tested against its
// own output.
//
// FastPartialGCD used to overshoot by exactly one Euclidean step for every input large
// enough to enter the half-GCD recursion (deg >= hgcdThreshold), because fastGCDRec
// asked hgcd to reach stopDegree-1 rather than stopDegree. One step too far still
// decodes while the error count leaves slack, so the decoder only misbehaved at exactly
// (n-k)/2 errors -- it silently corrected one fewer error than it advertised, at every
// codeword length from 256 up.
func TestFastPartialGCDMatchesClassical(t *testing.T) {
	_, pr := divTestField(t)
	rng := rand.New(rand.NewSource(7))

	// Straddle hgcdThreshold: below it both paths are the same code.
	for _, degA := range []int{32, 200, hgcdThreshold - 1, hgcdThreshold, hgcdThreshold + 1, 400} {
		a, b := gcdPair(pr, rng, degA)

		for stop := 1; stop < degA; stop += max(1, degA/8) {
			gFast, xFast, yFast := pr.FastPartialGCD(a, b, stop)
			gSlow, xSlow, ySlow := pr.PartialExtendedEuclidean(a, b, stop)

			require.True(t, gFast.Equals(gSlow),
				"degA=%d stop=%d: gcd differs (fast deg %d, classical deg %d)",
				degA, stop, gFast.Degree(), gSlow.Degree())
			require.True(t, xFast.Equals(xSlow), "degA=%d stop=%d: Bezout x differs", degA, stop)
			require.True(t, yFast.Equals(ySlow), "degA=%d stop=%d: Bezout y differs", degA, stop)
		}
	}
}

// FuzzGCDAgreement is the same property over arbitrary inputs and stop degrees.
func FuzzGCDAgreement(fz *testing.F) {
	fz.Add(uint64(1), uint16(32), uint16(4))
	fz.Add(uint64(2), uint16(hgcdThreshold), uint16(hgcdThreshold/2))
	fz.Add(uint64(3), uint16(hgcdThreshold+1), uint16(1))
	fz.Add(uint64(4), uint16(400), uint16(399)) // the tightest stop there is
	fz.Add(uint64(5), uint16(300), uint16(0))
	// stop == deg(a): the answer is b itself, reached in zero Euclidean steps. The two
	// implementations returned the same zero Bezout coefficient in different shapes --
	// an empty coefficient slice against a single zero -- which Polynomial.Equals used
	// to call unequal.
	fz.Add(uint64(49), uint16(11), uint16(97))

	fz.Fuzz(func(t *testing.T, seed uint64, rawDegA, rawStop uint16) {
		_, pr := divTestField(t)

		// Bounded, but wide enough to enter the recursion.
		degA := int(rawDegA)%(2*hgcdThreshold) + 2
		stop := int(rawStop) % (degA + 1)

		a, b := gcdPair(pr, rand.New(rand.NewSource(int64(seed))), degA)

		gFast, xFast, yFast := pr.FastPartialGCD(a, b, stop)
		gSlow, xSlow, ySlow := pr.PartialExtendedEuclidean(a, b, stop)

		require.True(t, gFast.Equals(gSlow),
			"degA=%d stop=%d: gcd differs (fast deg %d, classical deg %d)",
			degA, stop, gFast.Degree(), gSlow.Degree())
		require.True(t, xFast.Equals(xSlow), "degA=%d stop=%d: Bezout x differs", degA, stop)
		require.True(t, yFast.Equals(ySlow), "degA=%d stop=%d: Bezout y differs", degA, stop)
	})
}
