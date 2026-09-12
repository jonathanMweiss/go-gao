// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/bits"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// rootTestPrimes span the range of 2-adicities that matter: 157 admits only a 4-point
// transform, 929 is the PDF417 field, 65537 is the Fermat prime this package's examples
// use, and the last is a 63-bit prime, where Pow runs on full-width operands.
var rootTestPrimes = []uint64{157, 929, 12289, 65537, 9191248642791733759}

// TestRootOfUnityHasExactOrder pins what a primitive n-th root of unity is, rather than
// which one this package happens to return: w^n = 1 and w^(n/2) != 1, so the order is n
// exactly and not a proper divisor of it.
//
// The distinction is the whole content of the derivation. Raising any element to the
// power (p-1)/n lands it somewhere in the subgroup of order n, which satisfies w^n = 1
// for free; only half of that subgroup generates it. An implementation that skipped the
// second check would return 1 for every n whenever it started from a square, and the NTT
// built on it would evaluate every point at the same place.
//
// Every admissible n is swept, up to the largest power of two dividing p-1, and one past
// it to pin the rejection.
func TestRootOfUnityHasExactOrder(t *testing.T) {
	for _, p := range rootTestPrimes {
		f, err := NewPrimeField(p)
		require.NoError(t, err)

		twoAdicity := bits.TrailingZeros64(p - 1)

		for j := 1; j <= twoAdicity; j++ {
			n := uint64(1) << j

			w, err := RootOfUnity(f, n)
			require.NoError(t, err, "p=%d admits an %d-point transform", p, n)

			assert.Equal(t, uint64(1), f.Pow(w, n), "p=%d n=%d: w^n must be 1", p, n)
			assert.NotEqual(t, uint64(1), f.Pow(w, n/2),
				"p=%d n=%d: w^(n/2) must not be 1, or the order is smaller than n", p, n)
		}

		// One past the 2-adicity, no such element exists and none may be invented.
		beyond := uint64(1) << (twoAdicity + 1)
		_, err = RootOfUnity(f, beyond)
		assert.ErrorIs(t, err, errNotDivisible, "p=%d must reject n=%d", p, beyond)
	}
}

// TestRootOfUnityRejectsUnusableSizes: the NTT in this package is radix-2, so an n that
// is not a power of two has no transform here even when a root of that order exists in
// the field. 157-1 = 156 is divisible by 3, yet a 3-point root is of no use.
func TestRootOfUnityRejectsUnusableSizes(t *testing.T) {
	f, err := NewPrimeField(157)
	require.NoError(t, err)

	_, err = RootOfUnity(f, 3)
	assert.ErrorIs(t, err, errNotPowerOfTwo, "3 divides 156 but is not a power of two")

	for _, n := range []uint64{0, 1} {
		_, err = RootOfUnity(f, n)
		assert.ErrorIs(t, err, errNSTooSmall, "n=%d is not a transform size", n)
	}
}

// TestRootOfUnityIsDeterministic: the search walks candidates from 2 upward rather than
// sampling, so a given field and size always yield the same root.
func TestRootOfUnityIsDeterministic(t *testing.T) {
	for _, tc := range []struct {
		p    uint64
		n    uint64
		want uint64
	}{
		{65537, 4, 65281},
		{65537, 8, 4096},
		{65537, 16, 64},
		{157, 4, 129},
	} {
		f, err := NewPrimeField(tc.p)
		require.NoError(t, err)

		for range 3 {
			got, err := RootOfUnity(f, tc.n)
			require.NoError(t, err)
			assert.Equal(t, tc.want, got, "p=%d n=%d", tc.p, tc.n)
		}
	}
}
