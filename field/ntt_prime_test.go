// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/big"
	"math/bits"
	"testing"

	"github.com/stretchr/testify/require"
)

// TestNTTFriendlyPrimeProperties pins what the constant promises. Both properties are
// silent if they break: a non-prime modulus still computes, just not in a field, and a
// short 2-adicity only shows up as a construction error at some size nobody tested.
func TestNTTFriendlyPrimeProperties(t *testing.T) {
	const p uint64 = NTTFriendlyPrime

	require.True(t, (&big.Int{}).SetUint64(p).ProbablyPrime(20), "modulus must be prime")

	// Above 2^56, so a symbol carries seven whole bytes.
	require.Greater(t, p, uint64(1)<<56, "must pack at least 7 bytes per symbol")

	// Below the Shoup ceiling: mulShoup corrects with one conditional subtraction, which
	// is sufficient only while the pre-correction residue stays under 2p, so 2p must not
	// overflow a uint64.
	require.Less(t, p, uint64(1)<<63, "2p must fit in a uint64")

	// 2-adicity: the largest power of two dividing p-1 is what bounds transform length.
	adicity := bits.TrailingZeros64(p - 1)
	require.GreaterOrEqual(t, adicity, 20, "2-adicity too small for large transforms")

	f := newPrimeField(t, p)

	// Every power-of-two length up to the 2-adicity really does yield a root of that
	// exact order.
	for k := 1; k <= adicity && k <= 20; k++ {
		n := uint64(1) << k

		w, err := RootOfUnity(f, n)
		require.NoError(t, err, "n=2^%d", k)

		require.Equal(t, uint64(1), f.Pow(w, n), "w^n must be 1 at n=2^%d", k)
		require.NotEqual(t, uint64(1), f.Pow(w, n/2), "order must be exactly n at n=2^%d", k)
	}
}
