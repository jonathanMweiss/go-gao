// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"math/rand"
	"testing"

	"github.com/stretchr/testify/require"
)

// divTestField is large enough that the NTT convolution path inside Div is reachable at
// every size these tests use.
func divTestField(t testing.TB) (Field, PolyRing) {
	t.Helper()

	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	return f, NewDensePolyRing(f)
}

// TestDivDividendWithLowOrderZeros is a regression test.
//
// Div hands anything with a long enough quotient to the Newton-iteration path, which
// reverses the quotient series to recover the quotient. That reversal used to go around
// the series' true degree rather than the quotient length, so a dividend whose low-order
// coefficients were zero came back with a quotient shifted down by exactly that many
// degrees -- and a non-zero remainder. Dividing by the constant 1 did not return the
// dividend.
//
// The decoder hit this whenever the message it was recovering had zero low-order
// symbols, which is to say routinely.
func TestDivDividendWithLowOrderZeros(t *testing.T) {
	_, pr := divTestField(t)

	const degree = 16

	for _, lowZeros := range []int{0, 1, 2, 3, 7, degree} {
		coeffs := make([]uint64, degree+1)
		for i := lowZeros; i <= degree; i++ {
			coeffs[i] = uint64(i*7 + 1)
		}

		a := pr.NewPolynomial(coeffs, false)

		for _, divisor := range [][]uint64{{1}, {5}, {3, 1}, {0, 0, 1}} {
			b := pr.NewPolynomial(divisor, false)

			q, rem := pr.Div(a, b)

			// a == q*b + rem, and deg(rem) < deg(b).
			got := polyAdd(pr.(*DensePolyRing), polyMul(pr.(*DensePolyRing), q, b), rem)
			require.Equal(t, a.Degree(), got.Degree(),
				"lowZeros=%d divisor=%v: q*b+rem must reproduce a", lowZeros, divisor)
			require.True(t, a.Equals(got),
				"lowZeros=%d divisor=%v: q*b+rem must reproduce a", lowZeros, divisor)
			require.Less(t, rem.Degree(), b.Degree(),
				"lowZeros=%d divisor=%v: remainder must be smaller than the divisor", lowZeros, divisor)
		}
	}
}

// TestDivByOneIsIdentity is the smallest statement of the same bug.
func TestDivByOneIsIdentity(t *testing.T) {
	_, pr := divTestField(t)

	a := pr.NewPolynomial([]uint64{0, 0, 0, 14629, 3975}, false)
	one := pr.NewPolynomial([]uint64{1}, false)

	q, rem := pr.Div(a, one)

	require.True(t, rem.IsZero(), "dividing by 1 leaves no remainder")
	require.True(t, a.Equals(q), "dividing by 1 returns the dividend")
}

// FuzzDiv pins the defining property of division: a = q*b + rem with deg(rem) < deg(b),
// across both the schoolbook and the Newton paths that Div dispatches between.
func FuzzDiv(fz *testing.F) {
	fz.Add(uint64(1), uint8(16), uint8(1), uint8(0))
	fz.Add(uint64(2), uint8(16), uint8(1), uint8(3))  // low-order zeros in the dividend
	fz.Add(uint64(3), uint8(40), uint8(9), uint8(11)) // well past the NTT threshold
	fz.Add(uint64(4), uint8(20), uint8(20), uint8(5)) // equal degrees
	fz.Add(uint64(5), uint8(3), uint8(9), uint8(0))   // divisor larger than dividend

	fz.Fuzz(func(t *testing.T, seed uint64, aLen, bLen, lowZeros uint8) {
		_, pr := divTestField(t)
		ring := pr.(*DensePolyRing)

		// Keep the sizes bounded: this is about correctness, not throughput.
		na := int(aLen)%64 + 1
		nb := int(bLen)%32 + 1

		rng := rand.New(rand.NewSource(int64(seed)))

		aCoeffs := make([]uint64, na)
		for i := int(lowZeros) % na; i < na; i++ {
			aCoeffs[i] = rng.Uint64() % 65537
		}

		bCoeffs := make([]uint64, nb)
		for i := range bCoeffs {
			bCoeffs[i] = rng.Uint64() % 65537
		}

		// A zero divisor is a documented panic, not a case to check here.
		bCoeffs[nb-1] = rng.Uint64()%65536 + 1

		a := pr.NewPolynomial(aCoeffs, false)
		b := pr.NewPolynomial(bCoeffs, false)

		q, rem := pr.Div(a, b)

		require.Less(t, rem.Degree(), b.Degree(), "remainder must be smaller than the divisor")

		got := polyAdd(ring, polyMul(ring, q, b), rem)
		require.True(t, a.Equals(got), "q*b + rem must reproduce a\n a  =%v\n got=%v", a, got)
	})
}
