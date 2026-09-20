// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// randomPoly returns a polynomial of exactly the given degree.
func randomPoly(t testing.TB, f Field, pr *PolyRing, rng *rand.Rand, degree int) *Polynomial {
	t.Helper()

	coeffs := make([]uint64, degree+1)
	for i := range coeffs {
		coeffs[i] = rng.Uint64() % f.Modulus()
	}

	// a polynomial of degree d must not have a zero leading coefficient.
	for f.Equals(coeffs[degree], 0) {
		coeffs[degree] = rng.Uint64() % f.Modulus()
	}

	return pr.NewPolynomial(coeffs, false)
}

// TestSeriesInverseTruncates pins the property the cache rests on: an inverse computed to
// precision k is already an inverse to every lower precision, so one inverse at the
// longest quotient length covers every shorter division.
func TestSeriesInverseTruncates(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(3))

	for _, k := range []int{1, 2, 5, 16, 64, 129} {
		for _, bLen := range []int{1, 2, 9, 40} {
			coeffs := make([]uint64, bLen)
			for i := range coeffs {
				coeffs[i] = rng.Uint64() % f.Modulus()
			}

			// the constant term must be invertible for the series to exist.
			coeffs[0] = rng.Uint64()%65536 + 1

			b := pr.NewPolynomial(coeffs, false)
			full := pr.seriesInverse(b, k).ToSlice()

			for j := 1; j <= k; j++ {
				want := pr.seriesInverse(b, j).ToSlice()

				got := full
				if j < len(got) {
					got = got[:j]
				}

				// ToSlice drops trailing zeros, so the shorter one bounds the comparison.
				n := min(len(want), len(got))
				require.Equal(t, want[:n], got[:n], "k=%d len(b)=%d j=%d", k, bLen, j)
			}
		}
	}
}

// TestDivByMatchesDiv is the whole correctness argument for the cache: reusing an
// inverse must give exactly what recomputing it gives, at every quotient length on
// either side of what the cache was built for.
func TestDivByMatchesDiv(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(7))

	for _, divisorDeg := range []int{1, 5, 16, 64, 257} {
		divisor := randomPoly(t, f, pr, rng, divisorDeg)

		const maxQuot = 128
		dc := pr.NewDivisorCache(divisor, maxQuot)

		// quotient lengths below, at, and beyond what the cache covers, plus the
		// degenerate case of a dividend smaller than the divisor.
		for _, quotLen := range []int{1, 2, 17, 64, maxQuot - 1, maxQuot, maxQuot + 1, maxQuot * 3} {
			t.Run(fmt.Sprintf("divisor=%d/quot=%d", divisorDeg, quotLen), func(t *testing.T) {
				dividend := randomPoly(t, f, pr, rng, divisorDeg+quotLen-1)

				wantQ, wantR := pr.Div(dividend, divisor)
				gotQ, gotR := pr.DivBy(dividend, dc)

				assert.Equal(t, wantQ.ToSlice(), gotQ.ToSlice(), "quotient")
				assert.Equal(t, wantR.ToSlice(), gotR.ToSlice(), "remainder")
			})
		}

		t.Run(fmt.Sprintf("divisor=%d/dividend smaller", divisorDeg), func(t *testing.T) {
			if divisorDeg == 0 {
				t.Skip("nothing is smaller than a constant")
			}

			dividend := randomPoly(t, f, pr, rng, divisorDeg-1)

			wantQ, wantR := pr.Div(dividend, divisor)
			gotQ, gotR := pr.DivBy(dividend, dc)

			assert.Equal(t, wantQ.ToSlice(), gotQ.ToSlice(), "quotient")
			assert.Equal(t, wantR.ToSlice(), gotR.ToSlice(), "remainder")
		})
	}
}

// TestDivByExactDivision covers the case the decoder actually takes: the divisor goes in
// a whole number of times, so the remainder must be zero and the quotient exact.
func TestDivByExactDivision(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(11))

	divisor := randomPoly(t, f, pr, rng, 40)
	dc := pr.NewDivisorCache(divisor, 256)

	for _, quotDeg := range []int{0, 1, 30, 200, 255} {
		quotient := randomPoly(t, f, pr, rng, quotDeg)

		product := pr.newDst()
		pr.Mul(quotient, divisor, product)

		gotQ, gotR := pr.DivBy(product, dc)

		assert.True(t, gotR.IsZero(), "quotDeg=%d: remainder must vanish", quotDeg)
		assert.Equal(t, quotient.ToSlice(), gotQ.ToSlice(), "quotDeg=%d", quotDeg)
	}
}

// TestDivByLeavesInputsAlone: the cache is shared between dividends and across
// goroutines, so a division must not write through to it or to the dividend.
func TestDivByLeavesInputsAlone(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(13))

	divisor := randomPoly(t, f, pr, rng, 20)
	dc := pr.NewDivisorCache(divisor, 128)

	cachedBefore := dc.t.ToSlice()
	divisorBefore := divisor.ToSlice()

	dividend := randomPoly(t, f, pr, rng, 100)
	dividendBefore := dividend.ToSlice()

	pr.DivBy(dividend, dc)

	assert.Equal(t, cachedBefore, dc.t.ToSlice(), "cached inverse")
	assert.Equal(t, divisorBefore, divisor.ToSlice(), "divisor")
	assert.Equal(t, dividendBefore, dividend.ToSlice(), "dividend")
}

// TestDivisorCacheCopiesDivisor: a caller that reuses its polynomial buffer must not be
// able to corrupt a cache built earlier.
func TestDivisorCacheCopiesDivisor(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(17))

	divisor := randomPoly(t, f, pr, rng, 12)
	dc := pr.NewDivisorCache(divisor, 64)

	dividend := randomPoly(t, f, pr, rng, 60)
	wantQ, wantR := pr.DivBy(dividend, dc)

	// scribble over the caller's copy.
	inner := divisor.NoCopySlice()
	for i := range inner {
		inner[i] = f.Reduce(uint64(i + 1))
	}

	gotQ, gotR := pr.DivBy(dividend, dc)

	assert.Equal(t, wantQ.ToSlice(), gotQ.ToSlice())
	assert.Equal(t, wantR.ToSlice(), gotR.ToSlice())
}

// TestDivByIsConcurrencySafe drives one cache through many goroutines, which is what an
// ErasureSet shared across a batch does. Run under -race.
func TestDivByIsConcurrencySafe(t *testing.T) {
	f, pr := divTestField(t)

	rng := rand.New(rand.NewSource(19))

	divisor := randomPoly(t, f, pr, rng, 32)
	dc := pr.NewDivisorCache(divisor, 256)

	dividend := randomPoly(t, f, pr, rng, 200)
	wantQ, wantR := pr.Div(dividend, divisor)

	done := make(chan struct{}, 8)

	for range cap(done) {
		go func() {
			defer func() { done <- struct{}{} }()

			// each goroutine needs its own ring: a PolyRing carries scratch buffers.
			ring := NewPolyRing(f)

			gotQ, gotR := ring.DivBy(dividend, dc)
			assert.Equal(t, wantQ.ToSlice(), gotQ.ToSlice())
			assert.Equal(t, wantR.ToSlice(), gotR.ToSlice())
		}()
	}

	for range cap(done) {
		<-done
	}
}

// TestNewDivisorCacheRejectsBadDivisors mirrors what Div panics on.
func TestNewDivisorCacheRejectsBadDivisors(t *testing.T) {
	f, pr := divTestField(t)

	assert.Panics(t, func() { pr.NewDivisorCache(nil, 8) })
	assert.Panics(t, func() { pr.NewDivisorCache(pr.NewPolynomial(nil, false), 8) })

	ntt := pr.NewPolynomial(make([]uint64, 16), false)
	require.NoError(t, pr.NttForward(ntt))
	assert.Panics(t, func() { pr.NewDivisorCache(ntt, 8) })

	_ = f
}

// FuzzDivByMatchesDiv is the differential test widened: whatever the shapes, reusing a
// cached inverse must agree with recomputing it, including where the cache declines and
// falls back.
func FuzzDivByMatchesDiv(fz *testing.F) {
	fz.Add(uint64(1), uint8(40), uint8(9), uint8(16), uint8(0))
	fz.Add(uint64(2), uint8(64), uint8(3), uint8(4), uint8(7))   // quotient past maxQuot
	fz.Add(uint64(3), uint8(20), uint8(20), uint8(64), uint8(0)) // equal degrees
	fz.Add(uint64(4), uint8(3), uint8(9), uint8(32), uint8(0))   // divisor larger
	fz.Add(uint64(5), uint8(48), uint8(2), uint8(1), uint8(11))  // low-order zeros

	fz.Fuzz(func(t *testing.T, seed uint64, aLen, bLen, quotCap, lowZeros uint8) {
		_, pr := divTestField(t)

		na := int(aLen)%96 + 1
		nb := int(bLen)%32 + 1
		maxQuot := int(quotCap)%96 + 1

		rng := rand.New(rand.NewSource(int64(seed)))

		aCoeffs := make([]uint64, na)
		for i := int(lowZeros) % na; i < na; i++ {
			aCoeffs[i] = rng.Uint64() % NTTFriendlyPrime
		}

		bCoeffs := make([]uint64, nb)
		for i := range bCoeffs {
			bCoeffs[i] = rng.Uint64() % NTTFriendlyPrime
		}

		// a zero divisor is a documented panic, not a case to check here.
		bCoeffs[nb-1] = rng.Uint64()%65536 + 1

		a := pr.NewPolynomial(aCoeffs, false)
		b := pr.NewPolynomial(bCoeffs, false)

		wantQ, wantR := pr.Div(a, b)
		gotQ, gotR := pr.DivBy(a, pr.NewDivisorCache(b, maxQuot))

		require.Equal(t, wantQ.ToSlice(), gotQ.ToSlice(), "quotient")
		require.Equal(t, wantR.ToSlice(), gotR.ToSlice(), "remainder")
	})
}
