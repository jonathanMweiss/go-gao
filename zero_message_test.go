// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// zeroMessageCases covers both evaluators, since the degenerate path is reached
// identically through decodeNTT and decodeGeneric.
var zeroMessageCases = []testCase{
	{"pointwise", Pointwise(), 18, 5},
	{"ntt", nil, 16, 4},
}

// TestZeroMessageDecodes pins zeroCodewordIsNearest.
//
// The all-zero message is degenerate for Gao's algorithm: with f = 0 the Berlekamp-Welch
// product Q = E*f is identically zero, so it carries no degree for the partial GCD to
// stop at. zeroCodewordIsNearest settles the case up front instead.
//
// The degeneracy belongs to the message, not to the received word. Every word inside the
// decoding radius of the all-zero codeword decodes to f = 0, and every one of them used
// to fail: the exact-zero word by panicking on a division by the zero polynomial, and
// all the rest by returning ErrDecoding. So the check is 2*errors+erasures <= n-k
// measured against the zero codeword, not merely "is the received word all zero" --
// narrowing it to the latter reintroduces the ErrDecoding half of the bug.
//
// What is pinned here is the whole ball, with the boundary spent in each of its three
// ways: all errors, all erasures, and a mixture.
func TestZeroMessageDecodes(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range zeroMessageCases {
		code, err := NewCode(f, tc.n, tc.k, tc.opt)
		require.NoError(t, err)

		zero := make([]uint64, tc.k)
		budget := tc.n - tc.k

		// corrupt returns the all-zero codeword damaged by errs errors and erasures
		// erasures. Every value written differs from the zero it replaces, so the error
		// count zeroCodewordIsNearest sees is exactly errs.
		// t is a parameter rather than captured: require's FailNow must act on the
		// subtest that is running, not on the parent.
		corrupt := func(t *testing.T, errs, erasures int) ([]uint64, []int) {
			codeword, err := code.Encode(zero)
			require.NoError(t, err)

			return codeword, damageCodeword(f, testRNG(t), codeword, errs, erasures)
		}

		t.Run(tc.name+"/untouched", func(t *testing.T) {
			codeword, erased := corrupt(t, 0, 0)

			decoded, err := code.Decode(codeword, erased...)
			require.NoError(t, err)
			assert.Equal(t, zero, decoded, "must return k zeros, not an empty slice")
		})

		t.Run(tc.name+"/errors at the budget", func(t *testing.T) {
			codeword, erased := corrupt(t, code.MaxErrors(), 0)

			decoded, err := code.Decode(codeword, erased...)
			require.NoError(t, err, "%d errors is exactly MaxErrors", code.MaxErrors())
			assert.Equal(t, zero, decoded)
		})

		t.Run(tc.name+"/erasures at the budget", func(t *testing.T) {
			codeword, erased := corrupt(t, 0, budget)

			decoded, err := code.Decode(codeword, erased...)
			require.NoError(t, err, "%d erasures is exactly n-k", budget)
			assert.Equal(t, zero, decoded)
		})

		t.Run(tc.name+"/mixed at the budget", func(t *testing.T) {
			erasures := budget / 2
			errs := (budget - erasures) / 2

			codeword, erased := corrupt(t, errs, erasures)

			decoded, err := code.Decode(codeword, erased...)
			require.NoError(t, err, "2*%d+%d is within %d", errs, erasures, budget)
			assert.Equal(t, zero, decoded)
		})

		// One past the boundary the zero message is no longer guaranteed, so the
		// shortcut must decline and let the ordinary decoder answer. The ordinary
		// decoder cannot decode f = 0 either -- what matters is that it reports that
		// rather than panicking on the zero-polynomial division underneath.
		t.Run(tc.name+"/past the budget declines", func(t *testing.T) {
			codeword, erased := corrupt(t, code.MaxErrors()+1, 0)

			require.NotPanics(t, func() {
				_, err := code.Decode(codeword, erased...)
				assert.ErrorIs(t, err, ErrDecoding)
			})
		})
	}
}

// TestNearZeroMessagesRoundTrip: a message that is mostly but not entirely zero is not
// degenerate, and must go through the ordinary decoder untouched by the shortcut. The
// low-order case is the one that used to break division, and the high-order case the one
// that used to come back shorter than k.
func TestNearZeroMessagesRoundTrip(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range zeroMessageCases {
		code, err := NewCode(f, tc.n, tc.k, tc.opt)
		require.NoError(t, err)

		messages := map[string][]uint64{
			"only the lowest symbol":  append([]uint64{7}, make([]uint64, tc.k-1)...),
			"only the highest symbol": append(make([]uint64, tc.k-1), 7),
			"one in the middle":       append(append(make([]uint64, tc.k-2), 7), 0),
		}

		for name, data := range messages {
			t.Run(tc.name+"/"+name, func(t *testing.T) {
				codeword, err := code.Encode(data)
				require.NoError(t, err)

				corruptCodeword(f, testRNG(t), codeword, code.MaxErrors())

				decoded, err := code.Decode(codeword)
				require.NoError(t, err)
				assert.Equal(t, data, decoded)
			})
		}
	}
}
