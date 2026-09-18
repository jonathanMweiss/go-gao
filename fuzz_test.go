// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/require"
)

// fuzzParams are the two shapes FuzzRoundTrip alternates between: the NTT path and the
// pointwise one, so a difference between them shows up as a failure rather than as an
// untested branch.
var fuzzParams = map[bool]struct {
	opt  Option
	n, k int
}{
	true:  {RequireNTT(), 16, 4},
	false: {Pointwise(), 18, 5},
}

// fuzzMessage builds the message for a seed. Roughly a quarter of the symbols are zero,
// and seed 0 is the all-zero message, so both the high-order-zero case and the
// degenerate all-zero case are reached often rather than by luck.
func fuzzMessage(rng *rand.Rand, seed uint64, k int, modulus uint64) []uint64 {
	data := make([]uint64, k)
	if seed == 0 {
		return data
	}

	for i := range data {
		if rng.Intn(4) == 0 {
			continue
		}

		data[i] = rng.Uint64() % modulus
	}

	return data
}

// FuzzRoundTrip is the regression net this package was missing. Every hand-written test
// here encodes a message of distinct non-zero symbols, which is exactly the shape that
// hides two real bugs: a message whose high-order symbols are zero used to decode to a
// slice shorter than k, and the all-zero message used to panic as soon as an erasure was
// declared.
//
// The property: inside the budget 2*errors+erasures <= n-k, what comes out of the
// decoder equals what went into the encoder — same values, same length — through either
// API.
func FuzzRoundTrip(fz *testing.F) {
	seeds := []struct {
		msgSeed, corruptSeed uint64
		numErrors, numErases uint8
		useNTT               bool
	}{
		{1, 1, 0, 0, true},  // clean round trip
		{1, 1, 0, 0, false}, // clean round trip, pointwise
		{0, 1, 0, 0, true},  // all-zero message, untouched
		{0, 1, 0, 1, true},  // all-zero message + erasure: used to panic
		{0, 1, 1, 0, true},  // all-zero message + error: used to return ErrDecoding
		{0, 1, 0, 1, false}, // same, pointwise
		{2, 9, 6, 0, true},  // errors only, at the budget
		{2, 9, 0, 12, true}, // erasures only, at the budget
		{3, 4, 3, 6, true},  // mixed, at the budget
		{3, 4, 4, 5, false}, // mixed, at the budget, pointwise
		{7, 11, 2, 3, true}, // comfortably inside the budget
		{255, 255, 255, 255, true},
	}

	for _, s := range seeds {
		fz.Add(s.msgSeed, s.corruptSeed, s.numErrors, s.numErases, s.useNTT)
	}

	fz.Fuzz(func(t *testing.T, msgSeed, corruptSeed uint64, numErrors, numErases uint8, useNTT bool) {
		f := newfield(t, field.NTTFriendlyPrime)

		p := fuzzParams[useNTT]

		code, err := NewCode(f, p.n, p.k, p.opt)
		require.NoError(t, err)
		require.Equal(t, useNTT, code.UsesNTT())

		budget := p.n - p.k

		// Spend the budget: each erasure costs one, each error two.
		erasures := int(numErases) % (budget + 1)

		errors := 0
		if room := (budget - erasures) / 2; room > 0 {
			errors = int(numErrors) % (room + 1)
		}

		data := fuzzMessage(rand.New(rand.NewSource(int64(msgSeed))), msgSeed, p.k, f.Modulus())

		codeword, err := code.Encode(data)
		require.NoError(t, err)

		rng := rand.New(rand.NewSource(int64(corruptSeed)))
		erasedAt := damageCodeword(f, rng, codeword, errors, erasures)

		decoded, err := code.Decode(codeword, erasedAt...)
		require.NoError(t, err, "%d errors + %d erasures is within the budget of %d", errors, erasures, budget)
		require.Equal(t, data, decoded, "Decode must return the message unchanged, at length k")

		// Decoding is a pure function of the codeword and the erasure list: the same
		// input decodes the same way, and the caller's slice survives it.
		again, err := code.Decode(codeword, erasedAt...)
		require.NoError(t, err)
		require.Equal(t, decoded, again, "Decode must not depend on or disturb its input")
	})
}
