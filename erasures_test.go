// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"math/rand"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/jonathanmweiss/go-gao/field"
)

// TestErasureSetReusedMatchesFreshOne is the correctness argument for sharing a set: a
// set built once and used for many words must give exactly what a set built per word
// would have given.
func TestErasureSetReusedMatchesFreshOne(t *testing.T) {
	const n, k = 64, 16

	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	erased := []int{3, 17, 18, 42}
	shared := mustErasures(t, code, erased...)

	rng := rand.New(rand.NewSource(11))

	for round := range 32 {
		msg := make([]uint64, k)
		for i := range msg {
			msg[i] = rng.Uint64() % code.PrimeField().Modulus()
		}

		word, err := code.Encode(msg)
		require.NoError(t, err)

		for _, i := range erased {
			word[i] = rng.Uint64() % code.PrimeField().Modulus()
		}

		fromShared, err := code.Decode(word, shared)
		require.NoError(t, err, "round %d", round)

		fromFresh, err := code.Decode(word, mustErasures(t, code, erased...))
		require.NoError(t, err, "round %d", round)

		require.Equal(t, msg, fromShared, "round %d", round)
		require.Equal(t, fromFresh, fromShared, "round %d", round)
	}
}

// TestErasureSetCrossesEquivalentCodes: the set describes the code's parameters, not the
// particular Code value, so the sender and receiver of a codeword can each build their
// own code and still share a set.
func TestErasureSetCrossesEquivalentCodes(t *testing.T) {
	const n, k = 32, 8

	newCode := func() *Code {
		c, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
		require.NoError(t, err)

		return c
	}

	sender, receiver := newCode(), newCode()
	require.NotSame(t, sender, receiver)

	msg := make([]uint64, k)
	for i := range msg {
		msg[i] = uint64(i*7 + 1)
	}

	word, err := sender.Encode(msg)
	require.NoError(t, err)

	erased := []int{1, 2, 3}
	for _, i := range erased {
		word[i] = 0xDEAD
	}

	// built on one code, spent on the other.
	got, err := receiver.Decode(word, mustErasures(t, sender, erased...))
	require.NoError(t, err)
	require.Equal(t, msg, got)
}

// TestErasureSetRejectsForeignCode pins what the parameters have to cover. The strategy
// matters as much as the numbers: the two evaluators place their points differently, so
// a set built for one is meaningless to the other even when n, k and the modulus agree.
func TestErasureSetRejectsForeignCode(t *testing.T) {
	const n, k = 16, 4

	base, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	for name, other := range map[string]func() (*Code, error){
		"different n": func() (*Code, error) {
			return NewCode(newfield(t, field.NTTFriendlyPrime), n*2, k, RequireNTT())
		},
		"different k": func() (*Code, error) {
			return NewCode(newfield(t, field.NTTFriendlyPrime), n, k+1, RequireNTT())
		},
		"different modulus": func() (*Code, error) {
			return NewCode(newfield(t, 65537), n, k, RequireNTT())
		},
		"same numbers, other strategy": func() (*Code, error) {
			return NewCode(newfield(t, field.NTTFriendlyPrime), n, k, Pointwise())
		},
	} {
		t.Run(name, func(t *testing.T) {
			foreign, err := other()
			require.NoError(t, err)

			_, err = foreign.Decode(make(Codeword, foreign.N()), mustErasures(t, base, 1, 2))
			assert.ErrorIs(t, err, ErrForeignErasureSet)
		})
	}
}

// TestPointwiseAndNTTSetsDiffer is the reason the strategy is part of the identity: the
// same indices over the same field produce different locator evaluations under the two
// evaluators.
func TestPointwiseAndNTTSetsDiffer(t *testing.T) {
	const n, k = 16, 4

	ntt, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	slow, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, Pointwise())
	require.NoError(t, err)

	assert.NotEqual(t, mustErasures(t, ntt, 1, 2).sVals, mustErasures(t, slow, 1, 2).sVals)
}

// TestZeroErasureSetDecodesAnywhere: the zero value declares nothing, so it carries no
// parameters to clash with.
func TestZeroErasureSetDecodesAnywhere(t *testing.T) {
	const k = 8

	for name, opt := range map[string]Option{"ntt": RequireNTT(), "pointwise": Pointwise()} {
		t.Run(name, func(t *testing.T) {
			code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 32, k, opt)
			require.NoError(t, err)

			msg := make([]uint64, k)
			for i := range msg {
				msg[i] = uint64(i + 1)
			}

			word, err := code.Encode(msg)
			require.NoError(t, err)

			got, err := code.Decode(word, ErasureSet{})
			require.NoError(t, err)
			require.Equal(t, msg, got)
		})
	}
}

// TestErasureSetLen reports the erasures declared, so a caller can check its own budget.
func TestErasureSetLen(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 32, 8, RequireNTT())
	require.NoError(t, err)

	assert.Equal(t, 0, ErasureSet{}.Len())
	assert.Equal(t, 3, mustErasures(t, code, 4, 5, 6).Len())

	// ranges collapse onto the symbols they touch, so Len counts symbols, not bytes.
	bc := code.Bytes()
	assert.Equal(t, 2, mustByteErasures(t, bc, ByteRange{Off: 0, Len: bc.encodedSymbolSize() + 1}).Len())
}

// TestErasureSetIsConcurrencySafe drives one shared set through many goroutines, which is
// the case the sharing exists for. Run under -race.
func TestErasureSetIsConcurrencySafe(t *testing.T) {
	const n, k = 64, 16

	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	msg := make([]uint64, k)
	for i := range msg {
		msg[i] = uint64(i*3 + 1)
	}

	word, err := code.Encode(msg)
	require.NoError(t, err)

	erased := []int{0, 9, 33}
	for _, i := range erased {
		word[i] = 7
	}

	shared := mustErasures(t, code, erased...)

	done := make(chan []uint64, 8)

	for range cap(done) {
		go func() {
			got, err := code.Decode(word, shared)
			assert.NoError(t, err)
			done <- got
		}()
	}

	for range cap(done) {
		require.Equal(t, msg, <-done)
	}
}
