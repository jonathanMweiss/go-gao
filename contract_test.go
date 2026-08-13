// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"maps"
	"slices"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewCodeRejectsBadSizes covers parameters that used to panic from inside Encode
// instead of being reported at construction time.
func TestNewCodeRejectsBadSizes(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range []struct {
		name string
		opts []Option
		n, k int
		want error
	}{
		{"n smaller than k", nil, 4, 16, ErrNSmallerThanK},
		{"zero k", nil, 16, 0, ErrNonPositiveK},
		{"negative k", nil, 16, -1, ErrNonPositiveK},
		{"n exceeds the field", nil, 1 << 20, 4, ErrUnsupportedSize},
		// RequireNTT turns the silent pointwise fallback into an error.
		{"require ntt, n not a power of two", []Option{RequireNTT()}, 20, 4, ErrUnsupportedSize},
		{"require ntt, n does not divide p-1", []Option{RequireNTT()}, 1 << 20, 4, ErrUnsupportedSize},
		{"require ntt, n too small", []Option{RequireNTT()}, 1, 1, ErrUnsupportedSize},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := NewCode(f, tc.n, tc.k, tc.opts...)
			assert.ErrorIs(t, err, tc.want)
		})
	}
}

// TestNewCodeSelectsStrategy: the evaluator is chosen for the caller, so the choice must
// be right — and visible, since falling back to the quadratic path is otherwise silent.
func TestNewCodeSelectsStrategy(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	// 929 is the prime field PDF417 barcodes are defined over. p-1 = 928 = 2^5 * 29, so
	// an NTT exists only up to n=32; n=100 must fall back.
	pdf417, err := field.NewPrimeField(929)
	require.NoError(t, err)

	for _, tc := range []struct {
		name    string
		f       field.Field
		opts    []Option
		n, k    int
		wantNTT bool
	}{
		{"prefers ntt when possible", f, nil, 16, 4, true},
		{"ntt at the field's limit", f, nil, 1 << 16, 4, true},
		{"falls back when n is not a power of two", f, nil, 18, 5, false},
		{"Pointwise overrides an available ntt", f, []Option{Pointwise()}, 16, 4, false},
		{"GF(929) small enough for ntt", pdf417, nil, 32, 8, true},
		{"GF(929) beyond its ntt range", pdf417, nil, 100, 60, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			code, err := NewCode(tc.f, tc.n, tc.k, tc.opts...)
			require.NoError(t, err)
			assert.Equal(t, tc.wantNTT, code.UsesNTT())
			assert.Equal(t, tc.n, code.N())
			assert.Equal(t, tc.k, code.K())
		})
	}
}

// TestRequireNTTErrorIsActionable: a user who set RequireNTT hit the error precisely
// because they did not know their parameters were unsuitable, so the message has to say
// what the constraint is.
func TestRequireNTTErrorIsActionable(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	_, err = NewCode(f, 20, 4, RequireNTT())
	require.Error(t, err)

	assert.Contains(t, err.Error(), "power of two")
	assert.Contains(t, err.Error(), "65537")
}

// TestDecodeDoesNotMutateInput pins the contract that decoding leaves the caller's data
// alone. Decode used to write zeros into the map for every erased point, and
// DecodeFromSlice used to reduce and transform the caller's slice in place.
func TestDecodeDoesNotMutateInput(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	} {
		code, err := NewCode(f, tc.n, tc.k, tc.opt)
		require.NoError(t, err)

		data := makeTestSlice(tc.k)

		t.Run("Decode/map", func(t *testing.T) {
			encoded, err := code.Encode(data)
			require.NoError(t, err)

			// Drop the maximum number of points, so the erasure path runs.
			xs := code.EvaluationPoints()
			for i := 0; i < code.N()-code.K(); i++ {
				delete(encoded, xs[i])
			}

			before := maps.Clone(encoded)

			decoded, err := code.Decode(encoded)
			require.NoError(t, err)
			assert.Equal(t, data, decoded)

			assert.Equal(t, before, encoded, "Decode must not modify the map it is given")
		})

		t.Run("DecodeFromSlice", func(t *testing.T) {
			codeword, err := code.EncodeToSlice(data)
			require.NoError(t, err)

			before := slices.Clone(codeword)

			decoded, err := code.DecodeFromSlice(codeword)
			require.NoError(t, err)
			assert.Equal(t, data, decoded)

			assert.Equal(t, before, codeword, "DecodeFromSlice must not modify the slice it is given")
		})
	}
}

// TestDecodeFromSliceRejectsWrongLength: the positional API cannot express erasures, so
// a short slice is a caller error rather than a set of missing points.
func TestDecodeFromSliceRejectsWrongLength(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	code, err := NewCode(f, 16, 4, RequireNTT())
	require.NoError(t, err)

	_, err = code.DecodeFromSlice(make([]uint64, 15))
	assert.ErrorIs(t, err, ErrMismatchedLengths)

	_, err = code.DecodeFromSlice(make([]uint64, 17))
	assert.ErrorIs(t, err, ErrMismatchedLengths)
}

// TestEvaluationPointsReturnsCopy: a Code holds its own points for its whole life, so a
// caller mutating what EvaluationPoints hands back must not be able to reach them.
func TestEvaluationPointsReturnsCopy(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for name, opt := range map[string]Option{
		"ntt":       nil,
		"pointwise": Pointwise(),
	} {
		t.Run(name, func(t *testing.T) {
			code, err := NewCode(f, 16, 4, opt)
			require.NoError(t, err)

			first := code.EvaluationPoints()
			original := slices.Clone(first)

			for i := range first {
				first[i] = 123456
			}

			assert.Equal(t, original, code.EvaluationPoints(),
				"mutating a returned slice must not reach the code's own points")

			// And the code still works after that mutation attempt.
			data := makeTestSlice(code.K())
			enc, err := code.EncodeToSlice(data)
			require.NoError(t, err)
			dec, err := code.DecodeFromSlice(enc)
			require.NoError(t, err)
			assert.Equal(t, data, dec)
		})
	}
}

// TestDecodeIsConcurrencySafe: a *Code is documented as safe to share, and shares g0,
// the poly ring's twiddle cache and the evaluator's point cache across goroutines.
// Meaningful under -race.
func TestDecodeIsConcurrencySafe(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	code, err := NewCode(f, 64, 16, RequireNTT())
	require.NoError(t, err)

	data := makeTestSlice(code.K())

	codeword, err := code.EncodeToSlice(data)
	require.NoError(t, err)

	// Corrupt up to the correction budget so every goroutine runs the full GCD path.
	corrupted := slices.Clone(codeword)
	for i := 0; i < code.MaxErrors(); i++ {
		corrupted[i] = f.Reduce(uint64(i) + 7777)
	}

	const goroutines = 8

	errs := make(chan error, goroutines)

	for range goroutines {
		go func() {
			for range 20 {
				decoded, err := code.DecodeFromSlice(corrupted)
				if err != nil {
					errs <- err

					return
				}

				if len(decoded) == 0 || decoded[0] != data[0] {
					errs <- errors.New("concurrent decode produced the wrong message")

					return
				}
			}

			errs <- nil
		}()
	}

	for range goroutines {
		assert.NoError(t, <-errs)
	}
}
