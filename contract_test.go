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

// TestNewCodeParametersRejectsBadSizes covers parameters that used to panic from inside
// Encode instead of being reported at construction time.
func TestNewCodeParametersRejectsBadSizes(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range []struct {
		name string
		eval EvaluationMap
		n, k int
		want error
	}{
		{"n not a power of two", NewNttEvaluator(f), 20, 4, ErrUnsupportedSize},
		{"n does not divide p-1", NewNttEvaluator(f), 1 << 20, 4, ErrUnsupportedSize},
		{"n too small for an NTT", NewNttEvaluator(f), 1, 1, ErrUnsupportedSize},
		{"n exceeds the field", NewSlowEvaluator(f), 1 << 20, 4, ErrUnsupportedSize},
		{"n smaller than k", NewNttEvaluator(f), 4, 16, ErrNSmallerThanK},
		{"zero k", NewNttEvaluator(f), 16, 0, ErrNonPositiveK},
		{"negative k", NewNttEvaluator(f), 16, -1, ErrNonPositiveK},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := NewCodeParameters(tc.eval, tc.n, tc.k)
			assert.ErrorIs(t, err, tc.want)
		})
	}
}

// TestNewCodeParametersAcceptsValidSizes guards against the validation being so strict
// it rejects the cases the library is meant to serve.
func TestNewCodeParametersAcceptsValidSizes(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	// 929 is the prime field PDF417 barcodes are defined over; n is not a power of two,
	// so only the slow evaluator can serve it.
	pdf417, err := field.NewPrimeField(929)
	require.NoError(t, err)

	for _, tc := range []struct {
		name string
		eval EvaluationMap
		n, k int
	}{
		{"ntt power of two", NewNttEvaluator(f), 16, 4},
		{"ntt at the field's limit", NewNttEvaluator(f), 1 << 16, 4},
		{"slow non power of two", NewSlowEvaluator(f), 18, 5},
		{"slow over GF(929)", NewSlowEvaluator(pdf417), 100, 60},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := NewCodeParameters(tc.eval, tc.n, tc.k)
			assert.NoError(t, err)
		})
	}
}

// TestDecodeDoesNotMutateInput pins the contract that decoding leaves the caller's data
// alone. Decode used to write zeros into the map for every erased point, and
// DecodeFromSlice used to reduce and transform the caller's slice in place.
func TestDecodeDoesNotMutateInput(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for _, tc := range []testCase{
		{NewSlowEvaluator(f), 18, 5},
		{NewNttEvaluator(f), 16, 4},
	} {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		require.NoError(t, err)

		code := NewCodeGao(prms)
		data := makeTestSlice(tc.k)

		t.Run("Decode/map", func(t *testing.T) {
			encoded, err := code.Encode(data)
			require.NoError(t, err)

			// Drop the maximum number of points, so the erasure path runs.
			xs := prms.EvaluationPoints(prms.N())
			for i := 0; i < prms.N()-prms.K(); i++ {
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

	prms, err := NewCodeParameters(NewNttEvaluator(f), 16, 4)
	require.NoError(t, err)

	code := NewCodeGao(prms)

	_, err = code.DecodeFromSlice(make([]uint64, 15))
	assert.ErrorIs(t, err, ErrMismatchedLengths)

	_, err = code.DecodeFromSlice(make([]uint64, 17))
	assert.ErrorIs(t, err, ErrMismatchedLengths)
}

// TestEvaluationPointsReturnsCopy: the points are cached and reused by every encode and
// decode, so handing out the internal slice let a caller corrupt the whole code.
func TestEvaluationPointsReturnsCopy(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	for name, eval := range map[string]EvaluationMap{
		"ntt":  NewNttEvaluator(f),
		"slow": NewSlowEvaluator(f),
	} {
		t.Run(name, func(t *testing.T) {
			const n = 16

			first := eval.EvaluationPoints(n)
			original := slices.Clone(first)

			for i := range first {
				first[i] = 123456
			}

			assert.Equal(t, original, eval.EvaluationPoints(n),
				"mutating a returned slice must not disturb the cache")
		})
	}
}

// TestDecodeIsConcurrencySafe: a *Code is documented as safe to share, and shares g0,
// the poly ring's twiddle cache and the evaluator's point cache across goroutines.
// Meaningful under -race.
func TestDecodeIsConcurrencySafe(t *testing.T) {
	f, err := field.NewPrimeField(65537)
	require.NoError(t, err)

	prms, err := NewCodeParameters(NewNttEvaluator(f), 64, 16)
	require.NoError(t, err)

	code := NewCodeGao(prms)
	data := makeTestSlice(prms.K())

	codeword, err := code.EncodeToSlice(data)
	require.NoError(t, err)

	// Corrupt up to the correction budget so every goroutine runs the full GCD path.
	corrupted := slices.Clone(codeword)
	for i := 0; i < prms.MaxErrors(); i++ {
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
