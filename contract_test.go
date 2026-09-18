// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"fmt"
	"math/rand"
	"slices"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestNewCodeRejectsBadSizes covers parameters that used to panic from inside Encode
// instead of being reported at construction time.
func TestNewCodeRejectsBadSizes(t *testing.T) {
	for _, tc := range []struct {
		name string
		// prime overrides the default. Two cases below are about a size the field
		// cannot serve, which needs a field whose limits are reachable: over
		// NTTFriendlyPrime, n = 2^20 is both well inside the modulus and a divisor of
		// p-1, so neither rejection would fire and the cases would pass vacuously.
		prime uint64
		opts  []Option
		n, k  int
		want  error
	}{
		{"n smaller than k", 0, nil, 4, 16, ErrNSmallerThanK},
		{"zero k", 0, nil, 16, 0, ErrNonPositiveK},
		{"negative k", 0, nil, 16, -1, ErrNonPositiveK},
		{"n exceeds the field", 65537, nil, 1 << 20, 4, ErrUnsupportedSize},
		// RequireNTT turns the silent pointwise fallback into an error.
		{"require ntt, n not a power of two", 0, []Option{RequireNTT()}, 20, 4, ErrUnsupportedSize},
		{"require ntt, n does not divide p-1", 65537, []Option{RequireNTT()}, 1 << 20, 4, ErrUnsupportedSize},
		{"require ntt, n too small", 0, []Option{RequireNTT()}, 1, 1, ErrUnsupportedSize},
	} {
		t.Run(tc.name, func(t *testing.T) {
			prime := tc.prime
			if prime == 0 {
				prime = field.NTTFriendlyPrime
			}

			f, err := field.NewPrimeField(prime)
			require.NoError(t, err)

			_, err = NewCode(f, tc.n, tc.k, tc.opts...)
			assert.ErrorIs(t, err, tc.want)
		})
	}
}

func newfield(t testing.TB, prime uint64) *field.PrimeField {
	t.Helper()

	f, err := field.NewPrimeField(prime)
	if err != nil {

		t.Fatal(err)
	}

	return f
}
func ringOver(t testing.TB, prime uint64) *field.PolyRing {
	t.Helper()

	return field.NewPolyRing(newfield(t, prime))
}

// TestNewCodeSelectsStrategy: the evaluator is chosen for the caller, so the choice must
// be right — and visible, since falling back to the quadratic path is otherwise silent.
func TestNewCodeSelectsStrategy(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	// 929 is the prime field PDF417 barcodes are defined over. p-1 = 928 = 2^5 * 29, so
	// a transform exists only up to 32 points -- and since decoding needs a 2n-point one,
	// the usable n stops at 16. n=32 and n=100 must both fall back.
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
		{"ntt at the field's limit", f, nil, 1 << 15, 4, true},
		{"falls back when n is not a power of two", f, nil, 18, 5, false},
		{"Pointwise overrides an available ntt", f, []Option{Pointwise()}, 16, 4, false},
		{"GF(929) small enough for ntt", pdf417, nil, 16, 8, true},
		{"GF(929) evaluates at 32 but cannot decode there", pdf417, nil, 32, 8, false},
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
	f := newfield(t, 65537)

	_, err := NewCode(f, 20, 4, RequireNTT())
	require.Error(t, err)

	// The shape requirement and the field, so the reader can tell which to change.
	assert.Contains(t, err.Error(), "powers of two")
	assert.Contains(t, err.Error(), "65537")
}

// TestDecodeDoesNotMutateInput pins the contract that decoding leaves the caller's
// codeword alone: Decode used to reduce and transform the slice it was given in place.
// Both paths through it are covered, since the erasure path rescales every value before
// interpolating and so has the most to overwrite.
func TestDecodeDoesNotMutateInput(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	for _, tc := range []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	} {
		code, err := NewCode(f, tc.n, tc.k, tc.opt)
		require.NoError(t, err)

		data := makeTestSlice(tc.k)

		t.Run("erasure path", func(t *testing.T) {
			codeword, err := code.Encode(data)
			require.NoError(t, err)

			// Erase the maximum number of positions, so the rescaling path runs.
			erased := make([]int, code.N()-code.K())
			for i := range erased {
				erased[i] = i
				codeword[i] = 7777
			}

			before := slices.Clone(codeword)

			decoded, err := code.Decode(codeword, erased...)
			require.NoError(t, err)
			assert.Equal(t, data, decoded)

			assert.Equal(t, before, codeword, "Decode must not modify the slice it is given")
		})

		t.Run("error path", func(t *testing.T) {
			codeword, err := code.Encode(data)
			require.NoError(t, err)

			before := slices.Clone(codeword)

			decoded, err := code.Decode(codeword)
			require.NoError(t, err)
			assert.Equal(t, data, decoded)

			assert.Equal(t, before, codeword, "Decode must not modify the slice it is given")
		})
	}
}

// TestDecodeRejectsWrongLength: the positional API cannot express erasures, so
// a short slice is a caller error rather than a set of missing points.
func TestDecodeRejectsWrongLength(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	code, err := NewCode(f, 16, 4, RequireNTT())
	require.NoError(t, err)

	_, err = code.Decode(make([]uint64, 15))
	assert.ErrorIs(t, err, ErrMismatchedLengths)

	_, err = code.Decode(make([]uint64, 17))
	assert.ErrorIs(t, err, ErrMismatchedLengths)
}

// TestDecodeErasures: naming erasures is the whole point of the parameter —
// an erasure costs half an error, so a codeword that is hopeless when its damage is
// treated as errors decodes cleanly once the positions are declared.
func TestDecodeErasures(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	const n, k = 16, 4 // n-k = 12, so 6 errors or 12 erasures

	code, err := NewCode(f, n, k, RequireNTT())
	require.NoError(t, err)

	data := makeTestSlice(k)

	clean, err := code.Encode(data)
	require.NoError(t, err)

	t.Run("n-k erasures decode", func(t *testing.T) {
		damaged := slices.Clone(clean)

		erased := make([]int, 0, n-k)
		for i := range n - k {
			// Distinct garbage. Filling every slot with one value would place the
			// received word a mere n-k-(n-k)=4 errors from the constant codeword
			// p(x)=c, which is inside the budget of 6 and decodes legitimately.
			damaged[i] = uint64(1000 + i*7919)
			erased = append(erased, i)
		}

		// Declared as erasures: 12 <= n-k, so it decodes.
		decoded, err := code.Decode(damaged, erased...)
		require.NoError(t, err)
		assert.Equal(t, data, decoded)

		// Undeclared, the same slice is 12 errors against a budget of 6 — past the
		// distance bound, so the decoder either reports failure or lands on a
		// different codeword. What it must not do is return the original message.
		got, err := code.Decode(damaged)
		if err == nil {
			assert.NotEqual(t, data, got,
				"12 undeclared errors are beyond the correction radius")
		}
	})

	t.Run("value at an erased index is ignored", func(t *testing.T) {
		zeroed, garbage := slices.Clone(clean), slices.Clone(clean)

		erased := []int{1, 4, 9}
		for _, i := range erased {
			zeroed[i] = 0
			garbage[i] = 65536
		}

		fromZero, err := code.Decode(zeroed, erased...)
		require.NoError(t, err)

		fromGarbage, err := code.Decode(garbage, erased...)
		require.NoError(t, err)

		assert.Equal(t, data, fromZero)
		assert.Equal(t, fromZero, fromGarbage, "the filler value must not affect the result")
	})

	t.Run("errors and erasures together", func(t *testing.T) {
		// 4 erasures + 4 errors: 2*4 + 4 = 12 <= n-k.
		damaged := slices.Clone(clean)

		erased := []int{0, 1, 2, 3}
		for _, i := range erased {
			damaged[i] = 111
		}

		for _, i := range []int{7, 9, 11, 13} {
			damaged[i] = 222 // undeclared: genuine errors
		}

		decoded, err := code.Decode(damaged, erased...)
		require.NoError(t, err)
		assert.Equal(t, data, decoded)
	})

	t.Run("an erased position ignores whatever it holds", func(t *testing.T) {
		erased := []int{2, 5}

		zeroed := slices.Clone(clean)
		garbage := slices.Clone(clean)

		for _, i := range erased {
			zeroed[i] = 0
			garbage[i] = 7777
		}

		fromZeroed, err := code.Decode(zeroed, erased...)
		require.NoError(t, err)

		fromGarbage, err := code.Decode(garbage, erased...)
		require.NoError(t, err)

		assert.Equal(t, fromZeroed, fromGarbage)
		assert.Equal(t, data, fromGarbage)
	})
}

// TestDecodeRejectsBadErasures: these indices come straight from the caller,
// unlike the map form where absent keys are well-formed by construction. A duplicate
// would give the erasure locator a repeated root and corrupt the decode silently.
func TestDecodeRejectsBadErasures(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	const n, k = 16, 4

	code, err := NewCode(f, n, k, RequireNTT())
	require.NoError(t, err)

	ys := make([]uint64, n)

	for _, tc := range []struct {
		name   string
		erased []int
		want   error
	}{
		{"negative index", []int{-1}, ErrErasureOutOfRange},
		{"index == n", []int{n}, ErrErasureOutOfRange},
		{"index beyond n", []int{999}, ErrErasureOutOfRange},
		{"duplicate", []int{3, 3}, ErrDuplicateErasure},
		{"duplicate among valid", []int{1, 5, 1}, ErrDuplicateErasure},
		// Distinct indices, so this exercises the budget rather than the duplicate
		// check. n-k+1 erasures leave fewer than k good points.
		{"more than n-k", distinctIndices(n - k + 1), ErrTooManyMissingPoints},
		// Repeats must be reported as repeats even when there are enough of them to
		// also breach the budget, since the count is meaningless until they are gone.
		{"too many, but all duplicates", make([]int, n-k+1), ErrDuplicateErasure},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := code.Decode(ys, tc.erased...)
			assert.ErrorIs(t, err, tc.want)
		})
	}
}

// distinctIndices returns 0, 1, ..., count-1.
func distinctIndices(count int) []int {
	out := make([]int, count)
	for i := range out {
		out[i] = i
	}

	return out
}

// TestEvaluationPointsReturnsCopy: a Code holds its own points for its whole life, so a
// caller mutating what EvaluationPoints hands back must not be able to reach them.
func TestEvaluationPointsReturnsCopy(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

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
			enc, err := code.Encode(data)
			require.NoError(t, err)
			dec, err := code.Decode(enc)
			require.NoError(t, err)
			assert.Equal(t, data, dec)
		})
	}
}

// TestDecodeIsConcurrencySafe: a *Code is documented as safe to share, and shares g0,
// the poly ring's twiddle cache and the evaluator's point cache across goroutines.
// Meaningful under -race.
func TestDecodeIsConcurrencySafe(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	code, err := NewCode(f, 64, 16, RequireNTT())
	require.NoError(t, err)

	data := makeTestSlice(code.K())

	codeword, err := code.Encode(data)
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
				decoded, err := code.Decode(corrupted)
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

// TestDecodesAtFullErrorBudget: MaxErrors is a promise, and the interesting point is the
// boundary itself -- one corruption fewer exercises a different, slacker path through
// the partial GCD.
//
// It sweeps past hgcdThreshold (256 in the field package) deliberately. Below it the
// half-GCD falls back to the classical algorithm; above it the recursion runs, and a
// one-step overshoot there used to cost exactly one error of capability at every size
// from 256 up, while every test that stayed under the budget kept passing.
func TestDecodesAtFullErrorBudget(t *testing.T) {
	f := newfield(t, field.NTTFriendlyPrime)

	for _, n := range []int{16, 64, 128, 256, 512, 1024, 2048} {
		k := n / 2

		code, err := NewCode(f, n, k, RequireNTT())
		require.NoError(t, err)

		data := makeTestSlice(k)
		rng := rand.New(rand.NewSource(int64(n)))

		t.Run(fmt.Sprintf("errors/n=%d", n), func(t *testing.T) {
			codeword, err := code.Encode(data)
			require.NoError(t, err)

			corruptCodeword(f, rng, codeword, code.MaxErrors())

			decoded, err := code.Decode(codeword)
			require.NoError(t, err, "must correct exactly MaxErrors=%d corruptions", code.MaxErrors())
			require.Equal(t, data, decoded)
		})

		// The mixed budget 2*errors+erasures = n-k, spent at its other corner.
		t.Run(fmt.Sprintf("mixed/n=%d", n), func(t *testing.T) {
			codeword, err := code.Encode(data)
			require.NoError(t, err)

			erasures := (n - k) / 2
			errs := (n - k - erasures) / 2

			erasedAt := damageCodeword(f, rng, codeword, errs, erasures)

			decoded, err := code.Decode(codeword, erasedAt...)
			require.NoError(t, err, "must handle %d errors + %d erasures", errs, erasures)
			require.Equal(t, data, decoded)
		})
	}
}

// TestRequireNTTRejectsHalfFastFields: a field can admit an n-point transform and not a
// 2n-point one, which is enough to evaluate quickly and not enough to decode quickly --
// the partial GCD asks for convolutions of 1.25n, rounding up to 2n. The NTT strategy
// therefore requires both, and the half-equipped case falls back like any other n it
// cannot serve.
//
// 929 is the cheap witness: p-1 = 928 = 2^5 * 29 evaluates at 32 points and decodes at
// 16. The same gap exists at p=65537 between n=65536 and n=32768, where exercising it
// would mean building a pointwise code over 65536 points.
func TestRequireNTTRejectsHalfFastFields(t *testing.T) {
	pdf417, err := field.NewPrimeField(929)
	require.NoError(t, err)

	t.Run("2n available", func(t *testing.T) {
		code, err := NewCode(pdf417, 16, 8, RequireNTT())
		require.NoError(t, err, "2n=32 divides p-1=928")
		assert.True(t, code.UsesNTT())
	})

	t.Run("only n available", func(t *testing.T) {
		_, err := NewCode(pdf417, 32, 8, RequireNTT())
		require.ErrorIs(t, err, ErrUnsupportedSize)

		// The message has to carry the numbers: the reader's next question is which knob
		// to turn, and seeing 2n against p-1 says it is the prime, not n.
		assert.Contains(t, err.Error(), "2n-point transform")
		assert.Contains(t, err.Error(), "2n=64")
		assert.Contains(t, err.Error(), "p-1=928")
	})

	t.Run("falls back without the option", func(t *testing.T) {
		code, err := NewCode(pdf417, 32, 8)
		require.NoError(t, err)
		assert.False(t, code.UsesNTT(), "no 2n transform, so not the NTT strategy")
	})

	// The same ceiling on the prime the README recommends. Asserted through RequireNTT
	// because that fails inside selectEvaluator, before anything is built: letting it
	// fall back would construct a pointwise code over 65536 points and cost seconds.
	t.Run("65537 stops at 32768, not 65536", func(t *testing.T) {
		f, err := field.NewPrimeField(65537)
		require.NoError(t, err)

		code, err := NewCode(f, 1<<15, 4, RequireNTT())
		require.NoError(t, err, "2n=65536 divides p-1")
		assert.True(t, code.UsesNTT())

		_, err = NewCode(f, 1<<16, 4, RequireNTT())
		require.ErrorIs(t, err, ErrUnsupportedSize)
	})
}
