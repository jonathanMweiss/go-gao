// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"math/rand"
	"runtime"
	"slices"
	"testing"
	"time"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
)

type testCase struct {
	name string
	opt  Option
	n, k int
}

func makeTestSlice(k int) []uint64 {
	poly := make([]uint64, k)
	for i := 0; i < k; i++ {
		poly[i] = uint64(i + 1)
	}

	return poly
}

// testRNG returns a randomly seeded generator, logging the seed so a failure can be
// replayed.
func testRNG(t *testing.T) *rand.Rand {
	seed := time.Now().UnixNano()
	t.Logf("random seed: %d", seed)

	return rand.New(rand.NewSource(seed))
}

// damageCodeword damages a clean codeword in place: it corrupts errs positions and
// declares erasures more, drawn from one permutation so that the two never overlap and
// the cost against the n-k budget really is 2*errs+erasures. It returns the erasure
// indices, ready to hand to Decode.
//
// Garbage is written at the erased positions too, which is what makes a passing decode
// evidence that a declared position's value is ignored rather than merely tolerated.
//
// errs+erasures must not exceed len(codeword).
func damageCodeword(f field.Field, rng *rand.Rand, codeword []uint64, errs, erasures int) []int {
	positions := rng.Perm(len(codeword))

	for _, idx := range positions[:errs+erasures] {
		codeword[idx] = differentElement(f, rng, codeword[idx])
	}

	// the positions of all the erasures are returned.
	return positions[errs : errs+erasures]
}

// differentElement returns a random field element that is not old.
// ensuring that the old value is replaced with a different one.
func differentElement(f field.Field, rng *rand.Rand, old uint64) uint64 {
	for {
		if v := f.Reduce(rng.Uint64()); v != old {
			return v
		}
	}
}

// corruptCodeword overwrites count distinct, randomly chosen positions of an encoded
// codeword with random field elements (count must be <= len(codeword)).
func corruptCodeword(f field.Field, rng *rand.Rand, codeword []uint64, count int) {
	damageCodeword(f, rng, codeword, count, 0)
}

func TestNoCorruptions(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	}

	for _, tc := range testCases {

		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)

		encoded, err := gao.Encode(makeTestSlice(tc.k))
		a.NoError(err)

		// no corruptions
		decoded, err := gao.Decode(encoded)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}

}

func TestErasures(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	}

	for _, tc := range testCases {
		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)

		encoded, err := gao.Encode(makeTestSlice(tc.k))
		a.NoError(err)

		// add erasures. We should be able to handle up to n-k erasures.
		erased := damageCodeword(f, testRNG(t), encoded, 0, gao.N()-gao.K())

		decoded, err := gao.Decode(encoded, erased...)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}
}

func TestMixedErasuresAndCorruptions(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5}, // n-k=13. 2t+e <= 13. e=5, t=4 => 5+8=13.
		{"ntt", nil, 16, 4},               // n-k=12. 2t+e <= 12. e=4, t=4 => 4+8=12.
	}

	for _, tc := range testCases {
		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)
		originalData := makeTestSlice(tc.k)

		encoded, err := gao.Encode(originalData)
		a.NoError(err)

		numErasures := 4
		if tc.n == 18 {
			numErasures = 5
		}
		numCorruptions := 4

		erased := damageCodeword(f, testRNG(t), encoded, numCorruptions, numErasures)

		decoded, err := gao.Decode(encoded, erased...)
		a.NoError(err)
		a.Equal(originalData, decoded)
	}
}

func TestCorruptions(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	}

	for _, tc := range testCases {
		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)

		encoded, err := gao.Encode(makeTestSlice(tc.k))
		a.NoError(err)

		corrupted := slices.Clone(encoded)

		// add corruptions
		corruptCodeword(f, testRNG(t), corrupted, gao.MaxErrors())

		a.Len(corrupted, gao.N())
		a.NotEqual(encoded, corrupted)

		decoded, err := gao.Decode(corrupted)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}
}

func TestSliceEncodeDecode(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	rng := rand.New(rand.NewSource(1337))

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5},
		{"ntt", nil, 16, 4},
	}

	for _, tc := range testCases {
		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)
		originalData := makeTestSlice(tc.k)

		// Test Encode and Decode with no corruptions
		encodedSlice, err := gao.Encode(originalData)
		a.NoError(err)
		a.Len(encodedSlice, tc.n)

		encodedCopy := make([]uint64, len(encodedSlice))
		copy(encodedCopy, encodedSlice)

		decodedSlice, err := gao.Decode(encodedCopy)
		a.NoError(err)
		a.Equal(originalData, decodedSlice)

		// Test with corruptions
		corruptedSlice := make([]uint64, len(encodedSlice))
		copy(corruptedSlice, encodedSlice)
		corruptCodeword(f, rng, corruptedSlice, gao.MaxErrors())

		decodedFromCorrupted, err := gao.Decode(corruptedSlice)
		a.NoError(err)
		a.Equal(originalData, decodedFromCorrupted)
	}
}

// TestOptimisticErrorFreePath exercises the optimistic no-error short-circuit added to
// decodeNTT/decodeGeneric.
//
// (a) Decode is correct across the whole tolerated sweep 0..MaxErrors: 0 errors takes
// the fast path, 1..MaxErrors fall back to PartialGCD, and every one recovers the
// original message — so the fast path never misfires within tolerance.
//
// (b) sum of codewords is a codeword, thus fires the fast path.
func TestOptimisticErrorFreePath(t *testing.T) {
	a := assert.New(t)
	f := newfield(t, field.NTTFriendlyPrime)

	rng := rand.New(rand.NewSource(1337))

	testCases := []testCase{
		{"pointwise", Pointwise(), 18, 5}, // generic (interpolation) path
		{"ntt", nil, 16, 4},               // NTT path
	}

	for _, tc := range testCases {
		gao, err := NewCode(f, tc.n, tc.k, tc.opt)
		a.NoError(err)

		msg := makeTestSlice(tc.k)
		enc1, err := gao.Encode(msg)
		a.NoError(err)

		// (a) correctness across 0..MaxErrors corruptions at random positions.
		for e := 0; e <= gao.MaxErrors(); e++ {
			work := make([]uint64, len(enc1))
			copy(work, enc1)

			corruptCodeword(f, rng, work, e)

			decoded, err := gao.Decode(work)
			a.NoError(err, "n=%d errors=%d", tc.n, e)
			a.Equal(msg, decoded, "n=%d errors=%d", tc.n, e)
		}

		// (b) sum of codewords is a codeword.
		msg2 := make([]uint64, tc.k)
		for i := range msg2 {
			msg2[i] = uint64(2*i + 1)
		}
		enc2, err := gao.Encode(msg2)
		a.NoError(err)

		summed := make([]uint64, tc.n)
		for i := range summed {
			summed[i] = f.Add(enc1[i], enc2[i])
		}

		decoded, err := gao.Decode(summed)
		a.NoError(err, "beyond-tolerance n=%d", tc.n)

		want := make([]uint64, tc.k)
		for i := range want {
			want[i] = f.Add(msg[i], msg2[i])
		}
		a.Equal(want, decoded, "beyond-tolerance n=%d: decode returns the sum codeword's message", tc.n)
	}
}

// BenchmarkDecode measures decoding across both evaluation strategies.
//
// The errors dimension is not decoration. With a clean codeword g1 has degree < k and
// Decode returns on the optimistic path without ever running the partial GCD, so the
// NTT arm measures an inverse transform and little else -- 9 allocations for the whole
// call. Only a corrupted codeword exercises the decoder.
//
// Pointwise evaluation is quadratic and gets a shorter size list: at n=32768 one
// error-free decode already runs 45s and allocates 26 GB, which is a number to know
// rather than one to re-measure on every run.
func BenchmarkDecode(b *testing.B) {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}

	evaluators := []struct {
		name string
		opt  Option
		ks   []int
	}{
		{"pointwise", Pointwise(), []int{1 << 9, 1 << 10, 1 << 11}},
		{"ntt", RequireNTT(), []int{1 << 9, 1 << 10, 1 << 12, 1 << 13}},
	}

	for _, ev := range evaluators {
		for _, k := range ev.ks {
			n := k * 4

			for _, errs := range []int{0, max(1, n/100)} {
				name := fmt.Sprintf("eval=%s/n=%d/k=%d/errors=%d", ev.name, n, k, errs)

				b.Run(name, func(b *testing.B) {
					// --- Setup (not timed) ---
					gao, err := NewCode(f, n, k, ev.opt)
					if err != nil {
						b.Fatal(err)
					}

					encoding, err := gao.Encode(makeTestSlice(k))
					if err != nil {
						b.Fatal(err)
					}

					if errs > 0 {
						corruptCodeword(f, rand.New(rand.NewSource(1337)), encoding, errs)
					}

					// Symbols are uint64, so the codeword is 8 bytes per position.
					b.SetBytes(int64(len(encoding) * 8))
					b.ReportAllocs()
					b.ResetTimer()

					// Decode does not modify its input, so one codeword serves every
					// iteration.
					for i := 0; i < b.N; i++ {
						if _, err := gao.Decode(encoding); err != nil {
							b.Fatal(err)
						}
					}
				})
			}
		}
	}
}

// BenchmarkDecodeParallel is where allocation work either pays or does not.
//
// Decode allocates per call, so its cost to a server is not only its own runtime but the
// GC cycles that allocation rate drives, which are charged to every other goroutine in the
// process. A single-threaded ns/op cannot see that. This runs the same decode on every
// core at once, where the allocator and the collector are actually contended.
func BenchmarkDecodeParallel(b *testing.B) {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}

	const k = 1 << 12

	n := 2 * k

	code, err := NewCode(f, n, k, RequireNTT())
	if err != nil {
		b.Fatal(err)
	}

	encoded, err := code.Encode(makeTestSlice(k))
	if err != nil {
		b.Fatal(err)
	}

	corruptCodeword(f, rand.New(rand.NewSource(1337)), encoded, n/100)

	for _, procs := range []int{1, 4, runtime.NumCPU()} {
		b.Run(fmt.Sprintf("procs=%d", procs), func(b *testing.B) {
			b.SetParallelism(procs)
			b.SetBytes(int64(len(encoded) * 8))
			b.ReportAllocs()
			b.ResetTimer()

			b.RunParallel(func(pb *testing.PB) {
				for pb.Next() {
					if _, err := code.Decode(encoded); err != nil {
						b.Fatal(err)
					}
				}
			})
		})
	}
}

func BenchmarkDecodeOnePercentCorruptionsNTT(b *testing.B) {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}

	// n = 2k, and decoding needs a 2n-point transform, so over p=65537 (p-1 = 2^16) the
	// largest usable k is 2^14: k = 2^15 would ask for a 131072-point transform and
	// NewCode rejects it.
	ks := []int{1 << 11, 1 << 12, 1 << 13, 1 << 14}
	rng := rand.New(rand.NewSource(1337))

	for _, k := range ks {
		n := 2 * k
		gao, err := NewCode(f, n, k, RequireNTT())
		if err != nil {
			b.Fatal(err)
		}
		slc := makeTestSlice(k)

		encoded, err := gao.Encode(slc)
		if err != nil {
			b.Fatal(err)
		}

		corrupted := make([]uint64, len(encoded))
		copy(corrupted, encoded)

		corruptions := n / 100
		if corruptions == 0 {
			corruptions = 1
		}

		corruptCodeword(f, rng, corrupted, corruptions)

		name := fmt.Sprintf("n=%d/k=%d/errors=%d(1%%)", n, k, corruptions)
		b.Run(name, func(b *testing.B) {
			b.SetBytes(int64(len(corrupted) * 8))
			b.ReportAllocs()
			b.ResetTimer()

			// Decode does not modify its input, so the codeword is reused as is rather
			// than re-copied inside the timed loop.
			for i := 0; i < b.N; i++ {
				if _, err := gao.Decode(corrupted); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
