// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
)

// BenchmarkErasureSetReuse decodes a batch of words that all lost the same positions,
// once with a set built for the batch and once with a set per word. The gap is the work
// Erasures hoists out of the per-word path.
func BenchmarkErasureSetReuse(b *testing.B) {
	const batch = 64

	for _, sz := range []struct{ n, k int }{{256, 128}, {2048, 1024}, {8192, 4096}} {
		f, err := field.NewPrimeField(field.NTTFriendlyPrime)
		if err != nil {
			b.Fatal(err)
		}

		code, err := NewCode(f, sz.n, sz.k, RequireNTT())
		if err != nil {
			b.Fatal(err)
		}

		erased := make([]int, sz.n/8)
		for i := range erased {
			erased[i] = i * 3
		}

		words := damagedWords(b, code, erased, batch)

		b.Run(fmt.Sprintf("n=%d/one set for the batch", sz.n), func(b *testing.B) {
			es, err := code.Erasures(erased...)
			if err != nil {
				b.Fatal(err)
			}

			b.ResetTimer()

			for b.Loop() {
				for _, w := range words {
					if _, err := code.Decode(w, es); err != nil {
						b.Fatal(err)
					}
				}
			}
		})

		b.Run(fmt.Sprintf("n=%d/one set per word", sz.n), func(b *testing.B) {
			for b.Loop() {
				for _, w := range words {
					es, err := code.Erasures(erased...)
					if err != nil {
						b.Fatal(err)
					}

					if _, err := code.Decode(w, es); err != nil {
						b.Fatal(err)
					}
				}
			}
		})
	}
}

// damagedWords encodes count random messages and overwrites the erased positions.
func damagedWords(b *testing.B, code *Code, erased []int, count int) []Codeword {
	b.Helper()

	rng := rand.New(rand.NewSource(1))
	f := code.PrimeField()
	out := make([]Codeword, count)

	for w := range out {
		msg := make([]uint64, code.K())
		for i := range msg {
			msg[i] = f.Reduce(rng.Uint64())
		}

		word, err := code.Encode(msg)
		if err != nil {
			b.Fatal(err)
		}

		for _, i := range erased {
			word[i] = differentElement(f, rng, word[i])
		}

		out[w] = word
	}

	return out
}
