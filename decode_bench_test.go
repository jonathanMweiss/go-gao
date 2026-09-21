// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
)

var benchShapes = []struct{ n, k int }{{256, 128}, {1024, 512}, {4096, 2048}}

func benchCode(b *testing.B, n, k int) *ByteCode {
	b.Helper()

	f, err := field.NewPrimeField(field.NTTFriendlyPrime)
	if err != nil {
		b.Fatal(err)
	}

	code, err := NewCode(f, n, k, RequireNTT())
	if err != nil {
		b.Fatal(err)
	}

	return code.Bytes()
}

func benchWord(b *testing.B, bc *ByteCode, seed int64) ([]byte, int) {
	b.Helper()

	rng := rand.New(rand.NewSource(seed))

	payload := make([]byte, bc.MaxBytes())
	rng.Read(payload)

	raw, err := bc.Encode(payload)
	if err != nil {
		b.Fatal(err)
	}

	return raw, len(raw) / bc.code.N()
}

// BenchmarkErasureDecode: positions known, set built once and reused.
func BenchmarkErasureDecode(b *testing.B) {
	for _, sz := range benchShapes {
		bc := benchCode(b, sz.n, sz.k)

		for _, frac := range []int{1, 4, 2} { // (n-k)/frac erasures
			s := (sz.n - sz.k) / frac

			b.Run(fmt.Sprintf("n=%d/erasures=%d", sz.n, s), func(b *testing.B) {
				raw, width := benchWord(b, bc, 1)

				ranges := make([]ByteRange, s)
				for i := range ranges {
					ranges[i] = ByteRange{Off: i * width, Len: width}
					for j := range width {
						raw[i*width+j] = 0xAA
					}
				}

				lost, err := bc.Erasures(ranges...)
				if err != nil {
					b.Fatal(err)
				}

				b.SetBytes(int64(bc.MaxBytes()))
				b.ResetTimer()

				for b.Loop() {
					if _, err := bc.Decode(raw, lost); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
	}
}

// BenchmarkErrorDecode: positions unknown, so the decoder has to find them.
func BenchmarkErrorDecode(b *testing.B) {
	for _, sz := range benchShapes {
		bc := benchCode(b, sz.n, sz.k)
		maxErr := (sz.n - sz.k) / 2

		for _, frac := range []int{8, 2, 1} {
			e := max(1, maxErr/frac)

			b.Run(fmt.Sprintf("n=%d/errors=%d", sz.n, e), func(b *testing.B) {
				raw, width := benchWord(b, bc, 2)

				rng := rand.New(rand.NewSource(3))
				for _, idx := range rng.Perm(sz.n)[:e] {
					raw[idx*width] ^= 0xFF
				}

				if _, err := bc.Decode(raw, ErasureSet{}); err != nil {
					b.Fatal(err)
				}

				b.SetBytes(int64(bc.MaxBytes()))
				b.ResetTimer()

				for b.Loop() {
					if _, err := bc.Decode(raw, ErasureSet{}); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
	}
}

// BenchmarkErasureSetConstruction is the cost the erasure numbers amortise.
func BenchmarkErasureSetConstruction(b *testing.B) {
	for _, sz := range benchShapes {
		bc := benchCode(b, sz.n, sz.k)
		s := (sz.n - sz.k) / 2

		b.Run(fmt.Sprintf("n=%d/erasures=%d", sz.n, s), func(b *testing.B) {
			_, width := benchWord(b, bc, 1)

			ranges := make([]ByteRange, s)
			for i := range ranges {
				ranges[i] = ByteRange{Off: i * width, Len: width}
			}

			b.ResetTimer()

			for b.Loop() {
				if _, err := bc.Erasures(ranges...); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
