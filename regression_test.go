// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
)

// TestLargeCodeWithFastGCD pins a case that once failed only at scale: with n=16384 and
// a tenth of the codeword corrupted, the half-GCD path returned an error where the
// classical one decoded. It is kept as a regression test because nothing smaller
// reproduces it -- FastPartialGCD only takes its recursive path on large inputs.
func TestLargeCodeWithFastGCD(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(144115188075593729)
	a.NoError(err)

	k := 8192
	n := 2 * k
	gao, err := NewCode(f, n, k, RequireNTT())
	a.NoError(err)

	// Create test data
	slc := make([]uint64, k)
	for i := 0; i < k; i++ {
		slc[i] = uint64(i + 1)
	}

	encoded, err := gao.Encode(slc)
	a.NoError(err)

	corrupted := make([]uint64, len(encoded))
	copy(corrupted, encoded)

	corruptions := n / 10
	if corruptions == 0 {
		corruptions = 1
	}
	t.Logf("num corruptions %d out of %d", corruptions, n)

	// Fixed seed: this test exists to reproduce one specific failure.
	corruptCodeword(f, rand.New(rand.NewSource(1337)), corrupted, corruptions)

	decoded, err := gao.Decode(corrupted)
	if err != nil {
		t.Fatalf("fast decoding failed on %d corruptions out of %d, within the budget of %d: %v",
			corruptions, n, gao.MaxErrors(), err)
	}

	a.Equal(slc, decoded)
}
