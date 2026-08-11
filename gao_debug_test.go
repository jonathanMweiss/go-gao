// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
)

func TestReproductionBenchmarkFailure(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(144115188075593729)
	a.NoError(err)

	k := 8192
	n := 2 * k
	eval := NewNttEvaluator(f)
	prms, err := NewCodeParameters(eval, n, k)
	a.NoError(err)

	gao := NewCodeGao(prms)

	// Create test data
	slc := make([]uint64, k)
	for i := 0; i < k; i++ {
		slc[i] = uint64(i + 1)
	}

	encoded, err := gao.EncodeToSlice(slc)
	a.NoError(err)

	corrupted := make([]uint64, len(encoded))
	copy(corrupted, encoded)

	corruptions := n / 10
	if corruptions == 0 {
		corruptions = 1
	}
	t.Logf("num corruptions %d out of %d", corruptions, n)

	rng := rand.New(rand.NewSource(1337))
	indices := rng.Perm(n)[:corruptions]
	for _, idx := range indices {
		corrupted[idx] = f.Reduce(uint64(rng.Uint32()))
	}

	// 1. Try decoding with FastPartialGCD (default)
	decoded, err := gao.DecodeFromSlice(corrupted)
	if err != nil {
		t.Fatalf("Reproduction successful: Fast decoding failed but it should have passed.")
	}

	a.Equal(slc, decoded)
}
