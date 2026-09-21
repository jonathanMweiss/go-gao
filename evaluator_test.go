// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/jonathanmweiss/go-gao/field"
)

// TestEvaluatorsRejectUnusableSizes pins where a bad n is caught: the constructor
// derives the points, so it has to report the sizes it cannot serve rather than
// panicking once someone asks for them.
func TestEvaluatorsRejectUnusableSizes(t *testing.T) {
	pr := field.NewPolyRing(newfield(t, 65537))

	for _, n := range []int{0, -1, 1} {
		_, err := newNttEvaluator(pr, n)
		assert.Error(t, err, "ntt n=%d", n)
	}

	// 65537-1 is 2^16, so transforms stop at 65536 and decoding needs 2n.
	_, err := newNttEvaluator(pr, 65536)
	assert.Error(t, err, "n beyond the 2n-point transform")

	for _, n := range []int{0, -1, 65537, 70000} {
		_, err := newSlowEvaluator(pr, n)
		assert.Error(t, err, "pointwise n=%d", n)
	}
}

// TestEvaluationPointsReturnsIndependentSlices pins that callers cannot reach the
// evaluator's own points, which are shared by every goroutine using the Code.
func TestEvaluationPointsReturnsIndependentSlices(t *testing.T) {
	const n = 32

	pr := field.NewPolyRing(newfield(t, field.NTTFriendlyPrime))

	ntt, err := newNttEvaluator(pr, n)
	require.NoError(t, err)

	slow, err := newSlowEvaluator(pr, n)
	require.NoError(t, err)

	for name, e := range map[string]evaluationMap{"ntt": ntt, "pointwise": slow} {
		t.Run(name, func(t *testing.T) {
			xs := e.EvaluationPoints()
			require.Len(t, xs, n)

			xs[0] = 12345
			assert.NotEqual(t, xs, e.EvaluationPoints())
		})
	}
}
