// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/jonathanmweiss/go-gao/field"
)

// TestEvaluationPointsAreDerivedOnce pins the memo: every call hands back the same
// slice, so an nttEvaluator runs its transform once rather than per call.
func TestEvaluationPointsAreDerivedOnce(t *testing.T) {
	const n = 32

	pr := field.NewPolyRing(newfield(t, field.NTTFriendlyPrime))

	ntt := newNttEvaluator(pr, n)
	assert.Same(t, &ntt.points()[0], &ntt.points()[0])

	slow := newSlowEvaluator(pr, n)
	assert.Same(t, &slow.points()[0], &slow.points()[0])

	// EvaluationPoints hands out a copy, so a caller cannot reach the shared slice.
	xs := ntt.EvaluationPoints()
	xs[0] = 12345
	assert.NotEqual(t, xs, ntt.EvaluationPoints())
}

// TestEvaluationPointsAreConcurrencySafe derives the points from many goroutines at
// once, which is the case a Code shared across goroutines can reach. Run under -race.
func TestEvaluationPointsAreConcurrencySafe(t *testing.T) {
	const n = 64

	pr := field.NewPolyRing(newfield(t, field.NTTFriendlyPrime))

	for name, e := range map[string]evaluationMap{
		"ntt":       newNttEvaluator(pr, n),
		"pointwise": newSlowEvaluator(pr, n),
	} {
		t.Run(name, func(t *testing.T) {
			var (
				wg  sync.WaitGroup
				mu  sync.Mutex
				got [][]uint64
			)

			for range 16 {
				wg.Add(1)

				go func() {
					defer wg.Done()

					xs := e.EvaluationPoints()

					mu.Lock()
					defer mu.Unlock()

					got = append(got, xs)
				}()
			}

			wg.Wait()
			require.Len(t, got, 16)

			for _, xs := range got {
				assert.Equal(t, got[0], xs)
			}
		})
	}
}
