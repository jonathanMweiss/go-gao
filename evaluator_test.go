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

// TestEvaluationPointsCacheIsConcurrencySafe drives a cold cache from many goroutines at
// once, which is the case a Code shared across goroutines can reach. Run under -race.
func TestEvaluationPointsCacheIsConcurrencySafe(t *testing.T) {
	const n = 64

	pr := field.NewPolyRing(newfield(t, field.NTTFriendlyPrime))

	for name, e := range map[string]evaluationMap{
		"ntt":       newNttEvaluator(pr),
		"pointwise": newSlowEvaluator(pr),
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

					xs := e.EvaluationPoints(n)

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

			// each caller owns its slice: writing one must not disturb the next call.
			got[0][0] = 12345
			assert.NotEqual(t, got[0], e.EvaluationPoints(n))
		})
	}
}

// TestEvaluationPointsAreDerivedOnce pins the cache itself: every call hands back the
// same slice, so an nttEvaluator runs its transform once per length rather than per call.
func TestEvaluationPointsAreDerivedOnce(t *testing.T) {
	const n = 32

	pr := field.NewPolyRing(newfield(t, field.NTTFriendlyPrime))

	ntt := newNttEvaluator(pr)
	assert.Same(t, &ntt.cachedPoints(n)[0], &ntt.cachedPoints(n)[0])

	slow := newSlowEvaluator(pr)
	assert.Same(t, &slow.cachedPoints(n)[0], &slow.cachedPoints(n)[0])

	// a different length gets its own entry.
	assert.NotEqual(t, len(ntt.cachedPoints(n)), len(ntt.cachedPoints(2*n)))
}
