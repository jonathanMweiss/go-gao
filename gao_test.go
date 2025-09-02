package gao

import (
	"math/rand"
	"testing"
	"time"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
)

type testCase struct {
	EvaluationMap
	n, k int
}

func makeTestSlice(k int) []uint64 {
	poly := make([]uint64, k)
	for i := 0; i < k; i++ {
		poly[i] = uint64(i + 1)
	}

	return poly
}

func TestNoCorruptions(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5},
		{NewNttEvaluator(f), 16, 4}, // checking non powers of 2.
	}

	for _, tc := range testCases {

		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)

		gao := NewCodeGao(prms)

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
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5},
		{NewNttEvaluator(f), 16, 4}, // checking non powers of 2.
	}

	for _, tc := range testCases {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)

		gao := NewCodeGao(prms)

		encoded, err := gao.Encode(makeTestSlice(tc.k))
		a.NoError(err)

		// add erasures
		shuffledXs := shuffle(prms.EvaluationPoints(prms.n))
		for i := 0; i < prms.MaxErrors(); i++ {
			delete(encoded, shuffledXs[i])
		}

		a.Greater(prms.N(), len(encoded))

		decoded, err := gao.Decode(encoded)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}
}

func shuffle(slc []uint64) []uint64 {
	rnd := rand.New(rand.NewSource(time.Now().Unix()))

	cpy := make([]uint64, len(slc))
	copy(cpy, slc)

	rnd.Shuffle(len(slc), func(i, j int) {
		cpy[i], cpy[j] = cpy[j], cpy[i]
	})

	return cpy
}

func TestCorruptions(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5},
		{NewNttEvaluator(f), 16, 4}, // checking non powers of 2.
	}

	for _, tc := range testCases {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)

		gao := NewCodeGao(prms)

		encoded, err := gao.Encode(makeTestSlice(tc.k))
		a.NoError(err)

		corrupted := make(map[uint64]uint64, len(encoded))
		for x, y := range encoded {
			corrupted[x] = y
		}

		// add corruptions
		shuffledXs := shuffle(prms.EvaluationPoints(prms.n))
		for i := 0; i < prms.MaxErrors(); i++ {
			corrupted[shuffledXs[i]] = rand.Uint64()
		}

		a.Len(corrupted, prms.N())
		a.NotEqual(encoded, corrupted)

		decoded, err := gao.Decode(corrupted)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}
}
