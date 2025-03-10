package gao

import (
	"math/rand"
	"testing"
	"time"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
)

func TestNoCorruptions(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(157)
	a.NoError(err)

	prms, err := NewCodeParameters(NewSlowEvaluator(f), 16, 4)
	a.NoError(err)

	gao := NewCodeGao(prms)

	encoded, err := gao.Encode([]uint64{1, 2, 3, 4})
	a.NoError(err)

	// no corruptions
	decoded, err := gao.Decode(encoded)
	a.NoError(err)

	a.Equal([]uint64{1, 2, 3, 4}, decoded)
}

func TestErasures(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(157)
	a.NoError(err)

	prms, err := NewCodeParameters(NewSlowEvaluator(f), 16, 4)
	a.NoError(err)

	gao := NewCodeGao(prms)

	encoded, err := gao.Encode([]uint64{1, 2, 3, 4})
	a.NoError(err)

	// add erasures
	shuffledXs := shuffle(prms.EvaluationPoints(prms.n))
	for i := 0; i < prms.MaxErrors(); i++ {
		delete(encoded, shuffledXs[i])
	}

	a.Greater(prms.N(), len(encoded))

	decoded, err := gao.Decode(encoded)
	a.NoError(err)

	a.Equal([]uint64{1, 2, 3, 4}, decoded)
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
	f, err := field.NewPrimeField(157)
	a.NoError(err)

	prms, err := NewCodeParameters(NewSlowEvaluator(f), 16, 4)
	a.NoError(err)

	gao := NewCodeGao(prms)

	encoded, err := gao.Encode([]uint64{1, 2, 3, 4})
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

	a.Equal([]uint64{1, 2, 3, 4}, decoded)
}
