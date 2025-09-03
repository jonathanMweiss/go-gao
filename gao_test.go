package gao

import (
	"fmt"
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

func BenchmarkDecode(b *testing.B) {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}

	ks := []int{1 << 9, 1 << 10, 1 << 12, 1 << 13}

	evaluators := []struct {
		name string
		eval EvaluationMap
	}{
		{"slow", NewSlowEvaluator(f)},
		{"ntt", NewNttEvaluator(f)},
	}

	for _, k := range ks {
		k := k // capture
		for _, ev := range evaluators {
			ev := ev // capture
			n := k * 4
			name := fmt.Sprintf("eval=%s/n=%d/k=%d", ev.name, n, k)
			b.Run(name, func(b *testing.B) {
				// --- Setup (not timed) ---

				prms, err := NewCodeParameters(ev.eval, n, k)
				if err != nil {
					b.Fatal(err)
				}

				gao := NewCodeGao(prms)

				slc := makeTestSlice(k)

				encoding, err := gao.Encode(slc)
				if err != nil {
					b.Fatal(err)
				}

				// If Decode mutates the input slice, uncomment to protect the source:
				// mkCopy := func(src []Elem) []Elem { dst := make([]Elem, len(src)); copy(dst, src); return dst }

				// Rough throughput metric (bytes per op) if Elem is a byte-like type.
				// Adjust if your element size differs.
				b.SetBytes(int64(len(encoding)))
				b.ReportAllocs()
				b.ResetTimer()

				for i := 0; i < b.N; i++ {
					// enc := mkCopy(encoding) // use if Decode modifies input
					if _, err := gao.Decode(encoding); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
	}
}
