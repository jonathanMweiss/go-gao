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

// corruptCodeword overwrites count distinct, randomly chosen positions of an encoded
// codeword with random field elements (count must be <= len(codeword)).
func corruptCodeword(f field.Field, rng *rand.Rand, codeword []uint64, count int) {
	for _, idx := range rng.Perm(len(codeword))[:count] {
		codeword[idx] = f.Reduce(rng.Uint64())
	}
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

		// add erasures. We should be able to handle up to n-k erasures.
		numErasures := prms.N() - prms.K()
		shuffledXs := shuffle(t, prms.EvaluationPoints(prms.n))
		for i := 0; i < numErasures; i++ {
			delete(encoded, shuffledXs[i])
		}

		a.Equal(prms.K(), len(encoded))

		decoded, err := gao.Decode(encoded)
		a.NoError(err)

		a.Equal(makeTestSlice(tc.k), decoded)
	}
}

func TestMixedErasuresAndCorruptions(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5}, // n-k=13. 2t+e <= 13. e=5, t=4 => 5+8=13.
		{NewNttEvaluator(f), 16, 4},  // n-k=12. 2t+e <= 12. e=4, t=4 => 4+8=12.
	}

	for _, tc := range testCases {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)

		gao := NewCodeGao(prms)
		originalData := makeTestSlice(tc.k)

		encoded, err := gao.Encode(originalData)
		a.NoError(err)

		xs := tc.EvaluationPoints(tc.n)
		shuffledXs := shuffle(t, xs)

		numErasures := 4
		if tc.n == 18 {
			numErasures = 5
		}
		numCorruptions := 4

		// Add erasures
		for i := 0; i < numErasures; i++ {
			delete(encoded, shuffledXs[i])
		}

		// Add corruptions
		for i := numErasures; i < numErasures+numCorruptions; i++ {
			encoded[shuffledXs[i]] = rand.Uint64()
		}

		decoded, err := gao.Decode(encoded)
		a.NoError(err)
		a.Equal(originalData, decoded)
	}
}

func shuffle(t *testing.T, slc []uint64) []uint64 {
	seed := time.Now().UnixNano()
	t.Logf("Shuffling with seed: %d", seed)

	rnd := rand.New(rand.NewSource(seed))

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
		shuffledXs := shuffle(t, prms.EvaluationPoints(prms.n))
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

func TestSliceEncodeDecode(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	rng := rand.New(rand.NewSource(1337))

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5},
		{NewNttEvaluator(f), 16, 4},
	}

	for _, tc := range testCases {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)

		gao := NewCodeGao(prms)
		originalData := makeTestSlice(tc.k)

		// Test EncodeToSlice and DecodeFromSlice with no corruptions
		encodedSlice, err := gao.EncodeToSlice(originalData)
		a.NoError(err)
		a.Len(encodedSlice, tc.n)

		encodedCopy := make([]uint64, len(encodedSlice))
		copy(encodedCopy, encodedSlice)

		decodedSlice, err := gao.DecodeFromSlice(encodedCopy)
		a.NoError(err)
		a.Equal(originalData, decodedSlice)

		// Test with corruptions
		corruptedSlice := make([]uint64, len(encodedSlice))
		copy(corruptedSlice, encodedSlice)
		corruptCodeword(f, rng, corruptedSlice, prms.MaxErrors())

		decodedFromCorrupted, err := gao.DecodeFromSlice(corruptedSlice)
		a.NoError(err)
		a.Equal(originalData, decodedFromCorrupted)
	}
}

// TestOptimisticErrorFreePath exercises the optimistic no-error short-circuit added to
// decodeNTT/decodeGeneric.
//
// (a) Decode is correct across the whole tolerated sweep 0..MaxErrors: 0 errors takes
// the fast path, 1..MaxErrors fall back to FastPartialGCD, and every one recovers the
// original message — so the fast path never misfires within tolerance.
//
// (b) sum of codewords is a codeword, thus fires the fast path.
func TestOptimisticErrorFreePath(t *testing.T) {
	a := assert.New(t)
	f, err := field.NewPrimeField(65537)
	a.NoError(err)

	rng := rand.New(rand.NewSource(1337))

	testCases := []testCase{
		{NewSlowEvaluator(f), 18, 5}, // generic (interpolation) path
		{NewNttEvaluator(f), 16, 4},  // NTT path
	}

	for _, tc := range testCases {
		prms, err := NewCodeParameters(tc.EvaluationMap, tc.n, tc.k)
		a.NoError(err)
		gao := NewCodeGao(prms)

		msg := makeTestSlice(tc.k)
		enc1, err := gao.EncodeToSlice(msg)
		a.NoError(err)

		// (a) correctness across 0..MaxErrors corruptions at random positions.
		for e := 0; e <= prms.MaxErrors(); e++ {
			work := make([]uint64, len(enc1))
			copy(work, enc1)

			corruptCodeword(f, rng, work, e)

			decoded, err := gao.DecodeFromSlice(work)
			a.NoError(err, "n=%d errors=%d", tc.n, e)
			a.Equal(msg, decoded, "n=%d errors=%d", tc.n, e)
		}

		// (b) sum of codewords is a codeword.
		msg2 := make([]uint64, tc.k)
		for i := range msg2 {
			msg2[i] = uint64(2*i + 1)
		}
		enc2, err := gao.EncodeToSlice(msg2)
		a.NoError(err)

		summed := make([]uint64, tc.n)
		for i := range summed {
			summed[i] = f.Add(enc1[i], enc2[i])
		}

		decoded, err := gao.DecodeFromSlice(summed)
		a.NoError(err, "beyond-tolerance n=%d", tc.n)

		want := make([]uint64, tc.k)
		for i := range want {
			want[i] = f.Add(msg[i], msg2[i])
		}
		a.Equal(want, decoded, "beyond-tolerance n=%d: decode returns the sum codeword's message", tc.n)
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

func BenchmarkDecodeFromSliceOnePercentCorruptionsNTT(b *testing.B) {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		b.Fatal(err)
	}

	ks := []int{1 << 12, 1 << 13, 1 << 14, 1 << 15}
	rng := rand.New(rand.NewSource(1337))

	for _, k := range ks {
		n := 2 * k
		prms, err := NewCodeParameters(NewNttEvaluator(f), n, k)
		if err != nil {
			b.Fatal(err)
		}

		gao := NewCodeGao(prms)
		slc := makeTestSlice(k)

		encoded, err := gao.EncodeToSlice(slc)
		if err != nil {
			b.Fatal(err)
		}

		corrupted := make([]uint64, len(encoded))
		copy(corrupted, encoded)

		corruptions := n / 100
		if corruptions == 0 {
			corruptions = 1
		}

		corruptCodeword(f, rng, corrupted, corruptions)

		name := fmt.Sprintf("n=%d/k=%d/errors=%d(1%%)", n, k, corruptions)
		b.Run(name, func(b *testing.B) {
			work := make([]uint64, len(corrupted))

			b.SetBytes(int64(len(corrupted) * 8))
			b.ReportAllocs()
			b.ResetTimer()

			for i := 0; i < b.N; i++ {
				copy(work, corrupted)
				if _, err := gao.DecodeFromSlice(work); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
