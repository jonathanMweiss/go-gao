package field

import (
	"fmt"
	"testing"
)

// blockedShapes covers both sides of blockedConvPlan's decision, so the agreement test
// exercises the blocked path and the plan test pins down which shapes reach it.
var blockedShapes = []struct{ ls, ll int }{
	{17, 1024},
	{33, 512},
	{163, 512},
	{163, 1024},
	{163, 16384},
	{4096, 16384},
	{8192, 16384},
	{70, 160},
	{130, 300},
	{1024, 1024},
	{12000, 16384},
	// Either side of the margin in blockedConvPlan: 65x3969 is predicted at 0.90 of the
	// standard cost and measures 1% slower, 129x3073 is predicted at 0.84 and measures
	// 8% faster. A wider or narrower margin moves one of them to the wrong path.
	{65, 3969},
	{129, 3073},
}

// ringOver builds a ring over an explicit modulus, for the one test below whose subject
// is the transform ceiling itself and so cannot use a field with room to spare.
func ringOver(t testing.TB, prime uint64) *PolyRing {
	t.Helper()

	f, err := NewPrimeField(prime)
	if err != nil {
		t.Fatal(err)
	}

	return NewPolyRing(f)
}

func rampPoly(r *PolyRing, n, seed int) *Polynomial {
	c := make([]uint64, n)
	for i := range c {
		c[i] = uint64((i+1)*(7919+seed)) % NTTFriendlyPrime
	}

	return r.NewPolynomial(c, false)
}

// TestBlockedConvMatchesWholeTransform checks the overlap-add result against the single
// transform it replaces, for every shape regardless of which one the plan prefers.
func TestBlockedConvMatchesWholeTransform(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	for _, s := range blockedShapes {
		t.Run(fmt.Sprintf("%dx%d", s.ls, s.ll), func(t *testing.T) {
			short, long := rampPoly(r, s.ls, 1), rampPoly(r, s.ll, 2)

			want := r.newDst()
			r.mulTruncInto(want, short, long, s.ls+s.ll-1)

			got := r.newDst()
			r.mulBlockedInto(got, short, long, nextPow2(s.ls))

			if len(got.inner) != len(want.inner) {
				t.Fatalf("length: got %d, want %d", len(got.inner), len(want.inner))
			}
			for i, v := range want.inner {
				if got.inner[i] != v {
					t.Fatalf("coefficient %d: got %d, want %d", i, got.inner[i], v)
				}
			}
		})
	}
}

// TestBlockedConvAliasedDestination covers the documented case where the destination is
// one of the operands, which the blocked path must not overwrite while reading it.
func TestBlockedConvAliasedDestination(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	for _, s := range []struct{ ls, ll int }{{163, 16384}, {33, 512}} {
		t.Run(fmt.Sprintf("%dx%d", s.ls, s.ll), func(t *testing.T) {
			short, long := rampPoly(r, s.ls, 1), rampPoly(r, s.ll, 2)

			want := r.newDst()
			r.mulBlockedInto(want, short, long, nextPow2(s.ls))

			intoLong := long.Copy()
			r.mulBlockedInto(intoLong, short, intoLong, nextPow2(s.ls))

			intoShort := short.Copy()
			r.mulBlockedInto(intoShort, intoShort, long, nextPow2(s.ls))

			for _, got := range []*Polynomial{intoLong, intoShort} {
				for i, v := range want.inner {
					if got.inner[i] != v {
						t.Fatalf("coefficient %d: got %d, want %d", i, got.inner[i], v)
					}
				}
			}
		})
	}
}

// TestBlockedConvPlan records which shapes the dispatch sends down the blocked path. The
// expectations are the measured outcome: blocking wins from a length ratio of two
// upwards, and loses where the last block is mostly padding or there is only one block.
func TestBlockedConvPlan(t *testing.T) {
	want := map[[2]int]bool{
		{17, 1024}:     true,
		{33, 512}:      true,
		{163, 512}:     true,
		{163, 1024}:    true,
		{163, 16384}:   true,
		{4096, 16384}:  true,
		{8192, 16384}:  true,
		{70, 160}:      false,
		{130, 300}:     false,
		{1024, 1024}:   false,
		{12000, 16384}: false,
		{65, 3969}:     false,
		{129, 3073}:    true,
	}

	for _, s := range blockedShapes {
		blockLen, ok := blockedConvPlan(s.ls, s.ll)
		expected, known := want[[2]int{s.ls, s.ll}]
		if !known {
			t.Errorf("%dx%d: shape has no expectation recorded", s.ls, s.ll)

			continue
		}
		if ok != expected {
			t.Errorf("%dx%d: blocked=%v, want %v", s.ls, s.ll, ok, expected)
		}
		if ok && blockLen < s.ls {
			t.Errorf("%dx%d: block size %d is shorter than the short operand", s.ls, s.ll, blockLen)
		}
	}
}

// TestBlockedConvReachedThroughMul confirms the plan is wired into Mul, not just
// reachable on its own, by checking Mul against a schoolbook reference on a shape the
// plan accepts.
func TestBlockedConvReachedThroughMul(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	short, long := rampPoly(r, 163, 1), rampPoly(r, 2048, 2)
	if _, ok := blockedConvPlan(163, 2048); !ok {
		t.Fatal("shape no longer takes the blocked path; pick another for this test")
	}

	want := r.newDst()
	r.mulSchoolbook(short, long, want)

	got := r.newDst()
	r.Mul(short, long, got)

	if len(got.inner) != len(want.inner) {
		t.Fatalf("length: got %d, want %d", len(got.inner), len(want.inner))
	}
	for i, v := range want.inner {
		if got.inner[i] != v {
			t.Fatalf("coefficient %d: got %d, want %d", i, got.inner[i], v)
		}
	}
}

// TestAddSubDefineEveryCoefficient pins the property that lets Add and Sub skip zeroing
// the tail ensureLen exposes: between the overlap loop and the tail copy they write every
// coefficient of the destination. A destination arriving with stale contents must not be
// able to show them through.
func TestAddSubDefineEveryCoefficient(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	for _, s := range []struct{ la, lb int }{{1, 1}, {3, 3}, {5, 2}, {2, 5}, {64, 9}, {9, 64}} {
		t.Run(fmt.Sprintf("%dx%d", s.la, s.lb), func(t *testing.T) {
			a, b := rampPoly(r, s.la, 1), rampPoly(r, s.lb, 2)

			clean := r.newDst()
			r.Add(a, b, clean)

			// Same operation into a destination pre-loaded with a longer run of junk.
			dirty := rampPoly(r, max(s.la, s.lb)+8, 3)
			dirty.inner = dirty.inner[:0]
			r.Add(a, b, dirty)

			if !dirty.Equals(clean) {
				t.Errorf("Add: stale destination contents showed through")
			}

			cleanSub, dirtySub := r.newDst(), rampPoly(r, max(s.la, s.lb)+8, 4)
			dirtySub.inner = dirtySub.inner[:0]
			r.Sub(a, b, cleanSub)
			r.Sub(a, b, dirtySub)

			if !dirtySub.Equals(cleanSub) {
				t.Errorf("Sub: stale destination contents showed through")
			}
		})
	}
}

// TestMulBeyondTransformCeiling covers what the blocked arm's feasibility check buys:
// p-1 is 2^16, so 65536 evaluation points is this field's ceiling and a product of 70162
// coefficients is past it. The blocks are 512 points, which fits, so blocking applies
// where a transform of the whole product cannot.
func TestMulBeyondTransformCeiling(t *testing.T) {
	r := ringOver(t, 65537)

	const ls, ll = 163, 70000
	if r.canUseNTTConvolutionLen(ls + ll - 1) {
		t.Fatal("shape no longer exceeds the ceiling; pick a longer one")
	}

	blockLen, lopsided := blockedConvPlan(ls, ll)
	if !lopsided || !r.canUseNTTConvolutionLen(2*blockLen) {
		t.Fatal("shape no longer reaches the blocked arm; the test proves nothing")
	}

	short, long := rampPoly(r, ls, 1), rampPoly(r, ll, 2)

	want := r.newDst()
	r.mulSchoolbook(short, long, want)

	got := r.newDst()
	r.Mul(short, long, got)

	if !got.Equals(want) {
		t.Fatalf("length %d vs %d", len(got.inner), len(want.inner))
	}
}
