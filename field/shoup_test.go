package field

import "testing"

// TestShoupMatchesFieldMul checks [PrimeField.mulShoup] against [fieldMul] over complete
// input spaces.
//
// The small primes are exhaustive, which is the strongest statement available: every
// (w, x) pair agrees. None of them reaches the conditional subtraction inside mulShoup,
// though. That fires only when the quotient can come out one short, which needs
// p*p > 2^64, so p > 2^32. FuzzShoupMatchesFieldMul carries that path.
func TestShoupMatchesFieldMul(t *testing.T) {
	for _, prime := range []uint64{2, 3, 929, 3329} {
		f, err := NewPrimeField(prime)
		if err != nil {
			t.Fatalf("prime %d: %v", prime, err)
		}

		for w := uint64(0); w < prime; w++ {
			wp := f.shoupFactor(w)
			for x := uint64(0); x < prime; x++ {
				if got, want := f.mulShoup(w, wp, x), f.Mul(w, x); got != want {
					t.Fatalf("p=%d w=%d x=%d: mulShoup=%d fieldMul=%d", prime, w, x, got, want)
				}
			}
		}
	}

	// The sizes actually in use, where w*x overflows 64 bits and the low-word
	// subtraction has to rely on the wrapping cancelling.
	for _, prime := range []uint64{65537, largePrime} {
		f, err := NewPrimeField(prime)
		if err != nil {
			t.Fatal(err)
		}

		edges := []uint64{0, 1, 2, prime / 2, prime - 2, prime - 1}
		for _, w := range edges {
			wp := f.shoupFactor(w)
			for _, x := range edges {
				if got, want := f.mulShoup(w, wp, x), f.Mul(w, x); got != want {
					t.Fatalf("p=%d w=%d x=%d: mulShoup=%d fieldMul=%d", prime, w, x, got, want)
				}
			}
		}
	}
}

// FuzzShoupMatchesFieldMul hunts for a (w, x) the two multiplications disagree on.
//
// The 63-bit prime is the interesting half: it is the only size here where mulShoup's
// quotient can come out one short, so it is the only one that exercises the conditional
// subtraction. 65537 never reaches it.
func FuzzShoupMatchesFieldMul(fz *testing.F) {
	// 1<<62 times 4 is exactly 2^64, so their product wraps to zero in a uint64 while
	// neither operand is zero -- the case that rules out testing a*b == 0.
	for _, tc := range [][2]uint64{
		{0, 0}, {1, 1}, {1, largePrime - 1}, {largePrime - 1, largePrime - 1},
		{54347, 4534523}, {1 << 62, 4}, {largePrime / 2, largePrime / 2},
	} {
		fz.Add(tc[0], tc[1], false)
		fz.Add(tc[0], tc[1], true)
	}

	small, err := NewPrimeField(65537)
	if err != nil {
		fz.FailNow()
	}

	large, err := NewPrimeField(largePrime)
	if err != nil {
		fz.FailNow()
	}

	fz.Fuzz(func(t *testing.T, w, x uint64, useSmall bool) {
		f := large
		if useSmall {
			f = small
		}

		p := f.Modulus()
		w %= p
		x %= p

		if got, want := f.mulShoup(w, f.shoupFactor(w), x), f.Mul(w, x); got != want {
			t.Fatalf("p=%d w=%d x=%d: mulShoup=%d fieldMul=%d", p, w, x, got, want)
		}
	})
}
