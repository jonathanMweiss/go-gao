// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestNTTForward(t *testing.T) {
	a := assert.New(t)
	f, err := NewPrimeField(3329)
	a.NoError(err)

	p := newPolynomial(f, []uint64{1, 2, 3, 4, 5, 6, 7, 8}, false)
	expected := []uint64{36, 3240, 3067, 427, 3325, 2894, 254, 81}

	pr := NewPolyRing(f)

	a.NoError(pr.NttForward(p))
	a.Equal(expected, p.ToSlice())
}

func TestNTTForwardBackward(t *testing.T) {
	// Test the forward and backward NTT transforms.
	a := assert.New(t)
	f := newPrimeField(t, NTTFriendlyPrime)

	pr := NewPolyRing(f)
	for i := range 8 {
		cappingDegree := 1 << (i + 1)

		p1 := randomPolynomial(f, 12345+uint64(i), cappingDegree)
		pcpy := p1.Copy()

		a.NoError(pr.NttForward(p1))

		a.NoError(pr.NttBackward(p1))

		a.True(pcpy.Equals(p1))
	}
}

func TestPolyMult(t *testing.T) {
	a := assert.New(t)
	f := newPrimeField(t, NTTFriendlyPrime)

	pr := NewPolyRing(f)

	for i := 0; i < 8; i++ {

		degree := 1 << (i + 1)

		p1 := randomPolynomial(f, 12345+uint64(i), degree)

		regMulRes := &Polynomial{}
		pr.mulSchoolbook(p1, p1, regMulRes)

		// padding p1 with zeros:
		p1.inner = append(p1.inner, make([]uint64, degree)...)

		a.NoError(pr.NttForward(p1))

		nttRes := &Polynomial{}
		pr.mulWholeNTT(p1, p1, nttRes)

		a.NoError(pr.NttBackward(nttRes))

		if !regMulRes.Equals(nttRes) {
			regMulRes.Equals(nttRes)
			t.Errorf("Mismatch between regular and NTT multiplication results")
		}
		a.True(regMulRes.Equals(nttRes))
	}
}

// TestNTTRecursiveMatchesIterative holds the in-place transforms to the recursive
// reference they were derived from. A reference nothing checks drifts from the code it
// is supposed to explain, so both directions are compared at every size the field
// admits. The two reach their roots by different routes -- the reference looks one up
// per node, the ring squares a cached psi down the stages -- so agreeing is a real
// check on both.
func TestNTTRecursiveMatchesIterative(t *testing.T) {
	a := assert.New(t)
	f := newPrimeField(t, NTTFriendlyPrime)

	pr := NewPolyRing(f)

	for i := range 8 {
		n := 1 << (i + 1)

		p := randomPolynomial(f, 4242+uint64(i), n)
		coeffs := append([]uint64(nil), p.NoCopySlice()...)

		// Forward: the reference evaluates at the powers of w, and the stage loops
		// must land on the same values in the same order.
		a.NoError(pr.NttForward(p))
		a.Equal(nttRecursive(f, coeffs), p.NoCopySlice(), "forward, n=%d", n)

		// Backward: both invert the same evaluations onto the coefficients they came
		// from.
		evals := append([]uint64(nil), p.NoCopySlice()...)

		a.NoError(pr.NttBackward(p))
		a.Equal(coeffs, nttRecursiveInverse(f, evals), "backward, n=%d", n)
		a.Equal(coeffs, p.NoCopySlice(), "backward, n=%d", n)
	}
}

func TestBitReverseInPlace(t *testing.T) {
	a := []uint64{0, 1, 2, 3, 4, 5, 6, 7}
	bitReverseInPlace(a)
	fmt.Println(a)
}
