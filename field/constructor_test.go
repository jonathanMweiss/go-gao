// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestEmptyCoefficientsYieldZeroPolynomial: an empty sum of terms is zero, so an empty
// coefficient slice is the zero polynomial rather than an error.
func TestEmptyCoefficientsYieldZeroPolynomial(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	for name, p := range map[string]*Polynomial{
		"nil":           pr.NewPolynomial(nil, false),
		"empty":         pr.NewPolynomial([]uint64{}, false),
		"explicit zero": pr.NewPolynomial([]uint64{0}, false),
	} {
		t.Run(name, func(t *testing.T) {
			assert.True(t, p.IsZero(), "should be the zero polynomial")
			assert.Equal(t, -1, p.Degree(), "zero polynomial has degree -1")
			assert.Equal(t, uint64(0), p.LeadCoeff())
		})
	}
}

func TestIsZero(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	for _, tc := range []struct {
		name string
		in   []uint64
		want bool
	}{
		{"empty", nil, true},
		{"single zero", []uint64{0}, true},
		{"all zero", []uint64{0, 0, 0}, true},
		{"constant", []uint64{5}, false},
		{"x", []uint64{0, 1}, false},
		{"x^3", []uint64{0, 0, 0, 1}, false},
		{"5x^2, zero below", []uint64{0, 0, 5}, false},
		{"5 + x", []uint64{5, 1}, false},
		{"trailing zeros above a term", []uint64{0, 7, 0, 0}, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			assert.Equal(t, tc.want, pr.NewPolynomial(tc.in, false).IsZero())
		})
	}
}

// TestZeroPolynomialIsUsableInRingOps: constructing the zero polynomial is only useful
// if the ring then accepts it wherever a zero divisor is not implied.
func TestZeroPolynomialIsUsableInRingOps(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	zero := pr.NewPolynomial(nil, false)
	p := pr.NewPolynomial([]uint64{1, 2, 3}, false)

	sum := &Polynomial{}
	pr.Add(p, zero, sum)
	assert.True(t, sum.Equals(p), "p + 0 == p")

	prod := &Polynomial{}
	pr.Mul(p, zero, prod)
	assert.True(t, prod.IsZero(), "p * 0 == 0")

	// Dividing BY zero remains a broken invariant.
	assert.Panics(t, func() { pr.Div(p, zero) }, "division by the zero polynomial must panic")

	// Dividing zero by p is well defined.
	q, rem := pr.Div(zero, p)
	assert.True(t, q.IsZero())
	assert.True(t, rem.IsZero())
}

// TestNilFieldPanics: a ring is meaningless without a field. Catching it here is what
// lets PolyRing.NewPolynomial construct without any field check at all — the ring is
// now the only route to a Polynomial, so a valid ring implies a valid field.
func TestNilFieldPanics(t *testing.T) {
	assert.PanicsWithValue(t, "NewDensePolyRing: nil field", func() {
		NewDensePolyRing(nil)
	})
}

// TestProductEmptyIsOne: an empty product is the multiplicative identity, mirroring an
// empty coefficient slice yielding the zero polynomial.
func TestProductEmptyIsOne(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	one := pr.Product(nil)
	assert.Equal(t, 0, one.Degree())
	assert.Equal(t, uint64(1), one.LeadCoeff())

	// Multiplying by the empty product leaves a polynomial unchanged.
	p := pr.NewPolynomial([]uint64{3, 1, 4}, false)
	got := &Polynomial{}
	pr.Mul(p, one, got)
	assert.True(t, got.Equals(p), "p * (empty product) == p")
}

// TestProductDoesNotAliasOrMutateInput: productTree hands leaf polynomials to Mul
// uncopied, which is only safe because Mul reads its operands. Pin both halves of that:
// the inputs survive untouched, and the result is never the caller's polynomial.
func TestProductDoesNotAliasOrMutateInput(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	// Enough factors to build a real tree, and to cross nttMulThreshold on the way up.
	const factors = 40

	polys := make([]*Polynomial, factors)
	want := make([][]uint64, factors)

	for i := range polys {
		// (x - i)
		polys[i] = pr.NewPolynomial([]uint64{f.Neg(uint64(i + 1)), 1}, false)
		want[i] = polys[i].ToSlice()
	}

	got := pr.Product(polys)
	assert.Equal(t, factors, got.Degree(), "product of %d monic linear factors", factors)

	for i, p := range polys {
		assert.Equal(t, want[i], p.ToSlice(), "Product must not modify input %d", i)
		assert.NotSame(t, p, got, "result must not alias an input")
	}

	// The single-element case must copy, not hand back the caller's polynomial.
	single := pr.NewPolynomial([]uint64{7, 1}, false)
	out := pr.Product([]*Polynomial{single})
	assert.NotSame(t, single, out, "single-element Product must return a copy")
	assert.True(t, out.Equals(single))
}

// TestRingNewPolynomialUsesRingField: the whole point of the method form is that the
// field cannot be mismatched or forgotten.
func TestRingNewPolynomialUsesRingField(t *testing.T) {
	f, err := NewPrimeField(65537)
	require.NoError(t, err)

	pr := NewDensePolyRing(f)

	p := pr.NewPolynomial([]uint64{1, 2, 3}, false)
	assert.Equal(t, f.Modulus(), p.f.Modulus())
}
