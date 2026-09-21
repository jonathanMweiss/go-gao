// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

// nttEvaluator evaluates polynomials at the n-th roots of unity using the number
// theoretic transform, which is quasi-linear rather than the quadratic pointwise
// evaluation of slowEvaluator.
//
// It needs n to be a power of two with 2n dividing p-1, checked at construction so
// that NewCode can fall back or report ErrUnsupportedSize.
type nttEvaluator struct {
	pr *field.PolyRing // safe for concurrent use.
	n  int
	xs []uint64 // the n-th roots of unity, read-only after construction.
}

// newNttEvaluator builds the evaluator for codeword length n, or reports why the field
// admits no NTT usable at that length. The points are the n-th roots of unity, derived
// here because each derivation costs a transform.
func newNttEvaluator(pr *field.PolyRing, n int) (*nttEvaluator, error) {
	if err := nttSupportsSize(pr.GetField(), n); err != nil {
		return nil, err
	}

	// The roots of unity are the NTT of p(x) = x.
	inner := make([]uint64, n)
	inner[1] = 1

	p := pr.NewPolynomial(inner, false)
	if err := pr.NttForward(p); err != nil {
		return nil, err
	}

	return &nttEvaluator{pr: pr, n: n, xs: p.NoCopySlice()}, nil
}

var errNTooSmallForNTT = errors.New("the NTT needs a codeword length of at least 2")

// nttSupportsSize reports whether the field admits both NTT transforms this strategy needs.
func nttSupportsSize(fld field.Field, n int) error {
	if n <= 0 {
		return errNonPositiveN
	}

	// The roots are derived by planting one at coefficient 1, which n=1 has no room
	// for. A one-symbol codeword carries no redundancy in any case.
	if n == 1 {
		return errNTooSmallForNTT
	}

	if _, err := field.RootOfUnity(fld, uint64(n)); err != nil {
		return err
	}

	if _, err := field.RootOfUnity(fld, uint64(2*n)); err != nil {
		return fmt.Errorf("decoding needs a 2n-point transform, and 2n=%d does not divide p-1=%d: %w",
			2*n, fld.Modulus()-1, err)
	}

	return nil
}

// EvaluationPoints returns the n-th roots of unity used as evaluation points.
// Each call returns a fresh slice, cloned from the evaluator's own, so the caller may
// modify it freely.
func (e *nttEvaluator) EvaluationPoints() []uint64 {
	return slices.Clone(e.xs)
}

func (e *nttEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

// EvaluateCoeffs transforms a buffer of its own: NttForward works in place and its
// length sets the transform's, so the one buffer both pads coeffs to n and keeps the
// caller's slice intact.
func (e *nttEvaluator) EvaluateCoeffs(coeffs []uint64) ([]uint64, error) {
	inner := make([]uint64, e.n)
	copy(inner, coeffs)

	work := e.pr.NewPolynomial(inner, false)
	if err := e.pr.NttForward(work); err != nil {
		return nil, err
	}

	return work.NoCopySlice(), nil
}

// Interpolate recovers the polynomial from its values with the inverse transform, which
// works in place: ys becomes the polynomial's own coefficients.
func (e *nttEvaluator) Interpolate(ys []uint64) (*field.Polynomial, error) {
	p := e.pr.NewPolynomial(ys, true)
	if err := e.pr.NttBackward(p); err != nil {
		return nil, err
	}

	return p, nil
}

func (e *nttEvaluator) GenerateLocatorPolynomial() *field.Polynomial {
	// The locator polynomial L(x) = (x - x_1)(x - x_2)...(x - x_n)
	// where x_1, x_2, ..., x_n are the n-th roots of unity
	// is L(x) = x^n - 1
	f := e.pr.GetField()
	inner := make([]uint64, e.n+1)
	inner[0] = f.Neg(1)
	inner[e.n] = 1
	return e.pr.NewPolynomial(inner, false)
}

// supports fast Gao: decoding can invert the transform instead of interpolating.
func (e *nttEvaluator) isNTT() bool {
	return true
}
