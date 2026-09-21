// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

// evaluationMap supplies the n points a codeword is evaluated at, and the machinery to
// evaluate a polynomial over them. Implementations can be fast transforms such as an
// NTT, or plain pointwise polynomial evaluation.
//
// The interface is deliberately sealed: it carries unexported methods, so it cannot be
// implemented outside this package. Use newNttEvaluator or newSlowEvaluator.
type evaluationMap interface {
	// has access to a specific prime field.
	PrimeField() field.Field
	// returns the n evaluation points, where n is the codeword length this map was
	// built for. The returned slice is owned by the caller and safe to modify.
	EvaluationPoints() (xs []uint64)

	// EvaluateCoeffs returns the polynomial with these coefficients evaluated at the n
	// evaluation points, zero-padding where coeffs is shorter. coeffs holds at most n
	// values, and is read, never written, so a caller may keep it.
	EvaluateCoeffs(coeffs []uint64) (ys []uint64, err error)

	// Interpolate is the inverse: it returns the polynomial taking these n values at
	// the n evaluation points.
	//
	// Unlike EvaluateCoeffs it takes ownership of ys, which it may modify or keep, so
	// a caller passes a slice of its own. The decoder's is scratch it has finished
	// with, and sparing it a copy of the whole codeword is worth the asymmetry.
	Interpolate(ys []uint64) (p *field.Polynomial, err error)

	// The locator polynomial for the evaluation points.
	// Namely, given the evaluation points x_1, ..., x_n, the locator polynomial is
	// L(x) = (x - x_1)(x - x_2)...(x - x_n)
	GenerateLocatorPolynomial() *field.Polynomial

	isNTT() bool
}

var errNonPositiveN = errors.New("codeword length `n` must be positive")

// slowEvaluator evaluates polynomials pointwise at the points 1, 2, ..., n using
// classical arithmetic. Unlike nttEvaluator it places no constraint on n beyond
// 0 < n < p, so it is the option for a codeword length that is not a power of two
// dividing p-1.
//
// It is markedly slower: encode and decode are quadratic in n rather than
// quasi-linear, so it is intended for small codes. Prefer nttEvaluator when the
// field and n permit.
//
// The ring is the Code's own, shared rather than duplicated; see nttEvaluator.
//
// An evaluator is built complete and never written to afterwards, and the ring it
// borrows is safe for concurrent use, so an evaluator is too.
type slowEvaluator struct {
	pr           *field.PolyRing
	interpolator *field.Interpolator
	n            int
	xs           []uint64 // the evaluation points, read-only after construction.
}

// newSlowEvaluator builds the evaluator for codeword length n, or reports why the field
// cannot serve that many points. The points are 1, 2, ..., n.
func newSlowEvaluator(pr *field.PolyRing, n int) (*slowEvaluator, error) {
	if n <= 0 {
		return nil, errNonPositiveN
	}

	if uint64(n) >= pr.GetField().Modulus() {
		return nil, errNTooLargeForField
	}

	xs := make([]uint64, n)
	for i := range xs {
		xs[i] = uint64(i + 1)
	}

	return &slowEvaluator{
		pr:           pr,
		interpolator: field.NewInterpolator(pr),
		n:            n,
		xs:           xs,
	}, nil
}

// EvaluationPoints returns the points 1, 2, ..., n used to evaluate a codeword.
// Each call returns a fresh slice, cloned from the evaluator's own.
func (e *slowEvaluator) EvaluationPoints() []uint64 {
	return slices.Clone(e.xs)
}

func (e *slowEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

// EvaluateCoeffs evaluates pointwise. NewPolynomial reduces its slice in place, so the
// coefficients are copied rather than wrapped.
func (e *slowEvaluator) EvaluateCoeffs(coeffs []uint64) ([]uint64, error) {
	p := e.pr.NewPolynomial(slices.Clone(coeffs), false)

	values := make([]uint64, e.n)

	pr := e.pr
	for i, x := range e.xs {
		values[i] = pr.Evaluate(p, x)
	}

	return values, nil
}

// Interpolate recovers the polynomial from its values by Lagrange interpolation. It
// happens to leave ys alone, but callers may not rely on that.
func (e *slowEvaluator) Interpolate(ys []uint64) (*field.Polynomial, error) {
	return e.interpolator.Interpolate(e.xs, ys)
}

func (e *slowEvaluator) GenerateLocatorPolynomial() *field.Polynomial {
	polys := make([]*field.Polynomial, e.n)

	f := e.pr.GetField()
	for i, x := range e.xs {
		// create m_i(x) = (x - x_i)
		coeffs := make([]uint64, 2)
		coeffs[1] = 1
		coeffs[0] = f.Neg(f.Reduce(x))

		polys[i] = e.pr.NewPolynomial(coeffs, false)
	}

	return e.pr.Product(polys)
}

var errNTooLargeForField = errors.New("n must be smaller than the field modulus")

// does not support fast Gao.
func (e *slowEvaluator) isNTT() bool {
	return false
}
