// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"slices"
	"sync"

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

	// The locator polynomial for the evaluation points.
	// Namely, given the evaluation points x_1, ..., x_n, the locator polynomial is
	// L(x) = (x - x_1)(x - x_2)...(x - x_n)
	GenerateLocatorPolynomial() *field.Polynomial

	// supportsSize reports whether this map can produce its n evaluation points over
	// its field, so callers can fail with an error instead of panicking later.
	supportsSize() error

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
// Its one piece of mutable state is the memo below, written once under a sync.Once, and
// the ring it borrows is safe for concurrent use, so an evaluator is too.
type slowEvaluator struct {
	pr *field.PolyRing
	n  int

	// xs is derived on first use and never rewritten; once guards that one write.
	once sync.Once
	xs   []uint64
}

func newSlowEvaluator(pr *field.PolyRing, n int) *slowEvaluator {
	return &slowEvaluator{pr: pr, n: n}
}

// EvaluationPoints returns the points 1, 2, ..., n used to evaluate a codeword.
// Each call returns a fresh slice, cloned from the shared one.
func (e *slowEvaluator) EvaluationPoints() []uint64 {
	return slices.Clone(e.points())
}

// points returns the shared evaluation points, which callers must not modify. They are
// derived once: an evaluator serves one codeword length for its whole life.
func (e *slowEvaluator) points() []uint64 {
	e.once.Do(func() {
		e.xs = make([]uint64, e.n)
		for i := range e.xs {
			e.xs[i] = uint64(i + 1)
		}
	})

	return e.xs
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
	for i, x := range e.points() {
		values[i] = pr.Evaluate(p, x)
	}

	return values, nil
}

func (e *slowEvaluator) GenerateLocatorPolynomial() *field.Polynomial {
	xs := e.points()
	polys := make([]*field.Polynomial, e.n)

	f := e.pr.GetField()
	for i, x := range xs {
		// create m_i(x) = (x - x_i)
		coeffs := make([]uint64, 2)
		coeffs[1] = 1
		coeffs[0] = f.Neg(f.Reduce(x))

		polys[i] = e.pr.NewPolynomial(coeffs, false)
	}

	return e.pr.Product(polys)
}

// supportsSize accepts any positive n: the points are simply 1..n, which requires
// nothing of the field beyond having at least n distinct non-zero elements.
func (e *slowEvaluator) supportsSize() error {
	if e.n <= 0 {
		return errNonPositiveN
	}

	if uint64(e.n) >= e.pr.GetField().Modulus() {
		return errNTooLargeForField
	}

	return nil
}

var errNTooLargeForField = errors.New("n must be smaller than the field modulus")

// does not support fast Gao.
func (e *slowEvaluator) isNTT() bool {
	return false
}
