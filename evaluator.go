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
	// returns the evaluation points for a polynomial of degree n.
	// The returned slice is owned by the caller and safe to modify.
	EvaluationPoints(n int) (xs []uint64)

	// EvaluatePolynomial returns p evaluated at the n evaluation points, zero-padding
	// p where it is shorter. It leaves p untouched, so a caller may keep it.
	EvaluatePolynomial(p *field.Polynomial, n int) (ys []uint64, err error)

	// The locator polynomial for the evaluation points.
	// Namely, given the evaluation points x_1, ..., x_n, the locator polynomial is
	// L(x) = (x - x_1)(x - x_2)...(x - x_n)
	GenerateLocatorPolynomial(n int) *field.Polynomial

	// supportsSize reports whether this map can produce n evaluation points over its
	// field, so callers can fail with an error instead of panicking later.
	supportsSize(n int) error

	isNTT() bool
}

var errNonPositiveN = errors.New("codeword length `n` must be positive")

// pointCache memoizes evaluation points per codeword length. An evaluator is shared by
// a Code, which is safe for concurrent use, and deriving the points costs a whole NTT
// for nttEvaluator.
type pointCache struct {
	mu sync.RWMutex
	m  map[int][]uint64
}

// get returns the cached points for n, deriving them with build on a miss. build runs
// outside the lock, so two racing misses may both derive; the points are deterministic,
// and the first stored wins for everyone.
func (c *pointCache) get(n int, build func() []uint64) []uint64 {
	c.mu.RLock()
	xs, ok := c.m[n]
	c.mu.RUnlock()

	if ok {
		return xs
	}

	xs = build()

	c.mu.Lock()
	defer c.mu.Unlock()

	if stored, ok := c.m[n]; ok {
		return stored
	}

	if c.m == nil {
		c.m = make(map[int][]uint64, 1)
	}

	c.m[n] = xs

	return xs
}

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
// A slowEvaluator adds no mutable state of its own, and the ring it borrows is safe for
// concurrent use, so an evaluator is too.
type slowEvaluator struct {
	pr     *field.PolyRing
	points pointCache
}

func newSlowEvaluator(pr *field.PolyRing) *slowEvaluator {
	return &slowEvaluator{pr: pr}
}

// EvaluationPoints returns the points 1, 2, ..., n used to evaluate a codeword.
// Each call returns a fresh slice, cloned from the cached one.
func (e *slowEvaluator) EvaluationPoints(n int) []uint64 {
	return slices.Clone(e.cachedPoints(n))
}

// cachedPoints returns the shared points for n, which callers must not modify.
func (e *slowEvaluator) cachedPoints(n int) []uint64 {
	return e.points.get(n, func() []uint64 {
		xs := make([]uint64, n)
		for i := range xs {
			xs[i] = uint64(i + 1)
		}

		return xs
	})
}

var errNotInCoefficientForm = errors.New("polynomial not in coefficient form")

func (e *slowEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

// EvaluatePolynomial evaluates p pointwise, which reads p without modifying it. A p
// shorter than n needs no padding here: the missing coefficients are zero either way.
func (e *slowEvaluator) EvaluatePolynomial(p *field.Polynomial, n int) ([]uint64, error) {
	if !p.IsCoeffMode() {
		return nil, errNotInCoefficientForm
	}

	values := make([]uint64, n)

	pr := e.pr
	for i, x := range e.cachedPoints(n) {
		values[i] = pr.Evaluate(p, x)
	}

	return values, nil
}

func (e *slowEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
	xs := e.cachedPoints(n)
	polys := make([]*field.Polynomial, n)

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
func (e *slowEvaluator) supportsSize(n int) error {
	if n <= 0 {
		return errNonPositiveN
	}

	if uint64(n) >= e.pr.GetField().Modulus() {
		return errNTooLargeForField
	}

	return nil
}

var errNTooLargeForField = errors.New("n must be smaller than the field modulus")

// does not support fast Gao.
func (e *slowEvaluator) isNTT() bool {
	return false
}
