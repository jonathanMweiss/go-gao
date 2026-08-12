// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"slices"
	"sync"

	"github.com/jonathanmweiss/go-gao/field"
)

// EvaluationMap supplies the n points a codeword is evaluated at, and the machinery to
// evaluate a polynomial over them. Implementations can be fast transforms such as an
// NTT, or plain pointwise polynomial evaluation.
//
// The interface is deliberately sealed: it carries unexported methods, so it cannot be
// implemented outside this package. Use NewNttEvaluator or NewSlowEvaluator.
type EvaluationMap interface {
	// has access to a specific prime field.
	PrimeField() field.Field
	// returns the evaluation points for a polynomial of degree n.
	// The returned slice is owned by the caller and safe to modify.
	EvaluationPoints(n int) (xs []uint64)

	// might change the polynomial
	EvaluatePolynomial(p *field.Polynomial) (ys []uint64, err error)

	// The locator polynomial for the evaluation points.
	// Namely, given the evaluation points x_1, ..., x_n, the locator polynomial is
	// L(x) = (x - x_1)(x - x_2)...(x - x_n)
	GenerateLocatorPolynomial(n int) *field.Polynomial

	// supportsSize reports whether this map can produce n evaluation points over its
	// field, so callers can fail with an error instead of panicking later.
	supportsSize(n int) error

	// evaluationPoints returns the cached points without copying. Callers must not
	// modify the result. Code uses this on the encode/decode path, where it only
	// reads, to avoid copying n elements per operation.
	evaluationPoints(n int) []uint64

	isNTT() bool
}

var errNonPositiveN = errors.New("codeword length `n` must be positive")

type evaluationCache struct {
	sync.Locker
	degreeToPoints map[int][]uint64
}

func (e *evaluationCache) storePoints(n int, points []uint64) {
	e.Lock()
	defer e.Unlock()

	if _, ok := e.degreeToPoints[n]; ok {
		return
	}

	e.degreeToPoints[n] = points
}

// SlowEvaluator evaluates polynomials pointwise at the points 1, 2, ..., n using
// classical arithmetic. Unlike NttEvaluator it places no constraint on n beyond
// 0 < n < p, so it is the option for a codeword length that is not a power of two
// dividing p-1.
//
// It is markedly slower: encode and decode are quadratic in n rather than
// quasi-linear, so it is intended for small codes. Prefer NttEvaluator when the
// field and n permit.
type SlowEvaluator struct {
	cache *evaluationCache

	pr field.PolyRing
}

func (e *evaluationCache) loadPoints(n int) []uint64 {
	e.Lock()
	defer e.Unlock()

	if points, ok := e.degreeToPoints[n]; ok {
		return points
	}

	return nil
}

func NewSlowEvaluator(f field.Field) *SlowEvaluator {
	return &SlowEvaluator{
		pr:    field.NewDensePolyRing(f),
		cache: newEvaluatorCache(),
	}
}

func newEvaluatorCache() *evaluationCache {
	return &evaluationCache{
		Locker:         &sync.Mutex{},
		degreeToPoints: make(map[int][]uint64),
	}
}

// EvaluationPoints returns the points 1, 2, ..., n used to evaluate a codeword.
// The returned slice is a copy: mutating it does not disturb the internal cache.
func (e *SlowEvaluator) EvaluationPoints(n int) []uint64 {
	return slices.Clone(e.evaluationPoints(n))
}

func (e *SlowEvaluator) evaluationPoints(n int) []uint64 {
	if points := e.cache.loadPoints(n); points != nil {
		return points
	}

	points := make([]uint64, n)
	for i := range points {
		points[i] = uint64(i + 1)
	}

	e.cache.storePoints(n, points)

	return points
}

var errNotInCoefficientForm = errors.New("polynomial not in coefficient form")

func (e *SlowEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

func (e *SlowEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if !p.IsCoeffMode() {
		return nil, errNotInCoefficientForm
	}

	points := e.EvaluationPoints(len(p.ToSlice()))
	values := make([]uint64, len(points))

	pr := e.pr
	for i, x := range points {
		values[i] = pr.Evaluate(p, x)
	}

	return values, nil
}

func (e *SlowEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
	xs := e.EvaluationPoints(n)
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
func (e *SlowEvaluator) supportsSize(n int) error {
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
func (e *SlowEvaluator) isNTT() bool {
	return false
}
