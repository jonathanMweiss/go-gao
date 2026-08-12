// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

type NttEvaluator struct {
	cache *evaluationCache

	pr field.PolyRing
}

func NewNttEvaluator(f field.Field) *NttEvaluator {
	return &NttEvaluator{
		pr:    field.NewDensePolyRing(f),
		cache: newEvaluatorCache(),
	}
}

// supportsSize reports whether the field admits an NTT of length n, which requires n to
// be a power of two of at least 2 that divides p-1. NewCodeParameters calls this so a
// bad n surfaces as an error rather than a panic from inside Encode.
func (e *NttEvaluator) supportsSize(n int) error {
	if n <= 0 {
		return errNonPositiveN
	}

	_, err := e.pr.GetField().GetRootOfUnity(uint64(n))

	return err
}

// EvaluationPoints returns the n-th roots of unity used as evaluation points.
// The returned slice is a copy: mutating it does not disturb the internal cache.
//
// It panics if the field does not admit an NTT of length n. Construct the code through
// NewCodeParameters, which rejects such an n with ErrUnsupportedSize.
func (e *NttEvaluator) EvaluationPoints(n int) []uint64 {
	return slices.Clone(e.evaluationPoints(n))
}

func (e *NttEvaluator) evaluationPoints(n int) []uint64 {
	if points := e.cache.loadPoints(n); points != nil {
		return points
	}

	if err := e.supportsSize(n); err != nil {
		panic(fmt.Sprintf("gao: NttEvaluator cannot evaluate at %d points: %v", n, err))
	}

	// make polynomial p(x) = x.
	// then attempt to compute its NTT.
	inner := make([]uint64, n)
	inner[1] = 1
	p := e.pr.NewPolynomial(inner, false)

	if err := e.pr.NttForward(p); err != nil {
		panic(fmt.Sprintf("gao: NTT of length %d failed: %v", n, err))
	}

	points := p.NoCopySlice()

	e.cache.storePoints(n, points)

	return points
}

func (e *NttEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

func (e *NttEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if err := e.pr.NttForward(p); err != nil {
		return nil, err
	}

	return p.NoCopySlice(), nil
}

func (e *NttEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
	// The locator polynomial L(x) = (x - x_1)(x - x_2)...(x - x_n)
	// where x_1, x_2, ..., x_n are the n-th roots of unity
	// is L(x) = x^n - 1
	f := e.pr.GetField()
	inner := make([]uint64, n+1)
	inner[0] = f.Neg(1)
	inner[n] = 1
	return e.pr.NewPolynomial(inner, false)
}

// supports fast Gao: decoding can invert the transform instead of interpolating.
func (e *NttEvaluator) isNTT() bool {
	return true
}
