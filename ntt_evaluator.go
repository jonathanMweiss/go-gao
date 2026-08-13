// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"

	"github.com/jonathanmweiss/go-gao/field"
)

// nttEvaluator evaluates polynomials at the n-th roots of unity using the number
// theoretic transform, which is quasi-linear rather than the quadratic pointwise
// evaluation of slowEvaluator.
//
// It requires n to be a power of two dividing p-1. NewCodeParameters checks this and
// reports ErrUnsupportedSize rather than letting the evaluator fail later.
//
// An nttEvaluator holds no mutable state and is safe for concurrent use.
type nttEvaluator struct {
	pr field.PolyRing
}

func newNttEvaluator(f field.Field) *nttEvaluator {
	return &nttEvaluator{pr: field.NewDensePolyRing(f)}
}

// supportsSize reports whether the field admits an NTT of length n, which requires n to
// be a power of two of at least 2 that divides p-1. NewCodeParameters calls this so a
// bad n surfaces as an error rather than a panic from inside Encode.
func (e *nttEvaluator) supportsSize(n int) error {
	if n <= 0 {
		return errNonPositiveN
	}

	_, err := e.pr.GetField().GetRootOfUnity(uint64(n))

	return err
}

// EvaluationPoints returns the n-th roots of unity used as evaluation points.
// Each call builds a fresh slice, so the caller may modify it freely.
//
// It panics if the field does not admit an NTT of length n. Construct the code through
// NewCodeParameters, which rejects such an n with ErrUnsupportedSize.
func (e *nttEvaluator) EvaluationPoints(n int) []uint64 {
	if err := e.supportsSize(n); err != nil {
		panic(fmt.Sprintf("gao: nttEvaluator cannot evaluate at %d points: %v", n, err))
	}

	// The roots of unity are the NTT of p(x) = x.
	inner := make([]uint64, n)
	inner[1] = 1
	p := e.pr.NewPolynomial(inner, false)

	if err := e.pr.NttForward(p); err != nil {
		panic(fmt.Sprintf("gao: NTT of length %d failed: %v", n, err))
	}

	return p.NoCopySlice()
}

func (e *nttEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

func (e *nttEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if err := e.pr.NttForward(p); err != nil {
		return nil, err
	}

	return p.NoCopySlice(), nil
}

func (e *nttEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
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
func (e *nttEvaluator) isNTT() bool {
	return true
}
