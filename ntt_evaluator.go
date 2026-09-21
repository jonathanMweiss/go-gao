// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"slices"
	"sync"

	"github.com/jonathanmweiss/go-gao/field"
)

// nttEvaluator evaluates polynomials at the n-th roots of unity using the number
// theoretic transform, which is quasi-linear rather than the quadratic pointwise
// evaluation of slowEvaluator.
//
// It requires n to be a power of two dividing p-1. NewCode checks this and
// reports ErrUnsupportedSize rather than letting the evaluator fail later.
type nttEvaluator struct {
	pr *field.PolyRing // safe for concurrent use.
	n  int

	// xs is derived on first use and never rewritten; once guards that one write.
	once sync.Once
	xs   []uint64
}

func newNttEvaluator(pr *field.PolyRing, n int) *nttEvaluator {
	return &nttEvaluator{pr: pr, n: n}
}

// supportsSize reports whether the field admits the transforms this strategy needs: an
// n-point one to evaluate with, and a 2n-point one to decode with. Both require a power
// of two of at least 2 dividing p-1, so in practice the second is the binding one.
//
// The 2n is not a margin. The decoder's partial GCD multiplies polynomials of degree up
// to n, and the longest convolution that asks for measures 1.25n, which rounds up to a
// 2n-point transform. A field offering only the n-point one evaluates quickly and then
// multiplies in schoolbook inside a recursion built to avoid it, which measures slower
// than never taking that recursion at all -- so it does not count as support.
//
// NewCode calls this, so an n the strategy cannot serve surfaces as a fallback or an
// error rather than as a panic from inside Encode.
func (e *nttEvaluator) supportsSize() error {
	if e.n <= 0 {
		return errNonPositiveN
	}

	fld := e.pr.GetField()

	if _, err := field.RootOfUnity(fld, uint64(e.n)); err != nil {
		return err
	}

	if _, err := field.RootOfUnity(fld, uint64(2*e.n)); err != nil {
		return fmt.Errorf("decoding needs a 2n-point transform, and 2n=%d does not divide p-1=%d: %w",
			2*e.n, fld.Modulus()-1, err)
	}

	return nil
}

// EvaluationPoints returns the n-th roots of unity used as evaluation points.
// Each call returns a fresh slice, cloned from the cached one, so the caller may modify
// it freely.
//
// It panics if the field does not admit an NTT of length n. Construct the code through
// NewCode, which rejects such an n with ErrUnsupportedSize.
func (e *nttEvaluator) EvaluationPoints() []uint64 {
	return slices.Clone(e.points())
}

// points returns the shared roots of unity, which callers must not modify. Deriving
// them is a full transform, so they are derived once.
//
// The size check sits outside the once, so an evaluator built with an n this field
// cannot serve panics on every call rather than yielding nil after the first.
func (e *nttEvaluator) points() []uint64 {
	if err := e.supportsSize(); err != nil {
		panic(fmt.Sprintf("gao: nttEvaluator cannot evaluate at %d points: %v", e.n, err))
	}

	e.once.Do(func() {
		// The roots of unity are the NTT of p(x) = x.
		inner := make([]uint64, e.n)
		inner[1] = 1
		p := e.pr.NewPolynomial(inner, false)

		if err := e.pr.NttForward(p); err != nil {
			panic(fmt.Sprintf("gao: NTT of length %d failed: %v", e.n, err))
		}

		e.xs = p.NoCopySlice()
	})

	return e.xs
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
