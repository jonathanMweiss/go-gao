// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

// An ErasureSet names the codeword positions a decode treats as unknown. It holds the
// erasure locator S(x) = product of (x - xi) over those positions, and S evaluated at
// every evaluation point.
//
// Building that is the expensive half of an erasure decode, and it depends on the
// positions alone. Words that lost the same positions share one set:
//
//	es, err := code.Erasures(3, 17, 42)
//	for _, word := range words {
//		msg, err := code.Decode(word, es)
//	}
//
// The zero value is the empty set: no erasures declared.
//
// A set is read-only once built, safe for concurrent use, and suits any code built with
// the same modulus, n, k and evaluation strategy. Decode reports ErrForeignErasureSet
// for any other.
type ErasureSet struct {
	params codeParams
	at     []int

	// s is the erasure locator, sVals is s evaluated at every point of the code.
	s     *field.Polynomial
	sVals []uint64

	stopDegree int
}

// Erasures builds the erasure set for positions at of this code's codewords.
//
// It returns ErrErasureOutOfRange or ErrDuplicateErasure for a malformed at, and
// ErrTooManyMissingPoints if more than n-k positions are named. With no arguments it
// returns the empty set.
func (gao *Code) Erasures(at ...int) (ErasureSet, error) {
	erased, err := gao.checkErasures(at)
	if err != nil {
		return ErasureSet{}, err
	}

	if len(erased) == 0 {
		return ErasureSet{}, nil
	}

	s := gao.createErasureLocator(erased, gao.xs)

	// EvaluatePolynomial transforms its argument in place and takes its length as the
	// point count, so s is padded to n and handed over as a copy: the set keeps s for
	// the divisions in erasureOnlyMessage and recoverMessage.
	inner := make([]uint64, gao.N())
	copy(inner, s.NoCopySlice())

	sVals, err := gao.eval.EvaluatePolynomial(gao.pr.NewPolynomial(inner, false))
	if err != nil {
		return ErasureSet{}, err
	}

	return ErasureSet{
		params: gao.params(),
		at:     erased,
		s:      s,
		sVals:  sVals,
		// each erasure raises the stop degree by half of what an error does.
		stopDegree: (gao.N() + gao.K() + len(erased)) / 2,
	}, nil
}

// Erasures builds the erasure set for the symbols the lost byte ranges cover. A symbol
// any range touches is erased whole, and ranges may overlap, repeat, or fall partly
// outside the codeword.
func (bc *ByteCode) Erasures(lost ...ByteRange) (ErasureSet, error) {
	return bc.code.Erasures(bc.erasedSymbols(lost)...)
}

// Len is the number of erased positions.
func (e ErasureSet) Len() int { return len(e.at) }

// empty reports whether the set declares no erasures, leaving s, sVals and stopDegree
// unset.
func (e ErasureSet) empty() bool { return e.s == nil }

// codeParams is what the contents of an ErasureSet depend on: the modulus and the
// evaluation points fix the locator and its values, n and k fix the stop degree. The two
// strategies place their points differently -- 1..n against powers of a root of unity --
// so the strategy belongs here too.
type codeParams struct {
	mod  uint64
	n, k int
	ntt  bool
}

func (gao *Code) params() codeParams {
	return codeParams{
		mod: gao.PrimeField().Modulus(),
		n:   gao.N(),
		k:   gao.K(),
		ntt: gao.eval.isNTT(),
	}
}

// validFor reports whether the set may be used with gao. The empty set suits any code.
func (e ErasureSet) validFor(gao *Code) error {
	if e.empty() || e.params == gao.params() {
		return nil
	}

	return ErrForeignErasureSet
}

// create the erasure locator polynomial S(x) = product of (x - xi) for xi an evaluation point corresponding to an erased index.
// This is similar to the locator Polynomial g0=product of (x - xi) for all evaluation points, but only for the erased indices.
// Note S(x) is distinct from the error locator E(x) of the README: E is never formed explicitly, it
// falls out of the partial GCD as the Bezout coefficient v.
func (gao *Code) createErasureLocator(erasedIndices []int, xs []uint64) *field.Polynomial {
	f := gao.pr.GetField()
	polys := make([]*field.Polynomial, len(erasedIndices))
	for i, idx := range erasedIndices {
		coeffs := make([]uint64, 2)
		coeffs[1] = 1
		coeffs[0] = f.Neg(f.Reduce(xs[idx]))
		polys[i] = gao.pr.NewPolynomial(coeffs, false)
	}

	// complexity: O(n log^2 n)
	return gao.pr.Product(polys)
}

// checkErasures validates caller-supplied erasure indices.
func (gao *Code) checkErasures(erasedAt []int) ([]int, error) {
	if len(erasedAt) == 0 {
		return nil, nil
	}

	seen := make(map[int]struct{}, len(erasedAt))

	for _, idx := range erasedAt {
		if idx < 0 || idx >= gao.N() {
			return nil, fmt.Errorf("%w: %d not in [0, %d)", ErrErasureOutOfRange, idx, gao.N())
		}

		if _, dup := seen[idx]; dup {
			return nil, fmt.Errorf("%w: %d", ErrDuplicateErasure, idx)
		}

		seen[idx] = struct{}{}
	}

	// dervied from 2e+s <= n-k where e=0 and s=len(erasedAt).
	if len(erasedAt) > gao.N()-gao.K() {
		return nil, ErrTooManyMissingPoints
	}

	return slices.Clone(erasedAt), nil
}
