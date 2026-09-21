// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
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

	sVals, err := gao.evaluateEverywhere(s)
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

// evaluateEverywhere returns p evaluated at each of the code's evaluation points: one
// forward transform on the NTT path, point by point otherwise.
func (gao *Code) evaluateEverywhere(p *field.Polynomial) ([]uint64, error) {
	if !gao.eval.isNTT() {
		out := make([]uint64, gao.N())
		for i, x := range gao.xs {
			out[i] = gao.pr.Evaluate(p, x)
		}

		return out, nil
	}

	inner := make([]uint64, gao.N())
	copy(inner, p.NoCopySlice())

	poly := gao.pr.NewPolynomial(inner, false)
	if err := gao.pr.NttForward(poly); err != nil {
		return nil, err
	}

	return poly.NoCopySlice(), nil
}
