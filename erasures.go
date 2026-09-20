// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"github.com/jonathanmweiss/go-gao/field"
)

// An ErasureSet names the codeword positions a decode should treat as unknown, and holds
// the work those positions imply.
//
// That work is the expensive half of an erasure decode and none of it depends on the
// received word: building the locator S(x) = product of (x - xi) over the erased points,
// evaluating it at every evaluation point, and inverting its reversal as a power series
// so that dividing a word by it costs one multiplication. Words that lost the same
// positions can therefore share one set:
//
//	es, err := code.Erasures(3, 17, 42)
//	for _, word := range words {
//		msg, err := code.Decode(word, es)
//	}
//
// The zero value is the empty set, which decodes a word with no erasures declared.
//
// A set suits any code built with the same parameters, not only the one that built it,
// so codes constructed separately on either side of a link can share one. Decode reports
// ErrForeignErasureSet for a set whose parameters do not match. A set is read-only once
// built, so it is safe for concurrent use.
type ErasureSet struct {
	params codeParams
	at     []int

	// sVals is the erasure locator evaluated at every point of the code, and div is that
	// same locator prepared as a divisor, since every word in a batch is divided by it.
	sVals []uint64
	div   *field.DivisorCache
	sDeg  int

	stopDegree int
}

// Erasures prepares at as the erased positions of a codeword of this code.
//
// It returns ErrErasureOutOfRange or ErrDuplicateErasure for a malformed at, and
// ErrTooManyMissingPoints if more than n-k positions are named. Calling it with no
// arguments returns the empty set.
func (gao *Code) Erasures(at ...int) (ErasureSet, error) {
	erased, err := gao.checkErasures(at)
	if err != nil {
		return ErasureSet{}, err
	}

	if len(erased) == 0 {
		return ErasureSet{}, nil
	}

	locator := gao.createErasureLocator(erased, gao.xs)

	sVals, err := gao.evaluateEverywhere(locator)
	if err != nil {
		return ErasureSet{}, err
	}

	// A decode divides by the locator once per word, and every such quotient is the
	// message, so k coefficients is the longest one a successful decode produces.
	return ErasureSet{
		params: gao.params(),
		at:     erased,
		sVals:  sVals,
		div:    gao.pr.NewDivisorCache(locator, gao.K()),
		sDeg:   locator.Degree(),
		// each erasure raises the stop degree by half of what an error does.
		stopDegree: (gao.N() + gao.K() + len(erased)) / 2,
	}, nil
}

// Erasures prepares the byte ranges lost as the erased positions of a codeword of this
// code. A symbol any range touches is erased whole, and ranges may overlap, repeat, or
// fall partly outside the codeword.
func (bc *ByteCode) Erasures(lost ...ByteRange) (ErasureSet, error) {
	return bc.code.Erasures(bc.erasedSymbols(lost)...)
}

// Len is the number of erased positions.
func (e ErasureSet) Len() int { return len(e.at) }

// empty reports whether the set declares no erasures, in which case none of the derived
// fields are set and the decode takes its error-only path.
func (e ErasureSet) empty() bool { return e.div == nil }

// codeParams is everything the contents of an ErasureSet depend on. The locator and its
// evaluations follow from the modulus and the evaluation points, and the stop degree
// from n and k as well. The two strategies place their points differently -- 1..n
// against powers of a root of unity -- so which one is in use is part of the identity.
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

// validFor reports whether the set may be used with gao. The empty set declares nothing
// and suits any code.
func (e ErasureSet) validFor(gao *Code) error {
	if e.empty() || e.params == gao.params() {
		return nil
	}

	return ErrForeignErasureSet
}

// evaluateEverywhere returns p evaluated at each of the code's evaluation points.
//
// The NTT path gets all n values from one forward transform; the pointwise path has to
// evaluate point by point.
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
