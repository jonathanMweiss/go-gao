// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

type CodeParams struct {
	EvaluationMap
	n         int
	k         int
	maxErrors int
}

type Code struct {
	CodeParams
	pr           field.PolyRing
	interpolator *field.Interpolator
	// g0 polynomial from the Gao code.
	// with fast EvaluationMaps like NTT, this polynomial can be used to do fast division.
	g0 *field.Polynomial

	stopDegree int
}

// N returns the codeword length: the number of evaluation points the encoder emits.
func (c *CodeParams) N() int {
	return c.n
}

// K returns the message length: the number of data symbols carried by a codeword.
func (c *CodeParams) K() int {
	return c.k
}

// MaxErrors returns the number of corrupted symbols the code can repair when their
// positions are unknown, (n-k)/2. Erasures are cheaper: each one costs half an error,
// so a decode succeeds while 2*errors+erasures <= n-k.
func (c *CodeParams) MaxErrors() int {
	return c.maxErrors
}

var (
	ErrNSmallerThanK   = errors.New("redundancy value `n` must be greater than or equal to data size `k`")
	ErrNonPositiveK    = errors.New("data size `k` must be positive")
	ErrUnsupportedSize = errors.New("evaluation map does not support the requested codeword length `n`")
)

// NewCodeParameters validates n and k against the evaluation map and returns the
// parameters for a code that carries k data symbols in an n-symbol codeword.
//
// It returns ErrNonPositiveK if k <= 0, ErrNSmallerThanK if n < k, and
// ErrUnsupportedSize if e cannot produce n evaluation points — which for
// NttEvaluator means n must be a power of two dividing p-1.
func NewCodeParameters(e EvaluationMap, n, k int) (CodeParams, error) {
	if k <= 0 {
		return CodeParams{}, ErrNonPositiveK
	}

	if n < k {
		return CodeParams{}, ErrNSmallerThanK
	}

	if err := e.supportsSize(n); err != nil { // ensuring support for size n
		return CodeParams{}, fmt.Errorf("%w: n=%d: %w", ErrUnsupportedSize, n, err)
	}

	return CodeParams{
		EvaluationMap: e,
		n:             n,
		k:             k,
		maxErrors:     (n - k) / 2,
	}, nil
}

// NewCodeGao builds a Reed-Solomon code that decodes with Gao's algorithm from
// parameters already validated by NewCodeParameters.
//
// The returned *Code is safe for concurrent use by multiple goroutines.
func NewCodeGao(c CodeParams) *Code {
	fld := c.EvaluationMap.PrimeField()
	pr := field.NewDensePolyRing(fld)
	// create g0(x) = (x - x_1)(x - x_2)...(x - x_n)

	return &Code{
		CodeParams:   c,
		pr:           pr,
		g0:           c.EvaluationMap.GenerateLocatorPolynomial(c.N()),
		interpolator: field.NewInterpolator(pr),
		stopDegree:   (c.N() + c.K()) / 2,
	}
}

var ErrDataTooLarge = errors.New("data too large")
var ErrDataElementsTooLarge = errors.New("data elements too large")

func (gao *Code) Encode(data []uint64) (map[uint64]uint64, error) {
	ys, err := gao.EncodeToSlice(data)
	if err != nil {
		return nil, err
	}

	// create map of points.
	xs := gao.EvaluationMap.evaluationPoints(gao.N())
	points := make(map[uint64]uint64, gao.N())

	for i, y := range ys {
		points[xs[i]] = y
	}

	return points, nil
}

var ErrTooManyMissingPoints = errors.New("too many missing points")
var ErrTooManyPoints = errors.New("too many evaluated points")
var ErrDecoding = errors.New("decoding error")

// Decode recovers the original message from received evaluation points, repairing both
// corrupted values and missing ones.
//
// Points omitted from the input map `received` are treated as erasures, whose positions are therefore
// known.
// Since this is a Reed-Solomon code, the amount of errors and
// erasures must satisfy $2*errors + erasures <= n-k$; beyond that, correction
// is not guaranteed: Decode
//
// returns ErrDecoding if it detects an inconsistency, but with enough errors a codeword
// can be pushed closer to a different valid codeword, in which case it returns a
// confidently wrong message.
//
// It returns ErrTooManyPoints if received holds more than n entries,
// ErrTooManyMissingPoints if more than n-k are absent, and ErrDecoding if no message is
// consistent with the points given.
func (gao *Code) Decode(received map[uint64]uint64) ([]uint64, error) {
	xs, ys, erased, err := gao.prepareDecoding(received)
	if err != nil {
		return nil, err
	}

	return gao.sliceDecode(xs, ys, erased)
}

var ErrMismatchedLengths = errors.New("mismatched lengths of xs and ys")

func (gao *Code) sliceDecode(xs []uint64, ys []uint64, erased []int) ([]uint64, error) {
	if len(xs) != len(ys) {
		return nil, ErrMismatchedLengths
	}

	gao.reduceSlice(ys)

	var err error
	var f, r *field.Polynomial
	if gao.EvaluationMap.isNTT() {
		f, r, err = gao.decodeNTT(ys, xs, erased)
	} else {
		f, r, err = gao.decodeGeneric(ys, xs, erased)
	}

	if err != nil {
		return nil, err
	}

	if !r.IsZero() || f.Degree() >= gao.K() {
		return nil, ErrDecoding
	}

	return f.ToSlice(), nil
}

func (gao *Code) reduceSlice(ys []uint64) {
	fld := gao.pr.GetField()
	for i := range ys {
		ys[i] = fld.Reduce(ys[i])
	}
}

/*
prepare the decoding process by filling in missing evaluated points with zeros.
*/
func (gao *Code) prepareDecoding(toDecode map[uint64]uint64) ([]uint64, []uint64, []int, error) {
	if len(toDecode) > gao.N() {
		return nil, nil, nil, ErrTooManyPoints
	}

	erasedIndices := make([]int, 0)

	xs := gao.EvaluationMap.evaluationPoints(gao.N())
	ys := make([]uint64, gao.N())

	// ys follows the order of the EvaluationMap's evaluation points. A point absent from
	// toDecode is an erasure: it stays zero here, which the decoder is free to do because
	// the erasure locator annihilates whatever value sits at an erased position.
	for i, x := range xs {
		y, ok := toDecode[x]
		if !ok {
			erasedIndices = append(erasedIndices, i)

			continue
		}

		ys[i] = y
	}

	if len(erasedIndices) > gao.N()-gao.K() {
		return nil, nil, nil, ErrTooManyMissingPoints
	}

	return xs, ys, erasedIndices, nil
}

// full intuitive explanation in README.md
func (gao *Code) decodeGeneric(ys []uint64, xs []uint64, erased []int) (*field.Polynomial, *field.Polynomial, error) {
	var S *field.Polynomial

	stopDegree := gao.stopDegree
	fld := gao.pr.GetField()

	if len(erased) > 0 {
		S = gao.createErasureLocator(erased, xs)
		// scale ys by S(xi). interpolating the scaled values yields g1*S mod g0 (see README.md for reason).
		for i, x := range xs {
			valS := gao.pr.Evaluate(S, x)
			ys[i] = fld.Mul(ys[i], valS)
		}
		// each erasure raises the stop degree by half of what an error does.
		stopDegree = (gao.N() + gao.K() + len(erased)) / 2
	}

	g1, err := gao.interpolator.Interpolate(xs, ys)
	if err != nil {
		return nil, nil, err
	}

	// Optimistic error-free path:
	// When g_1 has degree < K, it'll be the first polynomial in the Euclidean
	// remainder sequence below stopDegree, so FastPartialGCD would
	// thus, GCD returns g=g1, v=1 with r=0.
	// Since the return value `f` is defined f=g1/v (and in this case v=1), we return g1 directly.
	// This is true only when there are no erasures.
	if len(erased) == 0 && g1.Degree() < gao.K() {
		f, r := gao.codewordMessage(g1)
		return f, r, nil
	}

	pr := gao.pr

	g, _, v := pr.FastPartialGCD(gao.g0, g1, stopDegree)

	if len(erased) > 0 {
		// after FastPartialGCD, g = S*v*f. We remove v, then S.
		G, remG := pr.Div(g, v)
		if !remG.IsZero() {
			return nil, nil, ErrDecoding
		}
		// remove S:
		f, remF := pr.Div(G, S)
		return f, remF, nil
	}

	f, r := pr.Div(g, v)

	return f, r, nil
}

func (gao *Code) decodeNTT(ys, xs []uint64, erased []int) (*field.Polynomial, *field.Polynomial, error) {
	var S *field.Polynomial
	stopDegree := gao.stopDegree
	fld := gao.pr.GetField()

	if len(erased) > 0 {
		S = gao.createErasureLocator(erased, xs)

		// Evaluate S(x) at all points xs
		sInner := make([]uint64, gao.N())
		copy(sInner, S.NoCopySlice())
		Spoly := gao.pr.NewPolynomial(sInner, false)
		if err := gao.pr.NttForward(Spoly); err != nil {
			return nil, nil, err
		}
		sVals := Spoly.NoCopySlice()

		// (see README.md for reason)
		// scale ys by S(xi): the inverse NTT below then yields (g1*S mod g0).
		for i := range ys {
			ys[i] = fld.Mul(ys[i], sVals[i])
		}
		stopDegree = (gao.N() + gao.K() + len(erased)) / 2
	}

	g1 := gao.pr.NewPolynomial(ys, true)
	if err := gao.pr.NttBackward(g1); err != nil {
		return nil, nil, err
	}

	// Optimistic error-free path:
	// When g_1 has degree < K, it'll be the first polynomial in the Euclidean
	// remainder sequence below stopDegree, so FastPartialGCD would
	// thus, GCD returns g=g1, v=1 with r=0.
	// Since the return value `f` is defined f=g1/v (and in this case v=1), we return g1 directly.
	// This is true only when there are no erasures.
	if len(erased) == 0 && g1.Degree() < gao.K() {
		f, r := gao.codewordMessage(g1)
		return f, r, nil
	}

	pr := gao.pr

	g, _, v := pr.FastPartialGCD(gao.g0, g1, stopDegree)

	if len(erased) > 0 {
		// after FastPartialGCD, g = S*v*f. We remove v, then S.
		G, remG := pr.Div(g, v)
		if !remG.IsZero() {
			return nil, nil, ErrDecoding
		}
		// remove S:
		f, remF := pr.Div(G, S)
		return f, remF, nil
	}

	f, r := pr.Div(g, v)

	return f, r, nil
}

// codewordMessage returns the optimistic error-free decode result for g1 — the
// inverse-NTT / interpolant of the received points — when deg(g1) < K. `received` is
// then exactly a codeword and g1 is its message. The returned pair mirrors the normal
// FastPartialGCD path: the message trimmed to its degree and a zero remainder, so sliceDecode's guard and the
// callers see identical output. deg(g1) < 0 is the all-zero message.
func (gao *Code) codewordMessage(g1 *field.Polynomial) (f, r *field.Polynomial) {
	coeffs := []uint64{0}
	if deg := g1.Degree(); deg >= 0 {
		coeffs = g1.ToSlice()[:deg+1]
	}

	f = gao.pr.NewPolynomial(coeffs, false)
	r = gao.pr.NewPolynomial(nil, false)

	return f, r
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

func (gao *Code) EncodeToSlice(data []uint64) ([]uint64, error) {
	f := gao.PrimeField()

	q := f.Modulus()
	for _, d := range data {
		if d >= q {
			return nil, ErrDataElementsTooLarge
		}
	}

	// check data length.
	if len(data) > gao.K() {
		return nil, ErrDataTooLarge
	}

	// pad:
	paddedData := make([]uint64, gao.N())
	copy(paddedData, data)

	// create polynomial from data.
	p := gao.pr.NewPolynomial(paddedData, false)
	// evaluate polynomial at n points.

	ys, err := gao.EvaluationMap.EvaluatePolynomial(p)
	if err != nil {
		return nil, err
	}

	return ys, nil
}

// DecodeFromSlice is the positional form of Decode: ys holds one value per evaluation
// point, in the order EvaluationPoints returns them, so it cannot express erasures —
// every position carries a value. Use Decode when some points are missing.
//
// It returns ErrMismatchedLengths unless len(ys) == n. DecodeFromSlice never modifies ys.
func (gao *Code) DecodeFromSlice(ys []uint64) ([]uint64, error) {
	if len(ys) != gao.N() {
		return nil, ErrMismatchedLengths
	}

	// sliceDecode reduces and transforms ys in place, so hand it a copy.
	work := slices.Clone(ys)

	return gao.sliceDecode(gao.EvaluationMap.evaluationPoints(gao.N()), work, nil)
}
