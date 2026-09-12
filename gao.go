// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

// Code is a Reed-Solomon code that decodes with Gao's algorithm.
//
// A *Code is immutable after construction and safe for concurrent use.
type Code struct {
	// eval is a named field, not embedded: embedding would promote the evaluator's
	// methods onto Code and make decoding internals part of the public API.
	eval evaluationMap

	n         int
	k         int
	maxErrors int

	pr           field.PolyRing
	interpolator *field.Interpolator
	// g0 polynomial from the Gao code.
	// with fast EvaluationMaps like NTT, this polynomial can be used to do fast division.
	g0 *field.Polynomial

	// xs are the n evaluation points. n is fixed for the life of a Code, so these are
	// computed once here rather than re-derived — an NTT per call, for nttEvaluator.
	// Read-only after construction, which is what makes a *Code safe to share.
	xs []uint64

	stopDegree int
}

// N returns the codeword length: the number of evaluation points the encoder emits.
func (gao *Code) N() int {
	return gao.n
}

// K returns the message length: the number of data symbols carried by a codeword.
func (gao *Code) K() int {
	return gao.k
}

// MaxErrors returns the number of corrupted symbols the code can repair when their
// positions are unknown, (n-k)/2. Erasures are cheaper: each one costs half an error,
// so a decode succeeds while 2*errors+erasures <= n-k.
func (gao *Code) MaxErrors() int {
	return gao.maxErrors
}

// UsesNTT reports whether this code evaluates via the number theoretic transform
// (quasi-linear) rather than pointwise (quadratic in n).
//
// Worth checking after NewCode without RequireNTT, since the fallback is silent.
func (gao *Code) UsesNTT() bool {
	return gao.eval.isNTT()
}

// EvaluationPoints returns the n points a codeword is evaluated at, in the order Encode
// emits them and Decode expects them.
//
// Callers do not need this to encode or decode -- both APIs are positional. It is here
// for interoperating with another Reed-Solomon implementation, and for understanding
// what a code is doing.
//
// The returned slice is a copy and may be modified freely.
func (gao *Code) EvaluationPoints() []uint64 {
	return slices.Clone(gao.xs)
}

// PrimeField returns the field this code operates over.
func (gao *Code) PrimeField() field.Field {
	return gao.pr.GetField()
}

// Errors reported by NewCode, Encode and Decode. Test for them with errors.Is.
var (
	// Construction.
	ErrNSmallerThanK   = errors.New("redundancy value `n` must be greater than or equal to data size `k`")
	ErrNonPositiveK    = errors.New("data size `k` must be positive")
	ErrUnsupportedSize = errors.New("evaluation strategy does not support the requested codeword length `n`")

	// Encoding: the message does not fit the code, or an element is outside the field.
	ErrDataTooLarge         = errors.New("data too large")
	ErrDataElementsTooLarge = errors.New("data elements too large")

	// Decoding.
	ErrTooManyMissingPoints = errors.New("too many missing points")
	ErrMismatchedLengths    = errors.New("codeword length does not match the code's n")
	ErrErasureOutOfRange    = errors.New("erasure index out of range")
	ErrDuplicateErasure     = errors.New("duplicate erasure index")
	// ErrDecoding means no message is consistent with the points given, so the error
	// and erasure budget was exceeded.
	ErrDecoding = errors.New("decoding error")
)

// An Option adjusts how NewCode selects its evaluation strategy.
type Option func(*config)

type config struct {
	requireNTT bool
	forceSlow  bool
}

// RequireNTT makes NewCode fail rather than fall back to pointwise evaluation, which is
// quadratic in n.
//
// The NTT strategy needs two transforms from the field, not one: an n-point transform to
// evaluate a codeword, and a 2n-point transform for the products the decoder's partial
// GCD takes. So size the prime against 2n -- p=65537, whose p-1 is 2^16, evaluates at up
// to 65536 points but is only usable up to n=32768.
//
// A field short of either falls back silently by default, and at large n pointwise
// evaluation is the difference between milliseconds and minutes. Use this option
// whenever the quasi-linear path is a requirement rather than a preference, or check
// UsesNTT afterwards.
func RequireNTT() Option {
	return func(c *config) { c.requireNTT = true }
}

// Pointwise forces classical pointwise evaluation even where an NTT is available.
// Intended for benchmarking and differential testing against the fast path.
func Pointwise() Option {
	return func(c *config) { c.forceSlow = true }
}

// NewCode builds a Reed-Solomon code carrying k data symbols in an n-symbol codeword
// over the prime field f, decoding with Gao's algorithm.
//
// By default it evaluates via the number theoretic transform when f and n allow; which
// requires n to be a power of two dividing p-1. Otherwise it falls back to pointwise
// evaluation, which accepts any 0 < n < p but has O(n^2) runtime. Pass RequireNTT to turn
// that fallback into an error, or Pointwise to force the classical path. Check UsesNTT
// to see which was chosen.
//
// It returns ErrNonPositiveK if k <= 0, ErrNSmallerThanK if n < k, and ErrUnsupportedSize
// if no permitted strategy can evaluate at n points.
func NewCode(f field.Field, n, k int, opts ...Option) (*Code, error) {
	var cfg config

	for _, opt := range opts {
		// A nil Option is a no-op, so callers can pass one conditionally.
		if opt == nil {
			continue
		}

		opt(&cfg)
	}

	if k <= 0 {
		return nil, ErrNonPositiveK
	}

	if n < k {
		return nil, ErrNSmallerThanK
	}

	eval, err := selectEvaluator(f, n, cfg)
	if err != nil {
		return nil, err
	}

	pr := field.NewDensePolyRing(f)

	return &Code{
		eval:      eval,
		n:         n,
		k:         k,
		maxErrors: (n - k) / 2,
		pr:        pr,
		// g0(x) = (x - x_1)(x - x_2)...(x - x_n)
		g0:           eval.GenerateLocatorPolynomial(n),
		xs:           eval.EvaluationPoints(n),
		interpolator: field.NewInterpolator(pr),
		stopDegree:   (n + k) / 2,
	}, nil
}

// selectEvaluator resolves the strategy, preferring the NTT unless told otherwise.
func selectEvaluator(f field.Field, n int, cfg config) (evaluationMap, error) {
	slow := newSlowEvaluator(f)

	if cfg.forceSlow {
		if err := slow.supportsSize(n); err != nil {
			return nil, fmt.Errorf("%w: n=%d: %w", ErrUnsupportedSize, n, err)
		}

		return slow, nil
	}

	ntt := newNttEvaluator(f)

	nttErr := ntt.supportsSize(n)
	if nttErr == nil {
		return ntt, nil
	}

	if cfg.requireNTT {
		return nil, fmt.Errorf(
			"%w: n=%d: RequireNTT was set but this field admits no NTT usable at that length "+
				"(n and 2n must both be powers of two dividing p-1, p=%d): %w",
			ErrUnsupportedSize, n, f.Modulus(), nttErr)
	}

	if err := slow.supportsSize(n); err != nil {
		return nil, fmt.Errorf("%w: n=%d: %w", ErrUnsupportedSize, n, err)
	}

	return slow, nil
}

// Decode recovers the original message from a received codeword, repairing both
// corrupted values and missing ones.
//
// ys holds one value per evaluation point, in the order EvaluationPoints returns them.
// erasedAt lists the indices of positions known to be unusable; whatever ys holds at an
// erased index is ignored, so there is no need to blank those entries first.
//
// An erasure is cheaper than an error precisely because its position is known: decoding
// succeeds while 2*errors+erasures <= n-k.
//
// Beyond that budget correction is not guaranteed. Decode returns ErrDecoding when it
// detects an inconsistency, but with enough errors a received word can be pushed closer
// to a different valid codeword, and then it returns a confidently wrong message. That
// is inherent to Reed-Solomon codes, not to this implementation.
//
// ys is not modified, and the returned message always has length k, zero-padded if the
// message it recovers has high-order zero symbols.
//
// It returns ErrMismatchedLengths if ys is not n long, ErrErasureOutOfRange or
// ErrDuplicateErasure for a malformed erasedAt, ErrTooManyMissingPoints if more than n-k
// positions are erased, and ErrDecoding if no message is consistent with what it was
// given.
func (gao *Code) Decode(ys []uint64, erasedAt ...int) ([]uint64, error) {
	if len(ys) != gao.N() {
		return nil, ErrMismatchedLengths
	}

	erased, err := gao.checkErasures(erasedAt)
	if err != nil {
		return nil, err
	}

	// The decode reduces and transforms its values in place, so it works on a copy.
	work := slices.Clone(ys)
	gao.reduceSlice(work)

	// The all-zero message is degenerate for Gao's algorithm and has to be settled here,
	// before the partial GCD ever sees it. See zeroCodewordIsNearest.
	if gao.zeroCodewordIsNearest(work, erased) {
		return make([]uint64, gao.K()), nil
	}

	var f, r *field.Polynomial
	if gao.eval.isNTT() {
		f, r, err = gao.decodeNTT(work, gao.xs, erased)
	} else {
		f, r, err = gao.decodeGeneric(work, gao.xs, erased)
	}

	if err != nil {
		return nil, err
	}

	if !r.IsZero() || f.Degree() >= gao.K() {
		return nil, ErrDecoding
	}

	return gao.messageOf(f), nil
}

// messageOf returns f's coefficients at exactly length k.
//
// f.ToSlice() stops at f's degree, so a message whose high-order symbols are zero --
// [10, 20, 30, 0] -- would otherwise come back shorter than it went in, and a caller
// comparing what it encoded against what it decoded would see a spurious mismatch.
// Callers have already checked deg(f) < k, so nothing is truncated here.
func (gao *Code) messageOf(f *field.Polynomial) []uint64 {
	msg := make([]uint64, gao.K())
	copy(msg, f.ToSlice())

	return msg
}

// zeroCodewordIsNearest reports whether the all-zero codeword is the one nearest to ys,
// which is to say whether ys decodes to the all-zero message.
//
// the zero message is f(x) = 0, zero evaluated anywhere is zero, so it
// encodes to n zeros and the distance to it is free to measure.
//
// It is worth the O(n) scan because f = 0 cannot go through the partial GCD at all. The
// Berlekamp-Welch product Q = E*f is then identically zero, carrying no degree for the
// GCD to stop at.
func (gao *Code) zeroCodewordIsNearest(ys []uint64, erased []int) bool {
	isErased := make([]bool, len(ys))
	for _, idx := range erased {
		isErased[idx] = true
	}

	// Distance from ys to the all-zero codeword: a non-erased position holding anything
	// other than zero is one symbol of disagreement.
	distance := 0

	for i, y := range ys {
		if !isErased[i] && y != 0 {
			distance++
		}
	}

	// The ordinary Reed-Solomon condition, with that distance standing in for the errors.
	return 2*distance+len(erased) <= gao.N()-gao.K()
}

func (gao *Code) reduceSlice(ys []uint64) {
	fld := gao.pr.GetField()
	for i := range ys {
		ys[i] = fld.Reduce(ys[i])
	}
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

	return gao.recoverMessage(g1, S, stopDegree)
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

	return gao.recoverMessage(g1, S, stopDegree)
}

// recoverMessage runs the partial GCD and strips the locators from what it returns.
func (gao *Code) recoverMessage(g1, S *field.Polynomial, stopDegree int) (f, rem *field.Polynomial, err error) {
	pr := gao.pr

	g, _, v := pr.FastPartialGCD(gao.g0, g1, stopDegree)

	// A zero error locator would make the divisions below panic. zeroCodewordIsNearest has
	// already answered the one input that produces one -- an all-zero received word --
	// so reaching this means the received word admits no consistent message.
	if v.IsZero() {
		return nil, nil, ErrDecoding
	}

	if S == nil {
		f, rem = pr.Div(g, v)

		return f, rem, nil
	}

	G, remG := pr.Div(g, v)
	if !remG.IsZero() {
		return nil, nil, ErrDecoding
	}

	f, rem = pr.Div(G, S)

	return f, rem, nil
}

// codewordMessage returns the optimistic error-free decode result for g1 — the
// inverse-NTT / interpolant of the received points — when deg(g1) < K. `received` is
// then exactly a codeword and g1 is its message. The returned pair mirrors the normal
//
// The deg(g1) < 0 branch is defensive: an all-zero g1 means an all-zero received word,
// which zeroCodewordIsNearest has already answered.
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

// Encode encodes up to k data symbols into an n-symbol codeword.
//
// The returned values are positional: index i is the evaluation at EvaluationPoints()[i],
// and Decode expects them in that order. Fewer than k symbols are zero-padded.
//
// It returns ErrDataTooLarge if data holds more than k symbols, and
// ErrDataElementsTooLarge if any symbol is not less than the field modulus.
func (gao *Code) Encode(data []uint64) ([]uint64, error) {
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

	ys, err := gao.eval.EvaluatePolynomial(p)
	if err != nil {
		return nil, err
	}

	return ys, nil
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
