// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"errors"
	"fmt"
	"slices"

	"github.com/jonathanmweiss/go-gao/field"
)

// A Codeword is the n evaluations [Code.Encode] produces and [Code.Decode] consumes.
// each uint64 is a symbol in the code.
type Codeword []uint64

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

	pr           *field.PolyRing
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
	// ErrForeignErasureSet means the set was built for a code with other parameters.
	ErrForeignErasureSet = errors.New("erasure set belongs to a different code")
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

	pr := field.NewPolyRing(f)

	eval, err := selectEvaluator(pr, n, cfg)
	if err != nil {
		return nil, err
	}

	return &Code{
		eval:      eval,
		n:         n,
		k:         k,
		maxErrors: (n - k) / 2,
		pr:        pr,
		// g0(x) = (x - x_1)(x - x_2)...(x - x_n)
		g0:           eval.GenerateLocatorPolynomial(),
		xs:           eval.EvaluationPoints(),
		interpolator: field.NewInterpolator(pr),
		stopDegree:   (n + k) / 2,
	}, nil
}

// selectEvaluator resolves the strategy, preferring the NTT unless told otherwise.
func selectEvaluator(pr *field.PolyRing, n int, cfg config) (evaluationMap, error) {
	if !cfg.forceSlow {
		ntt, nttErr := newNttEvaluator(pr, n)
		if nttErr == nil {
			return ntt, nil
		}

		if cfg.requireNTT {
			return nil, fmt.Errorf(
				"%w: n=%d: RequireNTT was set but this field admits no NTT usable at that length "+
					"(n and 2n must both be powers of two dividing p-1, p=%d): %w",
				ErrUnsupportedSize, n, pr.GetField().Modulus(), nttErr)
		}
	}

	slow, err := newSlowEvaluator(pr, n)
	if err != nil {
		return nil, fmt.Errorf("%w: n=%d: %w", ErrUnsupportedSize, n, err)
	}

	return slow, nil
}

// Decode recovers the original message from a received codeword, repairing both
// corrupted values and missing ones.
//
// ys holds one value per evaluation point, in the order EvaluationPoints returns them.
// erasures names the positions known to be unusable, and [Code.Erasures] builds it.
// Whatever ys holds at an erased index is ignored, so there is no need to blank those
// entries first. Pass the zero ErasureSet when nothing is missing.
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
// It returns ErrMismatchedLengths if ys is not n long, ErrForeignErasureSet if erasures
// was built for other parameters, and ErrDecoding if no message is consistent with what
// it was given.
func (gao *Code) Decode(ys Codeword, erasures ErasureSet) ([]uint64, error) {
	if len(ys) != gao.N() {
		return nil, ErrMismatchedLengths
	}

	if err := erasures.validFor(gao); err != nil {
		return nil, err
	}

	// The decode reduces and transforms its values in place, so it works on a copy.
	work := slices.Clone(ys)
	gao.reduceSlice(work)

	// The all-zero message is degenerate for Gao's algorithm and has to be settled here,
	// before the partial GCD ever sees it. See zeroCodewordIsNearest.
	if gao.zeroCodewordIsNearest(work, erasures.at) {
		return make([]uint64, gao.K()), nil
	}

	var (
		f, r *field.Polynomial
		err  error
	)

	if gao.eval.isNTT() {
		f, r, err = gao.decodeNTT(work, erasures)
	} else {
		f, r, err = gao.decodeGeneric(work, erasures)
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
	msg := f.ToSlice()
	if len(msg) < gao.K() {
		msg = append(msg, make([]uint64, gao.K()-len(msg))...)
	}

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
func (gao *Code) decodeGeneric(ys []uint64, erasures ErasureSet) (*field.Polynomial, *field.Polynomial, error) {
	stopDegree := gao.stopDegree

	if !erasures.empty() {
		// scale ys by S(xi). interpolating the scaled values yields g1*S mod g0 (see README.md for reason).
		fld := gao.pr.GetField()
		for i := range ys {
			ys[i] = fld.Mul(ys[i], erasures.sVals[i])
		}

		stopDegree = erasures.stopDegree
	}

	g1, err := gao.interpolator.Interpolate(gao.xs, ys)
	if err != nil {
		return nil, nil, err
	}

	// Optimistic error-free path:
	// When g_1 has degree < K, it'll be the first polynomial in the Euclidean
	// remainder sequence below stopDegree, so PartialGCD would
	// thus, GCD returns g=g1, v=1 with r=0.
	// Since the return value `f` is defined f=g1/v (and in this case v=1), we return g1 directly.
	// This is true only when there are no erasures.
	if erasures.empty() && g1.Degree() < gao.K() {
		f, r := gao.codewordMessage(g1)
		return f, r, nil
	}

	if f, r, ok := gao.erasureOnlyMessage(g1, erasures); ok {
		return f, r, nil
	}

	return gao.recoverMessage(g1, erasures, stopDegree)
}

func (gao *Code) decodeNTT(ys []uint64, erasures ErasureSet) (*field.Polynomial, *field.Polynomial, error) {
	stopDegree := gao.stopDegree

	if !erasures.empty() {
		// (see README.md for reason)
		// scale ys by S(xi): the inverse NTT below then yields (g1*S mod g0).
		fld := gao.pr.GetField()
		for i := range ys {
			ys[i] = fld.Mul(ys[i], erasures.sVals[i])
		}

		stopDegree = erasures.stopDegree
	}

	g1 := gao.pr.NewPolynomial(ys, true)
	if err := gao.pr.NttBackward(g1); err != nil {
		return nil, nil, err
	}

	// Optimistic error-free path:
	// When g_1 has degree < K, it'll be the first polynomial in the Euclidean
	// remainder sequence below stopDegree, so PartialGCD would
	// thus, GCD returns g=g1, v=1 with r=0.
	// Since the return value `f` is defined f=g1/v (and in this case v=1), we return g1 directly.
	// This is true only when there are no erasures.
	if erasures.empty() && g1.Degree() < gao.K() {
		f, r := gao.codewordMessage(g1)
		return f, r, nil
	}

	if f, r, ok := gao.erasureOnlyMessage(g1, erasures); ok {
		return f, r, nil
	}

	return gao.recoverMessage(g1, erasures, stopDegree)
}

// erasureOnlyMessage recovers the message directly from g1 and S, skipping the partial
// GCD, when the word carries erasures but no errors.
//
// Suppose by contradiction that one of the points has an error and deg(g1) < K+s.
// the following proof will reach a contradiction (meaning that if deg(g1) < k+s, there are no errors).
// denote g1(x_i)= (f(x_i) + e_i)*S(x_i); so we have e_i*S(x_i) = g1(x_i)-f(x_i)*S(x_i).
// The right-hand side is one polynomial evaluated at x_i, so name it D = g1 - f*S:
//
//	D(x_i) = g1(x_i) - f(x_i)*S(x_i)
//	       = (f(x_i) + e_i)*S(x_i) - f(x_i)*S(x_i)
//	       = e_i*S(x_i)
//
// D vanishes wherever the word is clean, which is nearly everywhere:
// e_i=0 when there are no errors (n-t such locations) and S(x_i)=0 when there are erasures.
// So with t symbols in error, D has at least n-t distinct roots,
// and is therefore either zero or of degree >= n-t.
//
// Now suppose the test below passes, deg(g1) < K+s, writing s=deg(S) for the number of
// erasures. That inequality is strict and degrees are integers, so deg(g1) <= K+s-1. A
// message polynomial carries K coefficients, deg(f) <= K-1, so deg(f*S) <= K-1+s too.
// D is their difference, hence deg(D) <= K+s-1.
//
// from the decoder's budget,
// 2t+s <= n-K, we have t <= (n-K-s)/2 and hence
//
//	n-t >= n-(n-K-s)/2 = (2n-n+K+s)/2 = (n+K+s)/2 >= K+s
//
// the last step because n >= K+s, which the budget gives at t=0 and checkErasures
// enforces outright: (n+K+s)/2 >= (K+s+K+s)/2 = K+s.
// So a nonzero D would have deg(D) >= n-t >= K+s.
//
// The assumed error sits at a non-erased point: there e_j != 0 and S(x_j) != 0, so
// D(x_j) != 0 and D is nonzero; leaving deg(D) >= K+s and deg(D) <= K+s-1. Absurd,
// so no such error exists, and D = 0 gives g1 = f*S: the division below returns f.
//
// The contrapositive is the guard itself: a word carrying t >= 1 errors has
// deg(g1) >= K+s, fails the test, and goes on to the partial GCD. Nothing error-free is
// turned away either, since deg(f*S) <= K-1+s < n means the mod g0 never bites.
func (gao *Code) erasureOnlyMessage(g1 *field.Polynomial, erasures ErasureSet) (f, r *field.Polynomial, ok bool) {
	if erasures.empty() {
		return nil, nil, false
	}

	if g1.Degree() >= gao.K()+erasures.s.Degree() {
		return nil, nil, false
	}

	f, r = gao.pr.Div(g1, erasures.s)
	// Scaling by S zeroed the erased positions, so g1 vanishes there and S always divides it.
	// should never happen, but check anyway.
	if !r.IsZero() {
		return nil, nil, false
	}

	return f, r, true
}

// recoverMessage runs the partial GCD and strips the locators from what it returns.
func (gao *Code) recoverMessage(g1 *field.Polynomial, erasures ErasureSet, stopDegree int) (f, rem *field.Polynomial, err error) {
	pr := gao.pr

	g, _, v := pr.PartialGCD(gao.g0, g1, stopDegree)

	// can't divide by zero.
	if v.IsZero() {
		return nil, nil, ErrDecoding
	}

	G, remG := pr.Div(g, v)

	if erasures.empty() {
		return G, remG, nil
	}

	if !remG.IsZero() {
		return nil, nil, ErrDecoding
	}

	f, rem = pr.Div(G, erasures.s)

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

// Encode encodes up to k data symbols into an n-symbol codeword.
//
// The returned values are positional: index i is the evaluation at EvaluationPoints()[i],
// and Decode expects them in that order. Fewer than k symbols are zero-padded.
//
// It returns ErrDataTooLarge if data holds more than k symbols, and
// ErrDataElementsTooLarge if any symbol is not less than the field modulus.
func (gao *Code) Encode(data []uint64) (Codeword, error) {
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

	ys, err := gao.eval.EvaluateCoeffs(data)
	if err != nil {
		return nil, err
	}

	return ys, nil
}
