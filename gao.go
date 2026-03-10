package gao

import (
	"errors"

	"example.com/pir2peer/gao/field"
)

type Coder interface {
	EvaluationMap

	// redundancy value
	N() int

	// data size
	K() int

	// maximum number of errors that can be corrected.
	MaxErrors() int
}

type Decoder interface {
	Coder
	Decode(encodedData map[uint64]uint64) ([]uint64, error)
}

type Encoder interface {
	Coder
	Encode(data []uint64) (map[uint64]uint64, error)
}

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

func (c *CodeParams) N() int {
	return c.n
}

func (c *CodeParams) K() int {
	return c.k
}

func (c *CodeParams) MaxErrors() int {
	return c.maxErrors
}

var ErrNSmallerThanK = errors.New("redundancy value `n` must be greater than or equal to data size `k`")

func NewCodeParameters(e EvaluationMap, n, k int) (CodeParams, error) {
	if n < k {
		return CodeParams{}, ErrNSmallerThanK
	}

	return CodeParams{
		EvaluationMap: e,
		n:             n,
		k:             k,
		maxErrors:     (n - k) / 2,
	}, nil
}

func NewCodeGao(c CodeParams) *Code {
	fld := c.EvaluationMap.PrimeField()
	pr := field.NewDensePolyRing(fld)
	// create g0(x) = (x - x_1)(x - x_2)...(x - x_n)
	// TODO: for FastEvaluationMaps, we can skip this step, and create g0 without computing it.

	return &Code{
		CodeParams:   c,
		pr:           pr,
		g0:           c.EvaluationMap.GenerateLocatorPolynomial(c.N()),
		interpolator: field.NewInterpolator(pr),
		stopDegree:   (c.N() + c.K()) / 2,
	}
}

func (gao *Code) Copy() *Code {
	return &Code{
		CodeParams:   gao.CodeParams,
		pr:           gao.pr,
		interpolator: gao.interpolator,
		stopDegree:   gao.stopDegree,
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
	xs := gao.EvaluationMap.EvaluationPoints(gao.N())
	points := make(map[uint64]uint64, gao.N())

	for i, y := range ys {
		points[xs[i]] = y
	}

	return points, nil
}

var ErrTooManyMissingPoints = errors.New("too many missing points")
var ErrTooManyPoints = errors.New("too many evaluated points")
var ErrDecoding = errors.New("decoding error")

func (gao *Code) Decode(received map[uint64]uint64) ([]uint64, error) {
	// fill missing evaluated points with 0.
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

	numMissing := 0
	erasedIndices := make([]int, 0)

	xs := gao.EvaluationMap.EvaluationPoints(gao.N())
	for i, x := range xs {
		if _, ok := toDecode[x]; !ok {
			toDecode[x] = 0
			numMissing += 1
			erasedIndices = append(erasedIndices, i)
		}
	}

	if numMissing > gao.N()-gao.K() {
		return nil, nil, nil, ErrTooManyMissingPoints
	}

	ys := make([]uint64, gao.N())
	for i, x := range xs {
		ys[i] = toDecode[x] // according to the order of the EvaluationMap's EvaluationPoints.
	}

	return xs, ys, erasedIndices, nil
}

// full intuitive explanation in README.md
func (gao *Code) decodeGeneric(ys []uint64, xs []uint64, erased []int) (*field.Polynomial, *field.Polynomial, error) {
	var E *field.Polynomial
	stopDegree := gao.stopDegree
	fld := gao.pr.GetField()

	if len(erased) > 0 {
		E = gao.createErasureLocator(erased, xs)
		// scale ys by E(xi)
		for i, x := range xs {
			valE := gao.pr.Evaluate(E, x)
			ys[i] = fld.Mul(ys[i], valE)
		}
		stopDegree = (gao.N() + gao.K() + len(erased)) / 2
	}

	g1, err := gao.interpolator.Interpolate(xs, ys)
	if err != nil {
		return nil, nil, err
	}

	pr := gao.pr

	g, _, v := pr.FastPartialGCD(gao.g0, g1, stopDegree)

	if len(erased) > 0 {
		// G = g/v
		G, remG := pr.Div(g, v)
		if !remG.IsZero() {
			return nil, nil, ErrDecoding
		}
		// f = G/E
		f, remF := pr.Div(G, E)
		return f, remF, nil
	}

	f, r := pr.Div(g, v)

	return f, r, nil
}

func (gao *Code) decodeNTT(ys, xs []uint64, erased []int) (*field.Polynomial, *field.Polynomial, error) {
	var E *field.Polynomial
	stopDegree := gao.stopDegree
	fld := gao.pr.GetField()

	if len(erased) > 0 {
		E = gao.createErasureLocator(erased, xs)

		// Evaluate E(x) at all points xs
		eInner := make([]uint64, gao.N())
		copy(eInner, E.NoCopySlice())
		Epoly := field.NewPolynomial(fld, eInner, false)
		if err := gao.pr.NttForward(Epoly); err != nil {
			return nil, nil, err
		}
		eVals := Epoly.NoCopySlice()

		// scale ys by E(xi)
		for i := range ys {
			ys[i] = fld.Mul(ys[i], eVals[i])
		}
		stopDegree = (gao.N() + gao.K() + len(erased)) / 2
	}

	g1 := field.NewPolynomial(gao.pr.GetField(), ys, true)
	if err := gao.pr.NttBackward(g1); err != nil {
		return nil, nil, err
	}

	pr := gao.pr

	g, _, v := pr.FastPartialGCD(gao.g0, g1, stopDegree)

	if len(erased) > 0 {
		// G = g/v
		G, remG := pr.Div(g, v)
		if !remG.IsZero() {
			return nil, nil, ErrDecoding
		}
		// f = G/E
		f, remF := pr.Div(G, E)
		return f, remF, nil
	}

	f, r := pr.Div(g, v)

	return f, r, nil
}

// create the erasure locator polynomial E(x) = product of (x - xi) for xi an evaluation point corresponding to an erased index.
// This is similar to the locator Polynomial g0=product of (x - xi) for all evaluation points, but only for the erased indices.
func (gao *Code) createErasureLocator(erasedIndices []int, xs []uint64) *field.Polynomial {
	f := gao.pr.GetField()
	polys := make([]*field.Polynomial, len(erasedIndices))
	for i, idx := range erasedIndices {
		coeffs := make([]uint64, 2)
		coeffs[1] = 1
		coeffs[0] = f.Neg(f.Reduce(xs[idx]))
		polys[i] = field.NewPolynomial(f, coeffs, false)
	}

	// complexity: O(n log^2 n)
	return field.PolyProduct(gao.pr, polys)
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
	p := field.NewPolynomial(f, paddedData, false)
	// evaluate polynomial at n points.

	ys, err := gao.EvaluationMap.EvaluatePolynomial(p)
	if err != nil {
		return nil, err
	}

	return ys, nil
}

// Notice: This might change the input slice (depending on the EvaluationMap used).
// TODO: Change API to receive XS, YS so that we can infer erasures from missing points, and so that we can avoid modifying the input slice.
func (gao *Code) DecodeFromSlice(ys []uint64) ([]uint64, error) {
	return gao.sliceDecode(gao.EvaluationMap.EvaluationPoints(gao.N()), ys, nil)
}
