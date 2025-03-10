package gao

import (
	"errors"

	"github.com/jonathanmweiss/go-gao/field"
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

	xs := c.EvaluationMap.EvaluationPoints(c.N())

	// create g0(x) = (x - x_1)(x - x_2)...(x - x_n)
	// TODO: for FastEvaluationMaps, we can skip this step, and create g0 without computing it.

	polys := make([]*field.Polynomial, c.N())

	for i, x := range xs {
		// create m_i(x) = (x - x_i)
		coeffs := make([]field.Elem, 2)
		coeffs[1] = fld.ElemFromUint64(1)
		coeffs[0] = fld.ElemFromUint64(x).Neg()

		polys[i] = field.NewPolynomial(coeffs, false)
	}

	return &Code{
		CodeParams:   c,
		g0:           fld.PolyProduct(polys),
		interpolator: field.NewInterpolator(fld),
		stopDegree:   (c.N() + c.K()) / 2,
	}
}

var ErrDataTooLarge = errors.New("data too large")
var ErrDataElementsTooLarge = errors.New("data elements too large")

func (gao *Code) Encode(data []uint64) (map[uint64]uint64, error) {
	f := gao.PrimeField()

	q := f.Prime()
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
	coefficients := f.ElemSlice(paddedData)
	p := field.NewPolynomial(coefficients, false)
	// evaluate polynomial at n points.

	ys, err := gao.EvaluationMap.EvaluatePolynomial(p)
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
	// fill missing points with 0.
	if len(received) > gao.N() {
		return nil, ErrTooManyPoints
	}

	numMissing := 0

	xs := gao.EvaluationMap.EvaluationPoints(gao.N())
	for _, x := range xs {
		if _, ok := received[x]; !ok {
			received[x] = 0
			numMissing += 1
		}
	}

	if numMissing > gao.MaxErrors() {
		return nil, ErrTooManyMissingPoints
	}

	ys := make([]uint64, gao.N())
	for i, x := range xs {
		ys[i] = received[x] // according to the order of the EvaluationMap's EvaluationPoints.
	}

	g1, err := gao.interpolator.Interpolate(xs, ys)
	if err != nil {
		return nil, err
	}

	g, _, v := field.PartialExtendedEuclidean(gao.g0, g1, gao.stopDegree)
	f, r := g.LongDiv(v)

	if !r.IsZero() || f.Degree() > gao.K() {
		return nil, ErrDecoding
	}

	return f.ToSlice(), nil
}
