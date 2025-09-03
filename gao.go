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
		g0:           gao.g0.Copy(),
		interpolator: field.NewInterpolator(gao.pr),
		stopDegree:   gao.stopDegree,
	}
}

var ErrDataTooLarge = errors.New("data too large")
var ErrDataElementsTooLarge = errors.New("data elements too large")

func (gao *Code) Encode(data []uint64) (map[uint64]uint64, error) {
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
	xs, ys, err := gao.prepareDecoding(received)
	if err != nil {
		return nil, err
	}

	var f, r *field.Polynomial
	if gao.EvaluationMap.isNTT() {
		f, r, err = gao.decodeNTT(ys, xs)
	} else {
		f, r, err = gao.decodeGeneric(ys, xs)
	}

	if err != nil {
		return nil, err
	}

	if !r.IsZero() || f.Degree() > gao.K() {
		return nil, ErrDecoding
	}

	return f.ToSlice(), nil
}

/*
prepare the decoding process by filling in missing evaluated points with zeros.
*/
func (gao *Code) prepareDecoding(toDecode map[uint64]uint64) ([]uint64, []uint64, error) {
	if len(toDecode) > gao.N() {
		return nil, nil, ErrTooManyPoints
	}

	numMissing := 0

	xs := gao.EvaluationMap.EvaluationPoints(gao.N())
	for _, x := range xs {
		if _, ok := toDecode[x]; !ok {
			toDecode[x] = 0
			numMissing += 1
		}
	}

	if numMissing > gao.MaxErrors() {
		return nil, nil, ErrTooManyMissingPoints
	}

	ys := make([]uint64, gao.N())
	for i, x := range xs {
		ys[i] = toDecode[x] // according to the order of the EvaluationMap's EvaluationPoints.
	}

	return xs, ys, nil
}

func (gao *Code) decodeGeneric(ys []uint64, xs []uint64) (*field.Polynomial, *field.Polynomial, error) {
	g1, err := gao.interpolator.Interpolate(xs, ys)
	if err != nil {
		return nil, nil, err
	}

	pr := gao.pr

	g, _, v := pr.PartialExtendedEuclidean(gao.g0, g1, gao.stopDegree)
	f, r := pr.LongDiv(g, v)

	return f, r, nil
}

func (gao *Code) decodeNTT(ys []uint64, xs []uint64) (*field.Polynomial, *field.Polynomial, error) {
	g1 := field.NewPolynomial(gao.pr.GetField(), ys, true)
	if err := gao.pr.NttBackward(g1); err != nil {
		return nil, nil, err
	}

	pr := gao.pr

	g, _, v := pr.NttPartialExtendedEuclidean(gao.g0, g1, gao.stopDegree)
	f, r := pr.LongDivNTT(g, v)

	return f, r, nil
}
