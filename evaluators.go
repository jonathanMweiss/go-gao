package gao

import (
	"errors"
	"sync"

	"github.com/jonathanmweiss/go-gao/field"
)

// Can be fast evaluation, like NTT, negacyclic-NTT, or just plain polynomial evaluation.
type EvaluationMap interface {
	// has access to a specific prime field.
	PrimeField() field.Field
	// returns the evaluation points for a polynomial of degree n
	EvaluationPoints(n int) (xs []uint64)

	// might change the polynomial
	EvaluatePolynomial(p *field.Polynomial) (ys []uint64, err error)

	// The locator polynomial for the evaluation points.
	// Namely, given the evaluation points x_1, ..., x_n, the locator polynomial is
	// L(x) = (x - x_1)(x - x_2)...(x - x_n)
	GenerateLocatorPolynomial(n int) *field.Polynomial

	isNTT() bool
}

type evaluationCache struct {
	sync.Locker
	degreeToPoints map[int][]uint64
}

func (e evaluationCache) storePoints(n int, points []uint64) {
	e.Lock()
	defer e.Unlock()

	if _, ok := e.degreeToPoints[n]; ok {
		return
	}

	e.degreeToPoints[n] = points
}

// Not Tested. but can support non powers of 2 for $n$ param.
type SlowEvaluator struct {
	cache *evaluationCache

	pr field.PolyRing
}

func (e *evaluationCache) loadPoints(n int) []uint64 {
	e.Lock()
	defer e.Unlock()

	if points, ok := e.degreeToPoints[n]; ok {
		return points
	}

	return nil
}

func NewSlowEvaluator(f field.Field) *SlowEvaluator {
	return &SlowEvaluator{
		pr:    field.NewDensePolyRing(f),
		cache: newEvaluatorCache(),
	}
}

func newEvaluatorCache() *evaluationCache {
	return &evaluationCache{
		Locker:         &sync.Mutex{},
		degreeToPoints: make(map[int][]uint64),
	}
}

func (e *SlowEvaluator) EvaluationPoints(n int) []uint64 {
	points := e.cache.loadPoints(n)
	if points != nil {
		return points
	}

	points = make([]uint64, n)
	for i := range points {
		points[i] = uint64(i + 1)
	}

	e.cache.storePoints(n, points)

	return points
}

var errNotInCoefficientForm = errors.New("polynomial not in coefficient form")

func (e *SlowEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

func (e *SlowEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if p.IsCoeffMode() {
		return nil, errNotInCoefficientForm
	}

	points := e.EvaluationPoints(len(p.ToSlice()))
	values := make([]uint64, len(points))

	pr := e.pr
	for i, x := range points {
		values[i] = pr.Evaluate(p, x)
	}

	return values, nil
}

func (e *SlowEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
	xs := e.EvaluationPoints(n)
	polys := make([]*field.Polynomial, n)

	f := e.pr.GetField()
	for i, x := range xs {
		// create m_i(x) = (x - x_i)
		coeffs := make([]uint64, 2)
		coeffs[1] = 1
		coeffs[0] = f.Neg(f.Reduce(x))

		polys[i] = field.NewPolynomial(f, coeffs, false)
	}

	return field.PolyProduct(e.pr, polys)
}

// does not support fast Gao.
func (e *SlowEvaluator) isNTT() bool {
	return false
}
