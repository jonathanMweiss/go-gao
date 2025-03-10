package gao

import (
	"errors"
	"sync"

	"github.com/jonathanmweiss/go-gao/field"
)

// TODO: Implement FAST evaluation using Negacyclic NTT.

// Can be fast evaluation, like NTT, negacyclic-NTT, or just plain polynomial evaluation.
type EvaluationMap interface {
	// has access to a specific prime field.
	PrimeField() *field.PrimeField
	// returns the evaluation points for a polynomial of degree n
	EvaluationPoints(n int) (xs []uint64)
	EvaluatePolynomial(p *field.Polynomial) (ys []uint64, err error)
}

// can be fast, can be slow.

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

type SlowEvaluator struct {
	cache *evaluationCache

	field *field.PrimeField
}

func (e *evaluationCache) loadPoints(n int) []uint64 {
	e.Lock()
	defer e.Unlock()

	if points, ok := e.degreeToPoints[n]; ok {
		return points
	}

	return nil
}
func NewSlowEvaluator(field *field.PrimeField) *SlowEvaluator {
	return &SlowEvaluator{
		field: field,
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

func (e *SlowEvaluator) PrimeField() *field.PrimeField {
	return e.field
}

func (e *SlowEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if p.IsCoeffMode() {
		return nil, errNotInCoefficientForm
	}

	points := e.EvaluationPoints(len(p.ToSlice()))
	values := make([]uint64, len(points))

	for i, x := range points {
		values[i] = p.Eval(x).Value()
	}

	return values, nil
}
