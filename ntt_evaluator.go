package gao

import (
	"github.com/jonathanmweiss/go-gao/field"
)

type NttEvaluator struct {
	cache *evaluationCache

	pr field.PolyRing
}

func NewNttEvaluator(f field.Field) *NttEvaluator {
	return &NttEvaluator{
		pr:    field.NewDensePolyRing(f),
		cache: newEvaluatorCache(),
	}
}

func (e *NttEvaluator) EvaluationPoints(n int) []uint64 {
	points := e.cache.loadPoints(n)
	if points != nil {
		return points
	}

	// make polynomial p(x) = x.
	// then attempt to compute its NTT.
	inner := make([]uint64, n)
	inner[1] = 1
	p := field.NewPolynomial(e.pr.GetField(), inner, false)

	if err := e.pr.NttForward(p); err != nil {
		panic(err) //. TODO: change API.
	}

	points = p.NoCopySlice()

	e.cache.storePoints(n, points)

	return points
}

func (e *NttEvaluator) PrimeField() field.Field {
	return e.pr.GetField()
}

func (e *NttEvaluator) EvaluatePolynomial(p *field.Polynomial) ([]uint64, error) {
	if err := e.pr.NttForward(p); err != nil {
		return nil, err
	}

	return p.NoCopySlice(), nil
}

func (e *NttEvaluator) GenerateLocatorPolynomial(n int) *field.Polynomial {
	// The locator polynomial L(x) = (x - x_1)(x - x_2)...(x - x_n)
	// where x_1, x_2, ..., x_n are the evaluation points
	// is vanishing for the roots of unity: L(x)=1*x^n-1
	f := e.pr.GetField()
	inner := make([]uint64, n+1)
	inner[0] = 1
	inner[n] = f.Neg(1)
	return field.NewPolynomial(f, inner, false)
}

// does not support fast Gao.
func (e *NttEvaluator) isNTT() bool {
	return true
}
