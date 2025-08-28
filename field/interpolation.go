package field

import "errors"

type Interpolator struct {
	pr PolyRing
}

func NewInterpolator(pr PolyRing) *Interpolator {
	return &Interpolator{pr: pr}
}

var (
	errPointsSizeMismatch = errors.New("points size mismatch")
	errNonUniqueXs        = errors.New("non-unique x values")
)

// Interpolation code follows the Lagrange interpolation method
// https://en.wikipedia.org/wiki/Lagrange_polynomial
// This algorithm is optimise to save on operations. It is O(n^2) in total.
// The algorithm is as follows:
// 1. Create m(x) = \prod_{0\le i \le n} m_i(x) = \prod_{0\le i \le n} (x - x_i)
// 2. For each i, create q_i(x) = m(x) / m_i(x). This is done by removing m_i(x) from m(x) by dividing by m_i(x).
// 3. then from each q_i create l_i by multiplying q_i by the inverse of q_i(x_i).
// 4. Finally, sum all l_i* y_i to get the polynomial.
func (intr *Interpolator) Interpolate(xs, ys []uint64) (*Polynomial, error) {
	if err := validateInterpolationPoints(xs, ys); err != nil {
		return nil, err
	}

	// Creating m(x) = \prod_{0\le i \le n} m_i(x) = \prod_{0\le i \le n} (x - x_i)
	miSlice := intr.createMiSlice(xs)

	// O(n^2) total cost, since we are multiplying n polynomials of degree 1.
	m := PolyProduct(intr.pr, miSlice)

	liSlice := make([]Polynomial, len(xs))

	pr := intr.pr
	f := pr.GetField()
	for i, mi := range miSlice {
		qi := intr.mDivMi(m, mi) // O(n) fast division.
		s := pr.Evaluate(qi, xs[i])

		// this will be the denominator inside the product: \prod_{0\le j \le n, j\ne i} (x_i - u_j)/ (u_i-u_j)
		sinv := f.Inverse(s)

		// O(n):
		pr.MulScalar(qi, sinv, &liSlice[i])
		// liSlice[i] = qi.MulScalarInPlace(sinv) // l_i(x)
	}

	for i, li := range liSlice {
		pr.MulScalar(&li, ys[i], &li)
		// li.MulScalarInPlace(intr.pr.Reduce(ys[i]))
	}

	return intr.similarDegreePolySum(liSlice), nil
}

// PolyProduct multiplies a slice of polynomials.
func PolyProduct(pr PolyRing, miSlice []*Polynomial) *Polynomial {
	m := makeConstantPoly(pr.GetField(), 1)
	for _, mi := range miSlice {
		pr.MulPoly(m, mi, m)
	}

	return m
}

// similarDegreePolySum sums polynomials of the same degree.
func (intr *Interpolator) similarDegreePolySum(polys []Polynomial) *Polynomial {
	inner := make([]uint64, len(polys[0].inner))
	fld := intr.pr.GetField()
	for _, poly := range polys {
		for i, coef := range poly.inner {
			inner[i] = fld.Add(inner[i], coef)
		}
	}

	return NewPolynomial(fld, inner, false)

}

// createMiSlice creates the m_i(x) = (x - x_i) polynomials.
func (intr *Interpolator) createMiSlice(xs []uint64) []*Polynomial {
	miSlice := make([]*Polynomial, len(xs))
	f := intr.pr.GetField()
	for i, x := range xs {
		miInner := make([]uint64, 2)

		miInner[0] = f.Neg(f.Reduce(x))
		miInner[1] = 1

		miSlice[i] = NewPolynomial(f, miInner, false)
	}

	return miSlice
}

/*
mDivMi divides m by mi. This is quicker than the long division method since
we know that mi is of degree 1, and that we don't have a remainder.
*/
func (intr *Interpolator) mDivMi(m_, mi_ *Polynomial) *Polynomial {
	m := m_.Copy()
	qinner := make([]uint64, len(m.inner)-1)
	ui := mi_.inner[0]

	f := intr.pr.GetField()
	tmp := uint64(0)

	for i := len(m.inner) - 1; i > 0; i-- {
		qinner[i-1] = m.inner[i]
		// take m_i = x - ui
		// remove ui from m:

		tmp = f.Neg(f.Mul(m.inner[i], ui))
		m.inner[i-1] = f.Add(tmp, m.inner[i-1])
	}

	return NewPolynomial(f, qinner, false)
}

func validateInterpolationPoints(xs []uint64, ys []uint64) error {
	if len(xs) != len(ys) {
		return errPointsSizeMismatch
	}

	mapXs := make(map[uint64]struct{})
	for _, x := range xs {
		mapXs[x] = struct{}{}
	}

	if len(mapXs) != len(xs) {
		return errNonUniqueXs
	}

	return nil
}
