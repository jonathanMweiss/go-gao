package field

import "errors"

type Interpolator struct {
	Field *PrimeField
}

func NewInterpolator(field *PrimeField) *Interpolator {
	return &Interpolator{Field: field}
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
	m := intr.Field.PolyProduct(miSlice)

	liSlice := make([]*Polynomial, len(xs))

	for i, mi := range miSlice {
		qi := intr.mDivMi(m, mi) // O(n) fast division.
		s := qi.Eval(xs[i])

		// this will be the denominator inside the product: \prod_{0\le j \le n, j\ne i} (x_i - u_j)/ (u_i-u_j)
		sinv := s.Inverse()

		// O(n):
		liSlice[i] = qi.MulScalarInPlace(sinv) // l_i(x)
	}

	for i, li := range liSlice {
		li.MulScalarInPlace(intr.Field.ElemFromUint64(ys[i]))
	}

	return intr.similarDegreePolySum(liSlice), nil
}

// PolyProduct multiplies a slice of polynomials.
func (f *PrimeField) PolyProduct(miSlice []*Polynomial) *Polynomial {
	m := f.constantPolynomial(1)
	for _, mi := range miSlice {
		m = m.Mul(mi)
	}

	return m
}

// similarDegreePolySum sums polynomials of the same degree.
func (intr *Interpolator) similarDegreePolySum(polys []*Polynomial) *Polynomial {
	inner := make([]Elem, len(polys[0].inner))
	for _, poly := range polys {
		for i, coef := range poly.inner {
			inner[i] = inner[i].Add(coef)
		}
	}

	return NewPolynomial(inner, false)

}

// createMiSlice creates the m_i(x) = (x - x_i) polynomials.
func (intr *Interpolator) createMiSlice(xs []uint64) []*Polynomial {
	miSlice := make([]*Polynomial, len(xs))

	for i, x := range xs {
		miInner := make([]Elem, 2)

		miInner[0] = intr.Field.ElemFromUint64(x).Neg()
		miInner[1] = intr.Field.ElemFromUint64(1)

		miSlice[i] = NewPolynomial(miInner, false)
	}

	return miSlice
}

/*
mDivMi divides m by mi. This is quicker than the long division method since
we know that mi is of degree 1, and that we don't have a remainder.
*/
func (intr *Interpolator) mDivMi(m_, mi_ *Polynomial) *Polynomial {
	m := m_.Copy()
	qinner := make([]Elem, len(m.inner)-1)
	ui := mi_.inner[0]

	for i := len(m.inner) - 1; i > 0; i-- {
		qinner[i-1] = m.inner[i]
		// take m_i = x - ui
		// remove ui from m:
		tmp := m.inner[i].Mul(ui).Neg()
		m.inner[i-1] = tmp.Add(m.inner[i-1])
	}

	return NewPolynomial(qinner, false)
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
