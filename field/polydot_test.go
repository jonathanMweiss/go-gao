package field

import (
	"fmt"
	"testing"
)

// dotShapes spans both multiplication paths: the small entries stay on schoolbook, the
// lopsided ones reach the blocked convolution.
var dotShapes = []struct{ lx, lp, ly, lq int }{
	{1, 1, 1, 1},
	{4, 9, 4, 9},
	{9, 4, 4, 9},
	{163, 16384, 163, 16384},
	{163, 16384, 164, 16383},
	{164, 16384, 163, 512},
	{2048, 2048, 2048, 2048},
}

// TestDotNewMatchesNaive checks x*p + y*q against the three-buffer formulation it
// replaces, for shapes on either side of the multiplication dispatch.
func TestDotNewMatchesNaive(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	for _, s := range dotShapes {
		t.Run(fmt.Sprintf("%dx%d+%dx%d", s.lx, s.lp, s.ly, s.lq), func(t *testing.T) {
			x, p := rampPoly(r, s.lx, 1), rampPoly(r, s.lp, 2)
			y, q := rampPoly(r, s.ly, 3), rampPoly(r, s.lq, 4)

			want := r.addNew(r.mulNew(x, p), r.mulNew(y, q))
			got := r.dotNew(x, p, y, q)

			if !got.Equals(want) {
				t.Fatalf("length %d vs %d", len(got.inner), len(want.inner))
			}
		})
	}
}

// TestDotNewDoesNotMutateOperands guards the borrowed buffer against being handed out as
// one of the operands, which would corrupt a caller's polynomial.
func TestDotNewDoesNotMutateOperands(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	x, p := rampPoly(r, 163, 1), rampPoly(r, 16384, 2)
	y, q := rampPoly(r, 163, 3), rampPoly(r, 16384, 4)

	xc, pc, yc, qc := x.Copy(), p.Copy(), y.Copy(), q.Copy()

	r.dotNew(x, p, y, q)

	for _, pair := range []struct {
		name      string
		got, want *Polynomial
	}{{"x", x, xc}, {"p", p, pc}, {"y", y, yc}, {"q", q, qc}} {
		if !pair.got.Equals(pair.want) {
			t.Errorf("dotNew modified operand %s", pair.name)
		}
	}
}

// TestDotNewPoolReuse alternates a large shape with a small one so the borrowed buffer
// comes back longer than the next caller needs. A product that failed to define every
// coefficient it returns would surface here as leftovers from the previous call.
func TestDotNewPoolReuse(t *testing.T) {
	r := ringOver(t, NTTFriendlyPrime)

	big := []int{163, 16384}
	small := []int{3, 7}

	for i := range 8 {
		dims := small
		if i%2 == 0 {
			dims = big
		}

		x, p := rampPoly(r, dims[0], 1), rampPoly(r, dims[1], 2)
		y, q := rampPoly(r, dims[0], 3), rampPoly(r, dims[1], 4)

		want := r.addNew(r.mulNew(x, p), r.mulNew(y, q))
		if got := r.dotNew(x, p, y, q); !got.Equals(want) {
			t.Fatalf("iteration %d (%dx%d): result disagrees with the naive form", i, dims[0], dims[1])
		}
	}
}
