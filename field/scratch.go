package field

import "sync"

// This file holds every working buffer a ring's operations borrow, so that adding an
// algorithm does not mean adding another pool.
//
// One pool serves all of them, because what an operation needs is always a coefficient
// buffer and a header to pass it around by: the transform and convolution routines take
// *Polynomial, and lending the pair together saves allocating a throwaway header around a
// bare slice at every call site.
//
// The pool hangs off the ring rather than off the package: a ring is already the owner of
// what its operations need, buffers do not migrate between rings over different fields,
// and a ring that goes out of scope takes its buffers with it. It is usable from the zero
// value, so a ring built without [NewPolyRing] still works.

// scratchPools lends the polynomials borrowed during an operation and returned before it
// ends. A polynomial goes back with its array intact, so the next borrower usually finds
// one already large enough.
type scratchPools struct {
	polys sync.Pool // *Polynomial
}

// borrowPoly returns a scratch polynomial over this ring's field holding n coefficients
// of unspecified value: the caller must write every one it later reads. Use
// [PolyRing.borrowPolyZeroed] where the tail has to read as zero.
//
// The caller owns what it gets, and keeping it is an ordinary allocation. Returning it
// with [PolyRing.returnPoly] offers it to the next borrower instead, which is worth doing
// for anything that dies in the call.
func (r *PolyRing) borrowPoly(n int) *Polynomial {
	p := r.takePoly()
	ensureLenCheap(p, n)

	return p
}

// borrowPolyZeroed returns a scratch polynomial of exactly n zero coefficients, for the
// callers that copy a shorter operand in and need the rest to be padding.
func (r *PolyRing) borrowPolyZeroed(n int) *Polynomial {
	p := r.takePoly()
	p.inner = resizeZeroed(p.inner, n)

	return p
}

func (r *PolyRing) takePoly() *Polynomial {
	p, _ := r.scratch.polys.Get().(*Polynomial)
	if p == nil {
		p = &Polynomial{}
	}

	p.f = r.f
	p.isNTT = false

	return p
}

func (r *PolyRing) returnPoly(p *Polynomial) {
	r.scratch.polys.Put(p)
}
