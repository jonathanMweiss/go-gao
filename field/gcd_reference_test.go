// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

// The classical extended Euclidean algorithm, kept as the reference PartialGCD is
// measured and checked against.
//
// It was the original PartialGCD: it is the oracle in TestPartialGCDMatchesClassical
// and FuzzGCDAgreement, it carries the Bezout-identity property test in poly_test.go,
// and it is the quadratic baseline the half-GCD benchmarks are timed against.
// Being straightforward allows it to catch bugs and issues in the more complex half-GCD.

// partialExtendedEuclidean runs the extended Euclidean algorithm, stopping early.
//
// returns r= gcd(a,b), x, y such that ax + by = r.
// where r.Degree() < stopDegree. For full GCD, use stopDegree=0.
func (r *PolyRing) partialExtendedEuclidean(a, b *Polynomial, stopDegree int) (gcd, x, y *Polynomial) {
	// Work on local copies ensuring inputs aren't mutated.
	A := a.Copy()
	B := b.Copy()
	degA := A.Degree()
	degB := B.Degree()

	// M is the matrix that begins as identity, and each iteration it is updated by left-multiplying the Bézout matrix of the current division step:
	// M =  a00 a01 = 1 0
	//      a10 a11   0 1
	M := polyIdentity2x2(r.f)

	// Reusable temporaries (avoid allocations).
	tmp1 := &Polynomial{f: r.f} // holds q*M10 or q*M11
	tmp2 := &Polynomial{f: r.f} // holds M00 - q*M10 or M01 - q*M11

	// Bezout identity is GCD(A,B)= ax +by.
	// letsdive into the iterative algorithm.
	// Initially, we have A=a, B=b, M = I, so the invariant holds:
	// | A |   | 1  0 |   | a |
	// | B | = | 0  1 | * | b |
	// In each step of the Euclidean algorithm, we perform the division $a = q \cdot b + r$, which implies $r = a - q \cdot b$.
	// Thus,
	// | B |   | 0  1 |   | a |   | b |
	// | r | = |1  -q| *  | b | = | a - q*b |
	// Thus,
	// each step performs left-multiplication by the Bézout matrix of the current division step:
	// | 0   1 |   | M00  M01 |   | M10            M11           |
	// |1   -q | * | M10  M11 | = | M00 - q*M10    M01 - q*M11   |
	for degA >= stopDegree {
		// If B == 0, can't divide further.
		if degB < 0 {
			break
		}

		// A = q*B + r
		q, rrem := r.Div(A, B)
		A, B = B, rrem // GCD recursive step: gcd(A, B) = gcd(B,rrem)
		degA, degB = degB, B.Degree()

		// Performing the matrix update one part at a time, to reuse the same temporary
		// for both x and y updates and avoid extra allocations.

		// left matrix update: (M00, M10) = (M10, M00 - q*M10)
		r.Mul(q, M.a10, tmp1)    // tmp1 = q * M10
		r.Sub(M.a00, tmp1, tmp2) // tmp2 = M00 - q*M10
		M.a00, M.a10, tmp2 = M.a10, tmp2, M.a00

		// right update: (M01, M11) = (M11, M01 - q*M11)
		r.Mul(q, M.a11, tmp1)    // tmp1 = q * M11
		r.Sub(M.a01, tmp1, tmp2) // tmp2 = M01 - q*M11
		M.a01, M.a11, tmp2 = M.a11, tmp2, M.a01
	}

	// Notice that in each step we updated A,B using the GCD step, so now B=0, A=GCD.
	// we can extract the coeffiecnts x,y from the matrix M,
	// since we maintained the invariant that A = M.a00*a + M.a01*b.
	// Namely, A is the top-left position in our vector, which is the GCD,
	// and M00, M01 are what we multiplied against a,b to get the result A.
	return A, M.a00, M.a01
}
