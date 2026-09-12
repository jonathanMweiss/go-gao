// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

/*
Package field provides prime-field and polynomial arithmetic over uint64.

It exists to support the Reed-Solomon decoder in the parent gao package, but
stands on its own for any work over a prime field: NTT-based multiplication,
division with remainder, a half-GCD extended Euclidean algorithm, and
multipoint interpolation.

# Fields

[NewPrimeField] builds a [Field] for a prime modulus below 2^63. Elements are
plain uint64 values in [0, p), arithmetic is modular, and [Field.Reduce]
normalizes a value that may fall outside the range. Beyond the ring
operations, a Field exposes its [Field.Generator], the prime [Field.Factors]
of p-1, and [Field.GetRootOfUnity], which reports an error when the requested
order admits no root of unity — that is, unless the order is a power of two
dividing p-1.

	f, err := field.NewPrimeField(65537)
	if err != nil {
		return err
	}
	c := f.Mul(f.Add(2, 3), f.Inverse(7))

# Polynomials

A [Polynomial] holds its coefficients in one of two representations:
coefficient form, or point-value (NTT) form. Most operations require
coefficient form; [PolyRing.NttForward] and [PolyRing.NttBackward] convert
between them, in place, and [Polynomial.IsCoeffMode] reports which form a
value is in.

[PolyRing.NewPolynomial] does not copy the slice it is given — the polynomial
aliases that memory, and in-place operations such as the NTT conversions will
overwrite it. Pass a copy when the caller still needs the original.
[Polynomial.ToSlice] returns a copy of the coefficients;
[Polynomial.NoCopySlice] returns the backing array itself.

# Rings

[NewDensePolyRing] returns a [PolyRing], the interface carrying the
arithmetic. Its methods write into a caller-supplied destination polynomial
rather than allocating, and several dispatch on size: [PolyRing.Mul] and
[PolyRing.Div] choose between schoolbook and NTT-based algorithms according to
operand size and whether the field supports the roots of unity they need.

[PolyRing.Product] multiplies a whole slice of polynomials through a
divide-and-conquer product tree, which is how locator polynomials of the form
(x-x_1)...(x-x_n) are built.

For extended GCD there are two entry points with the same signature.
[PolyRing.PartialExtendedEuclidean] is the classical quadratic algorithm;
[PolyRing.FastPartialGCD] is the half-GCD variant, asymptotically faster and
the better choice at decoder-sized inputs, but with enough constant overhead
that it loses on small operands. Both stop early once the remainder drops
below stopDegree, and both leave their inputs unmodified.

Arithmetic on malformed input panics rather than returning an error: a nil
polynomial, a mismatched field, a polynomial in the wrong representation, or
division by the zero polynomial are all programming mistakes, not conditions a
caller is expected to recover from.

# Interpolation

[NewInterpolator] returns an [Interpolator], which recovers the unique
polynomial of degree < n passing through n given (x, y) pairs.
*/
package field
