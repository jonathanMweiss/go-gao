// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

/*
Package field provides prime-field and polynomial arithmetic over uint64.

It exists to support the Reed-Solomon decoder in the parent gao package, but
stands on its own for any work over a prime field: NTT-based multiplication,
division with remainder, a half-GCD extended Euclidean algorithm, and
multipoint interpolation. It depends on nothing outside the standard library.

# Fields

[NewPrimeField] builds a [PrimeField] for a prime modulus below 2^63. Elements
are plain uint64 values in [0, p), arithmetic is modular, and [PrimeField.Reduce]
normalizes a value that may fall outside the range.

	f, err := field.NewPrimeField(65537)
	if err != nil {
		return err
	}
	c := f.Mul(f.Add(2, 3), f.Inverse(7))

[Field] names that arithmetic. It is an alias for the concrete [PrimeField].

[RootOfUnity] supplies that structure instead, deriving a primitive n-th root
of unity from the field. It accepts only orders that are powers of two
dividing p-1, because those are the only sizes this package transforms at. The
largest such order — the largest power of two dividing p-1, not the size of p
— is what bounds the code lengths a field can serve, and it is worth checking
before choosing a prime: 65537 admits a 65536-point transform, while 929
admits only a 32-point one.

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

[NewPolyRing] returns a [PolyRing], which carries the polynomial arithmetic
and caches the NTT twiddle factors it derives, so one ring is worth reusing
across many operations over the same field.

Its methods write into a caller-supplied destination polynomial rather than
allocating, and several dispatch on size: [PolyRing.Mul] and [PolyRing.Div]
choose between schoolbook and NTT-based algorithms according to operand size
and whether the field admits the roots of unity they need. [PolyRing.PartialGCD]
chooses likewise between a classical quadratic Euclidean sequence and the
half-GCD recursion, which is asymptotically faster but carries enough constant
overhead to lose on small operands. A caller picks a ring, never an algorithm.

[PolyRing.Product] multiplies a whole slice of polynomials through a
divide-and-conquer product tree, which is how locator polynomials of the form
(x-x_1)...(x-x_n) are built.

# Panics

Arithmetic on malformed input panics rather than returning an error: a nil
polynomial, a mismatched field, a polynomial in the wrong representation, or
division by the zero polynomial. These are the same conditions Go's own
arithmetic panics on — integer division by zero, and math/big alike — and a
caller that can reach one is expected to say what it means before the
operation rather than recover from it afterwards.

# Interpolation

[NewInterpolator] returns an [Interpolator], which recovers the unique
polynomial of degree < n passing through n given (x, y) pairs.
*/
package field
