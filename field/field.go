// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"errors"
	"fmt"
	"math/big"
	"math/bits"
)

// Field is the arithmetic this package computes over, an alias for the concrete
// [PrimeField]. Elements are uint64 values in [0, p).
type Field = *PrimeField

// PrimeField is the arithmetic of the integers modulo a prime p < 2^63.
// It is safe for concurrent use.
//
// Every operation returns an element reduced into [0, p).
type PrimeField struct {
	prime uint64
}

var (
	errPrimeTooLarge = errors.New("supporting up to 63-bit prime")
	errNotPrime      = errors.New("this package only supports prime fields, so the modulus must be prime")
)

const maxBitUsage = 63

// NewPrimeField returns the field of integers modulo prime.
//
// prime must be prime and below 2^63. Which sizes of NTT the field admits follows from
// the largest power of two dividing prime-1; see [RootOfUnity].
func NewPrimeField(prime uint64) (*PrimeField, error) {
	if prime > (1 << maxBitUsage) {
		return nil, errPrimeTooLarge
	}

	b := (&big.Int{}).SetUint64(prime)
	// Probably prime is 100% accurate for 64-bit numbers. Thus, we can use one base check.
	if !b.ProbablyPrime(1) {
		return nil, errNotPrime
	}

	return &PrimeField{prime: prime}, nil
}

var (
	errNotPowerOfTwo = errors.New("n must be a power of 2")
	errNotDivisible  = errors.New("n must divide p-1")
	errNSTooSmall    = errors.New("n must be >= 2")
)

// RootOfUnity returns a primitive n-th root of unity in f: an element of multiplicative
// order exactly n, which is what an n-point NTT transforms over.
//
// n must be a power of two dividing Modulus()-1, and no other n is accepted. So the
// largest power of two dividing p-1 (not the size of p) is what bounds the transform
// sizes a field can serve.
func RootOfUnity(f Field, n uint64) (uint64, error) {
	if n < 2 {
		return 0, errNSTooSmall
	}

	if !IsPowerOfTwo(n) {
		return 0, errNotPowerOfTwo
	}

	p := f.Modulus()
	if (p-1)%n != 0 {
		return 0, errNotDivisible
	}

	// The nth root of unity is the generator raised to the power of (p-1)/n
	// since g^(x) == 1 (mod p) iff x=p-1, then w=g^((p-1)/n) is not 1, and the following n powers of w != 1 too.
	// proof is by contradiction to g being the generator of the field.
	//
	// instead, of looking for g, we can use the fact that g is a generator of the field:
	// put an arbitrary `a` where `g` stood, and see what changes.
	// Every `a` is g^x for some x, which gives
	//
	//	a^((p-1)/n) = (g^x)^((p-1)/n) = w^x
	//
	// that is a^((p-1)/n) yields a candidate w^x for the primitive root.
	// this candidate is a primitive w shifted by x.
	// thus, we know that w^n = 1. However, we need to certify that w^x is not a lower order root of unity.
	// we inspect the order of w^x by testing (w^x)^(n/2) == 1. If it is, then the order is below n,
	// and we reject it. if it is not, the order can only be n: it divides n, and n being a
	// power of two the orders available are 1, 2, 4, ..., n -- every one of them below n
	// divides n/2, so any of them would have given 1 here.
	//
	// This happens when x is odd:
	// for x=2y+1 we have
	//
	// (w^x)^(n/2) = (w^(2y+1))^(n/2) = (w^n)^y * w^(n/2) = 1^y * w^(n/2) = w^(n/2) != 1
	//
	// the last step because w is primitive: no power of w below the n-th is 1.
	//
	// we want an `a` such that a=g^x for odd x. Since g is a generator f={g^0, g^1, ..., g^(p-2)},
	// and half of the elements are g^x for odd x, we can scan through the field until we find one.
	exp := (p - 1) / n

	// Scanned from 2 rather than sampled, so one field always yields one root.
	// Two codes built over the same field then agree on their evaluation order.
	for a := uint64(2); a < p; a++ {
		w := f.Pow(a, exp)
		if !f.Equals(f.Pow(w, n/2), 1) {
			return w, nil
		}
	}

	// Unreachable for a prime modulus: half the field passes the test above.
	return 0, fmt.Errorf("no primitive %d-th root of unity modulo %d", n, p)
}

// Modulus returns the field prime.
func (f *PrimeField) Modulus() uint64 {
	return f.prime
}

// IsPowerOfTwo reports whether n is a power of two.
func IsPowerOfTwo(n uint64) bool {
	// https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
	return n != 0 && (n&(n-1)) == 0
}

// Reduce returns val modulo the field prime.
func (f *PrimeField) Reduce(val uint64) uint64 {
	if val < f.prime {
		return val
	}

	return val % f.prime
}

// Add returns a + b modulo the field prime.
func (f *PrimeField) Add(a, b uint64) uint64 {
	tmp := a + b
	if tmp >= f.prime {
		tmp -= f.prime
	}

	return tmp
}

// Mul returns e * b (mod field prime).
func (f *PrimeField) Mul(a, b uint64) uint64 {
	if a == 0 || b == 0 {
		return 0
	}

	return fieldMul(a, b, f.prime)
}

func fieldMul(a, b uint64, mod uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	_, rem := bits.Div64(hi, lo, mod)

	return rem
}

// Pow returns base^exp modulo the field prime, by exponentiation by squaring.
// https://en.wikipedia.org/wiki/Exponentiation_by_squaring
func (f *PrimeField) Pow(base, exp uint64) uint64 {
	mod := f.prime

	x := uint64(1)
	for exp > 0 {
		if exp%2 == 1 { // If exponent is odd, multiply base with x
			x = fieldMul(x, base, mod)
			// x = x.Mul(base).Mod(mod)
		}

		base = fieldMul(base, base, mod) // Square the base
		exp /= 2                         // Halve the exponent
	}

	return x % mod
}

// Inverse returns the multiplicative inverse of e modulo the field prime.
//
// It panics if e is zero, which has no inverse. Guard the argument if it can be zero.
func (f *PrimeField) Inverse(e uint64) uint64 {
	// Fermat's little theorem: a^(p) = a (mod p)
	// thus:
	// a^(p-2)*a^p = a^(2p-2) = a^(p-1)^2 = 1*1=1 (mod p)
	// a^(p-2) is the inverse of a
	if e == 0 {
		panic("zero has no inverse")
	}

	return f.Pow(e, f.prime-2)
}

// Neg returns the additive inverse of e modulo the field prime.
func (f *PrimeField) Neg(e uint64) uint64 {
	res := f.prime - e
	if e == 0 {
		res = 0
	}
	return res
}

// Sub returns a - b modulo the field prime.
func (f *PrimeField) Sub(a, b uint64) uint64 {
	if a < b {
		return f.prime - (b - a)
	}

	return a - b
}

// Equals reports whether a and b are the same field element.
func (f *PrimeField) Equals(a, b uint64) bool {
	mod := f.prime
	return (a % mod) == (b % mod)
}
