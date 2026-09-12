// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package field

import (
	"errors"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// Field is the arithmetic of a finite field, with elements represented as uint64.
type Field interface {
	Equals(a, b uint64) bool
	Add(a, b uint64) uint64
	Sub(a, b uint64) uint64
	Mul(a, b uint64) uint64
	Pow(base, exp uint64) uint64

	Neg(a uint64) uint64
	Inverse(a uint64) uint64
	Reduce(a uint64) uint64

	Modulus() uint64
	GetRootOfUnity(n uint64) (uint64, error)
	Generator() uint64
	Factors() []uint64
}

// PrimeField implements Field over the integers modulo a prime p < 2^63.
type PrimeField struct {
	prime     uint64
	generator uint64
	factors   []uint64
}

var (
	errPrimeTooLarge = errors.New("supporting up to 63-bit prime")
	errNotPrime      = errors.New("this package only support prime fields. please use a prime order")
)

const maxBitUsage = 63

// NewPrimeField returns the field of integers modulo prime.
//
// prime must be prime and below 2^63.
func NewPrimeField(prime uint64) (Field, error) {
	if prime > (1 << maxBitUsage) {
		return nil, errPrimeTooLarge
	}

	b := (&big.Int{}).SetUint64(prime)
	// Probably prime is 100% accurate for 64-bit numbers. Thus, we can use one base check.
	if !b.ProbablyPrime(1) {
		return nil, errNotPrime
	}

	// TODO: write my own function to find a primitive root, thus dropping the dependency on lattigo altogether.
	g, factors, err := ring.PrimitiveRoot(prime, nil)
	if err != nil {
		return nil, err
	}

	return &PrimeField{
		prime:     prime,
		generator: g,
		factors:   factors,
	}, nil
}

var (
	errNotPowerOfTwo = errors.New("n must be a power of 2")
	errNotDivisible  = errors.New("n must divide p-1")
	errNSTooSmall    = errors.New("n must be >= 2")
)

// Modulus implements Field.
func (f *PrimeField) Modulus() uint64 {
	return f.prime
}

// GetRootOfUnity returns a primitive n-th root of unity, which exists only when n is
// a power of two dividing p-1.
func (f *PrimeField) GetRootOfUnity(n uint64) (uint64, error) {
	if n == 0 || n == 1 {
		return 0, errNSTooSmall
	}

	if !IsPowerOfTwo(n) {
		return 0, errNotPowerOfTwo
	}

	if (f.prime-1)%n != 0 {
		return 0, errNotDivisible
	}

	// The nth root of unity is the generator raised to the power of (prime-1)/n
	// since g^(x) == 1 (mod p) iff x=p-1, then w=g^((p-1)/n) is not 1, and the following n powers of w != 1 too.
	// proof is by contradiction to g being the generator of the field.
	return f.Pow(f.generator, (f.prime-1)/n), nil

}

// IsPowerOfTwo reports whether n is a power of two.
func IsPowerOfTwo(n uint64) bool {
	// https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
	return n != 0 && (n&(n-1)) == 0
}

// Generator returns a primitive root of the field: F_p^* = {1, g, g^2, ..., g^(p-2)}.
func (f *PrimeField) Generator() uint64 {
	return f.generator
}

// Factors returns the prime factorization of p-1.
func (f *PrimeField) Factors() []uint64 {
	return f.factors
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
