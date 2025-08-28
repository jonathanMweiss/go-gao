package field

import (
	"errors"
	"math/big"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

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
}

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

/*
Assumes you are using a prime. Will not check for validity.
*/
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

	bgint := &big.Int{}
	bgint.SetUint64(prime)

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

func (f *PrimeField) ElemSlice(vals []uint64) []uint64 {
	mod := f.prime
	for i, v := range vals {
		vals[i] = v % mod
	}

	return vals
}

func IsPowerOfTwo(n uint64) bool {
	// https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
	return n != 0 && (n&(n-1)) == 0
}

func (f *PrimeField) Prime() uint64 {
	return f.prime
}

func (f *PrimeField) Generator() uint64 {
	return f.generator
}

func (f *PrimeField) Factors() []uint64 {
	return f.factors
}

func (f *PrimeField) constantPolynomial(val uint64) *Polynomial {
	return &Polynomial{
		inner:            []uint64{f.Reduce(val)},
		isCoefficientMod: false,
		f:                f,
	}
}

func (f *PrimeField) Reduce(val uint64) uint64 {
	return val % f.prime
}

func (f *PrimeField) Add(a, b uint64) uint64 {
	if a == 0 {
		// used for Elem{}, or Elem{0,nilField}.
		return b
	}

	tmp := a + b // can't overflow since adding two integers smaller than 2^63.
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

func (f *PrimeField) Neg(e uint64) uint64 {
	if e == 0 {
		return 0
	}

	return (f.prime - e)
}

func (f *PrimeField) Sub(a, b uint64) uint64 {
	if a < b {
		return f.prime - (b - a)
	}

	return a - b
}

func (f *PrimeField) Equals(a, b uint64) bool {
	mod := f.prime
	return (a % mod) == (b % mod)
}
