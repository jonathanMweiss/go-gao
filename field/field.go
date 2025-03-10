package field

import (
	"errors"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/ring"
	bigint "lukechampine.com/uint128"
)

type PrimeField struct {
	prime     uint64
	generator uint64
	factors   []uint64
	asBigInt  *big.Int
}

var (
	errPrimeTooLarge = errors.New("supporting up to 63-bit prime")
	errNotPrime      = errors.New("this package only support prime fields. please use a prime order")
)

const maxBitUsage = 63

/*
Assumes you are using a prime. Will not check for validity.
*/
func NewPrimeField(prime uint64) (*PrimeField, error) {
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
		asBigInt:  bgint,
	}, nil
}

var (
	errNotPowerOfTwo = errors.New("n must be a power of 2")
	errNotDivisible  = errors.New("n must divide p-1")
	errNSTooSmall    = errors.New("n must be >= 2")
)

func (f *PrimeField) GetRootOfUnity(n uint64) (Elem, error) {
	if n == 0 || n == 1 {
		return Elem{}, errNSTooSmall
	}

	if !isPowerOfTwo(n) {
		return Elem{}, errNotPowerOfTwo
	}

	if bigint.From64(f.prime-1).Mod64(n) != 0 {
		return Elem{}, errNotDivisible
	}

	// The nth root of unity is the generator raised to the power of (prime-1)/n
	// since g^(x) == 1 (mod p) iff x=p-1, then w=g^((p-1)/n) is not 1, and the following n powers of w != 1 too.
	// proof is by contradiction to g being the generator of the field.
	return f.ElemFromUint64(f.generator).Pow((f.prime - 1) / n), nil
}

func (f *PrimeField) ElemSlice(vals []uint64) []Elem {
	elems := make([]Elem, len(vals))

	for i, v := range vals {
		elems[i] = f.ElemFromUint64(v)
	}

	return elems
}

func isPowerOfTwo(n uint64) bool {
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
		inner:            []Elem{f.ElemFromUint64(val)},
		isCoefficientMod: false,
		f:                f,
	}
}

type Elem struct {
	field *PrimeField
	value uint64
}

func (f *PrimeField) ElemFromUint64(val uint64) Elem {
	return Elem{
		field: f,
		value: val % f.prime,
	}
}

func (e Elem) Value() uint64 {
	return e.value % e.field.prime
}

func (e Elem) Copy() Elem {
	return Elem{
		field: e.field,
		value: e.value,
	}
}

func (e Elem) IsZero() bool {
	return e.value == 0
}

func (e Elem) Add(b Elem) Elem {
	if e.IsZero() {
		return b // empty elem is 0.
	}

	if e.field.prime != b.field.prime {
		panic("prime mismatch")
	}

	return e.field.ElemFromUint64(e.value + b.value)
}

func (e Elem) Mul(b Elem) Elem {
	if e.IsZero() || b.IsZero() {
		return Elem{b.field, 0} // empty elem is 0
	}

	if e.field.prime != b.field.prime {
		panic("prime mismatch")
	}

	n := bigint.From64(e.value).Mul64(b.value)

	return e.field.ElemFromUint64(n.Mod64(e.field.prime))
}

// https://en.wikipedia.org/wiki/Exponentiation_by_squaring
func (e Elem) Pow(exp uint64) Elem {
	mod := bigint.From64(e.field.prime)
	base := bigint.From64(e.value % e.field.prime)

	x := bigint.From64(1)

	for exp > 0 {
		if exp%2 == 1 { // If exponent is odd, multiply base with x
			x = x.Mul(base).Mod(mod)
		}

		base = base.Mul(base).Mod(mod) // Square the base
		exp /= 2                       // Halve the exponent
	}

	return e.field.ElemFromUint64(x.Mod64(e.field.prime))
}

func (e Elem) Inverse() Elem {
	// Fermat's little theorem: a^(p) = a (mod p)
	// thus:
	// a^(p-2)*a^p = a^(2p-2) = a^(p-1)^2 = 1*1=1 (mod p)
	// a^(p-2) is the inverse of a
	if e.IsZero() {
		panic("zero has no inverse")
	}

	return e.Pow(e.field.prime - 2)
}

func (e Elem) Neg() Elem {
	return e.field.ElemFromUint64(e.field.prime - e.value)
}

func (e Elem) Sub(b Elem) Elem {
	return e.Add(b.Neg())
}

func (e Elem) Equals(o Elem) bool {
	return e.field == o.field && e.value == o.value
}
