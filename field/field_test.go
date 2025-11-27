package field

import (
	"math/big"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestRootsOfUnity(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(65537)
	a.NoError(err)

	root, err := f.GetRootOfUnity(4)
	a.NoError(err)
	a.Equal(uint64(65281), root)

	root, err = f.GetRootOfUnity(8)
	a.NoError(err)
	a.Equal(uint64(4096), root)

	f, err = NewPrimeField(157)
	a.NoError(err)

	root, err = f.GetRootOfUnity(4)
	a.NoError(err)
	a.Equal(uint64(129), root)
}

func TestCorrectOps(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(9191248642791733759) // p > 2^62
	a.NoError(err)

	n := uint64((1 << 63) - 1)

	e1 := f.Reduce(n)

	e2 := &big.Int{}
	e2.SetUint64(n)
	e2.Mul(e2, e2)
	e2.Mod(e2, new(big.Int).SetUint64(f.Modulus()))

	a.Equal(e2.Uint64(), f.Mul(e1, e1))
	// ((1<<63 - 1)^2) %8651551326164947003

	res := f.Mul(e1, f.Inverse(e1))
	a.Equal(uint64(1), res)
}

func FuzzInverse(f *testing.F) {
	testcases := []uint64{1, 54347, 4534523, 021310, 1<<63 - 1}
	for _, tc := range testcases {
		f.Add(tc) // Use f.Add to provide a seed corpus
	}

	fld, err := NewPrimeField(9191248642791733759)
	if err != nil {
		f.FailNow()
	}

	f.Fuzz(func(t *testing.T, num uint64) {

		e1 := fld.Reduce(num)
		e2 := fld.Inverse(e1)

		res := fld.Mul(e1, e2)
		if res != 1 {
			t.Fatalf("expected 1, got %d", res)
		}

		ne1 := fld.Neg(e1)
		if fld.Add(ne1, e1) != uint64(0) {
			t.Fatalf("expected 0, got %d", fld.Add(ne1, e1))
		}
	})
}

func FuzzNegate(f *testing.F) {
	testcases := []uint64{1, 54347, 4534523, 021310, 1<<63 - 1}
	for _, tc := range testcases {
		f.Add(tc) // Use f.Add to provide a seed corpus
	}

	fld, err := NewPrimeField(9191248642791733759)
	if err != nil {
		f.FailNow()
	}

	f.Fuzz(func(t *testing.T, num uint64) {

		e1 := fld.Reduce(num)
		ne1 := fld.Neg(e1)

		if fld.Add(ne1, e1) != uint64(0) {
			t.Fatalf("expected 0, got %d", fld.Add(ne1, e1))
		}
	})
}

func BenchmarkMulMod(b *testing.B) {
	//BigMul:  42.08 ns/op
	//elemMul: 4.949 ns/op
	f, err := NewPrimeField(9191248642791733759)
	if err != nil {
		b.FailNow()
	}

	e1 := f.Reduce((1 << 63) - 2)
	e2 := f.Reduce((1 << 60) + 312)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		f.Mul(e1, e2)
	}
}

func (f *PrimeField) PowSlow(e, p uint64) uint64 {
	if p == 0 {
		return 1
	}

	eBigInt := &big.Int{}
	eBigInt.SetUint64(e)

	pBigInt := &big.Int{}
	pBigInt.SetUint64(p)

	modAsBig := new(big.Int).SetUint64(f.Modulus())
	return f.Reduce(
		eBigInt.Exp(eBigInt, pBigInt, modAsBig).Uint64(),
	)
}

func BenchmarkPowMod(b *testing.B) {
	f, err := NewPrimeField(9191248642791733759)
	if err != nil {
		b.FailNow()
	}

	e1 := f.Reduce((1 << 63) - 2)

	b.ResetTimer()
	b.Run("Pow", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			f.Pow(e1, 1<<62)
		}
	})

	fp, ok := f.(*PrimeField)
	if !ok {
		b.FailNow()
	}

	b.Run("PowBig", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			fp.PowSlow(e1, 1<<62)
		}
	})

}

func BenchmarkMulModBig(b *testing.B) {
	f, err := NewPrimeField(9191248642791733759)
	if err != nil {
		b.FailNow()
	}

	b1 := big.NewInt((1 << 63) - 2)
	b2 := big.NewInt((1 << 60) + 312)

	fAsBigint := new(big.Int).SetUint64(f.Modulus())
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		b1.Mul(b1, b2)
		b1.Mod(b1, fAsBigint)
	}
}

func FuzzSub_Simple(fz *testing.F) {
	// Use a single, fixed 63-bit prime (2^61 - 1, a known Mersenne prime).
	const p = uint64(157)

	pf, err := NewPrimeField(p)
	if err != nil {
		fz.Fatalf("failed to init PrimeField: %v", err)
	}

	// A couple of tiny seeds; the fuzzer will generate the rest.
	fz.Add(uint64(0), uint64(0))
	fz.Add(uint64(1), uint64(2))
	fz.Add(p-1, p-1)

	fz.Fuzz(func(t *testing.T, aSeed, bSeed uint64) {
		// Reduce inputs into field domain
		a := pf.Reduce(aSeed)
		b := pf.Reduce(bSeed)

		got := pf.Sub(a, b)
		want := pf.Add(a, pf.Neg(b))

		if got != want {
			t.Fatalf("Sub mismatch: got=%d, want=%d (a=%d, b=%d, p=%d)", got, want, a, b, p)
		}
	})
}

func TestRootsOfUnityGeneration(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(65537) // many roots of unity...
	a.NoError(err)

	for i := range 8 {
		N := uint64(1 << (i + 1))
		root, err := f.GetRootOfUnity(N)
		a.NoError(err)
		a.True(isRootOfUnityOfOrderN(f, root, N))
	}
}

func isRootOfUnityOfOrderN(field Field, root, n uint64) bool {
	mp := make(map[uint64]int)
	for i := uint64(0); i < n; i++ {
		tmp := field.Pow(root, i)
		mp[tmp]++
	}
	// Check if all powers are distinct
	return len(mp) == int(n) && mp[1] == 1
}
