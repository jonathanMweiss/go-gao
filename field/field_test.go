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
	a.Equal(uint64(65281), root.value)

	root, err = f.GetRootOfUnity(8)
	a.NoError(err)
	a.Equal(uint64(4096), root.value)

	f, err = NewPrimeField(157)
	a.NoError(err)

	root, err = f.GetRootOfUnity(4)
	a.NoError(err)
	a.Equal(uint64(129), root.value)
}

func TestCorrectOps(t *testing.T) {
	a := assert.New(t)

	f, err := NewPrimeField(9191248642791733759) // p > 2^62
	a.NoError(err)

	n := uint64((1 << 63) - 1)

	e1 := f.ElemFromUint64(n)

	e2 := &big.Int{}
	e2.SetUint64(n)
	e2.Mul(e2, e2)
	e2.Mod(e2, f.asBigInt)

	a.Equal(e2.Uint64(), e1.Mul(e1).Value())
	// ((1<<63 - 1)^2) %8651551326164947003

	res := e1.Mul(e1.Inverse())
	a.Equal(uint64(1), res.Value())
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

		e1 := fld.ElemFromUint64(num)
		e2 := e1.Inverse()

		res := e1.Mul(e2)
		if res.Value() != 1 {
			t.Fatalf("expected 1, got %d", res.Value())
		}

		ne1 := e1.Neg()
		if ne1.Add(e1).Value() != uint64(0) {
			t.Fatalf("expected 0, got %d", ne1.Add(e1).Value())
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

		e1 := fld.ElemFromUint64(num)
		ne1 := e1.Neg()

		if ne1.Add(e1).Value() != uint64(0) {
			t.Fatalf("expected 0, got %d", ne1.Add(e1).Value())
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

	e1 := f.ElemFromUint64((1 << 63) - 2)
	e2 := f.ElemFromUint64((1 << 60) + 312)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		e1.Mul(e2)
	}
}

func (e Elem) PowSlow(p uint64) Elem {
	if p == 0 {
		return e.field.ElemFromUint64(1)
	}

	eBigInt := &big.Int{}
	eBigInt.SetUint64(e.value)

	pBigInt := &big.Int{}
	pBigInt.SetUint64(p)

	return e.field.ElemFromUint64(
		eBigInt.Exp(eBigInt, pBigInt, e.field.asBigInt).Uint64(),
	)
}

func BenchmarkPowMod(b *testing.B) {
	f, err := NewPrimeField(9191248642791733759)
	if err != nil {
		b.FailNow()
	}

	e1 := f.ElemFromUint64((1 << 63) - 2)

	b.ResetTimer()
	b.Run("Pow", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			e1.Pow(1 << 62)
		}
	})

	b.Run("PowBig", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			e1.PowSlow(1 << 62)
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

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		b1.Mul(b1, b2)
		b1.Mod(b1, f.asBigInt)
	}
}
