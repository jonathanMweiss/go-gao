// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao_test

import (
	"fmt"

	"github.com/jonathanmweiss/go-gao"
	"github.com/jonathanmweiss/go-gao/field"
)

// Example shows a full encode/corrupt/decode round trip: a codeword of n=16
// symbols carries k=4 data symbols, which tolerates up to (n-k)/2 = 6 corruptions.
func Example() {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		panic(err)
	}

	const n, k = 16, 4

	// RequireNTT: fail rather than silently fall back to the quadratic path.
	code, err := gao.NewCode(f, n, k, gao.RequireNTT())
	if err != nil {
		panic(err)
	}

	data := []uint64{10, 20, 30, 40}

	codeword, err := code.Encode(data)
	if err != nil {
		panic(err)
	}

	// Corrupt 6 symbols -- the maximum this code can repair.
	for _, i := range []int{0, 3, 5, 9, 11, 14} {
		codeword[i] = 12345
	}

	decoded, err := code.Decode(codeword, gao.ErasureSet{})
	if err != nil {
		panic(err)
	}

	fmt.Println(decoded)
	fmt.Println("max repairable errors:", code.MaxErrors())
	// Output:
	// [10 20 30 40]
	// max repairable errors: 6
}

// ExampleByteCode shows the byte view of a code: encode a payload, lose a run of wire
// bytes, and decode what is left. Over p=65537 a symbol carries 2 payload bytes and
// occupies 3 on the wire, so k=4 symbols take 8 bytes of payload into a 48-byte codeword.
func ExampleByteCode() {
	f, err := field.NewPrimeField(65537)
	if err != nil {
		panic(err)
	}

	const n, k = 16, 4

	code, err := gao.NewCode(f, n, k, gao.RequireNTT())
	if err != nil {
		panic(err)
	}

	bc := code.Bytes()

	payload := []byte("attack")

	raw, err := bc.Encode(payload)
	if err != nil {
		panic(err)
	}

	// Nine wire bytes are lost, which erases the three symbols they touch.
	lost, err := bc.Erasures(gao.ByteRange{Off: 6, Len: 9})
	if err != nil {
		panic(err)
	}

	for i := 6; i < 15; i++ {
		raw[i] = 0xFF
	}

	got, err := bc.Decode(raw, lost)
	if err != nil {
		panic(err)
	}

	fmt.Println("max bytes:", bc.MaxBytes())
	fmt.Println("wire size:", len(raw))
	fmt.Println("erased symbols:", lost.Len())
	fmt.Printf("%q\n", got[:len(payload)])
	// Output:
	// max bytes: 8
	// wire size: 48
	// erased symbols: 3
	// "attack"
}
