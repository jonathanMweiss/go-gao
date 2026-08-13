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

	codeword, err := code.EncodeToSlice(data)
	if err != nil {
		panic(err)
	}

	// Corrupt 6 symbols -- the maximum this code can repair.
	for _, i := range []int{0, 3, 5, 9, 11, 14} {
		codeword[i] = 12345
	}

	decoded, err := code.DecodeFromSlice(codeword)
	if err != nil {
		panic(err)
	}

	fmt.Println(decoded)
	fmt.Println("max repairable errors:", code.MaxErrors())
	// Output:
	// [10 20 30 40]
	// max repairable errors: 6
}
