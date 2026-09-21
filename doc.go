// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

/*
Package gao implements Reed-Solomon error correction using Gao's decoder.

It repairs symbols that are silently wrong at positions you do not know: the
decoder finds the error locations itself. Where you do know the positions,
declare them as erasures, which cost half as much budget and decode faster.

In this library, symbols are elements of a prime field up to 63 bits wide, chosen by the
caller.

# Codes

The library works at two levels:

  - [ByteCode], the higher one, when you care about bytes only.
  - [Code], the lower one, when you want the symbols (field elements)
    directly.

If you are unsure, start with ByteCode and drop to Code when you need the
symbols.

Both are built over a prime field and fix two lengths: n symbols to a
codeword, k of them data. The difference n-k is the redundancy, so a larger n
buys budget and a larger k spends it on payload.

Building a Code:

	f, err := field.NewPrimeField(field.NTTFriendlyPrime) // or any other prime
	if err != nil {
		return err
	}

	// here n=16, k=4, but you can also try with
	// other values. choose powers of two for high performance with
	// FFT-like algorithms (NTT).
	code, err := gao.NewCode(f, 16, 4)
	if err != nil {
		return err
	}

	codeword, err := code.Encode([]uint64{10, 20, 30, 40})

[NewCode] reports bad parameters immediately, as [ErrNonPositiveK],
[ErrNSmallerThanK] or [ErrUnsupportedSize], rather than failing later.

A ByteCode is a view of a Code, so it starts from one:

	bc := code.Bytes()

	encoded, err := bc.Encode([]byte("My message"))

A symbol carries whole bytes of payload and occupies whole bytes on the wire,
both sized against the modulus. Over the 57-bit NTTFriendlyPrime that is 7
payload bytes in an 8-byte symbol, so the code above takes at most 28 bytes
([ByteCode.MaxBytes], k times 7) and produces 128 (n times 8). A short payload
is zero-padded, and the padding is indistinguishable from payload afterwards:
carry the original length and slice the result.

# The decoding budget

A corrupted symbol at an unknown position is an error; one at a position the
caller knows, and declares, is an erasure. Erasures are half the price,
because the decoder does not have to spend budget locating them. Decoding
succeeds while

	2*errors + erasures <= n-k

[Code.MaxErrors] reports (n-k)/2, the all-errors corner of that budget.

Past the budget [Code.Decode] usually returns [ErrDecoding], but it cannot
always tell: with enough errors a received word lands closer to a different
valid codeword, and the decoder returns a confidently wrong message. That is
inherent to Reed-Solomon codes, not to this implementation.

# Codewords

Codewords are positional, and reordering one makes it invalid. [Code.Encode]
returns n values, where index i is the evaluation at [Code.EvaluationPoints]
index i, and [Code.Decode] expects them back in that order. Store or split
them as you like, and reassemble them in the same order before decoding:

	codeword, err := code.Encode([]uint64{10, 20, 30, 40})

	// codeword[i] goes to disk i, and comes back at index i.
	msg, err := code.Decode(codeword, gao.ErasureSet{})

A codeword of any other length is rejected with [ErrMismatchedLengths], since
no amount of correction fixes framing. Decode returns a message of exactly
length k, zero-padded when the recovered message has high-order zero symbols.

# Erasures

As stated above, an erasure is a position the caller knows is unusable:
a missing symbol, or a range of bytes missing from a ByteCode codeword.
Whatever the codeword holds there is ignored,
so there is no need to blank it first, and the zero [ErasureSet] declares
nothing missing.

Building the set is the expensive half of an erasure decode and depends on the
positions alone, so words that lost the same positions should share the ErasureSet:

	// symbols 3 and 7 are unusable.
	lost, err := code.Erasures(3, 7)
	for _, word := range words {
		msg, err := code.Decode(word, lost)
	}

[ByteCode.Erasures] takes byte ranges instead and builds the same set. A
symbol any range touches is erased whole, and ranges may overlap, repeat, or
fall partly outside the codeword:

	lost, err := bc.Erasures(
		gao.ByteRange{Off: 6, Len: 9},
		gao.ByteRange{Off: 17, Len: 4},
	)

	for _, word := range byteWords {
		msg, err := bc.Decode(word, lost)
	}

A set suits any code built with the same modulus, n, k and evaluation
strategy (see next section).

# Evaluation strategies

The evaluation points dominate the cost of both operations, and NewCode picks
them:

  - The number theoretic transform evaluates at roots of unity and decodes
    with NTT-based polynomial arithmetic, in quasi-linear time. It needs n to
    be a power of two dividing p-1.
  - Pointwise evaluation at 1, 2, ..., n works over any prime field for any
    0 < n < p, but is quadratic in n.

The NTT is used whenever the field and n permit, and otherwise NewCode falls
back to pointwise evaluation silently. At large n that difference is
substantial, so if the fast path is a requirement rather than a preference,
say so with [RequireNTT] — or check [Code.UsesNTT] afterwards. [Pointwise]
forces the classical path.

To get the NTT, choose a prime with a large power of two dividing p-1, and
size it against 2n: evaluating needs an n-point transform, and decoding
needs a 2n-point one for the products inside the partial GCD, so the
strategy requires both.

If you are not sure; use field.NTTFriendlyPrime, which is 57 bits wide and
can support large n,k values.

# Input mutation and concurrency

[Code.Decode] does not modify its input: it works on a copy of the slice it is
given.

A [Code] is immutable after construction and safe for concurrent use by
multiple goroutines, as are [ErasureSet] and [ByteCode].
*/
package gao
