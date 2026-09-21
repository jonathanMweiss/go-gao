// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

/*
Package gao implements Reed-Solomon error correction using Gao's decoder.

Erasure coding repairs symbols only when you already know which ones are
missing. This package also repairs symbols that are silently wrong at
positions you do not know: the decoder finds the error locations itself.

Symbols are elements of a prime field up to 63 bits wide, chosen by the
caller. Byte-oriented Reed-Solomon fixes that alphabet at GF(2^8).

# Codes

A code is built over a prime field and fixes two lengths: n symbols to a
codeword, k of them data. Encoding treats the k data symbols as the
coefficients of a polynomial and evaluates it at n points; decoding recovers
those coefficients from the n possibly corrupted values.

	f, err := field.NewPrimeField(65537)
	if err != nil {
		return err
	}

	code, err := gao.NewCode(f, 16, 4)
	if err != nil {
		return err
	}

	codeword, err := code.Encode([]uint64{10, 20, 30, 40})

[NewCode] reports bad parameters immediately, as [ErrNonPositiveK],
[ErrNSmallerThanK] or [ErrUnsupportedSize], rather than failing later.

# The decoding budget

A corrupted symbol at an unknown position is an error; one at a position the
caller knows is an erasure. Erasures are half the price, because the decoder
does not have to spend budget locating them. Decoding succeeds while

	2*errors + erasures <= n-k

[Code.MaxErrors] reports (n-k)/2, the all-errors corner of that budget.

Past the budget [Code.Decode] usually returns [ErrDecoding], but it cannot
always tell: with enough errors a received word lands closer to a different
valid codeword, and the decoder returns a confidently wrong message. That is
inherent to the code, not to this implementation.

# Codewords

Codewords are positional. [Code.Encode] returns n values, where index i is the
evaluation at [Code.EvaluationPoints] index i, and [Code.Decode] expects them
back in that order.

Erasures are named by index through [Code.Erasures]: a set built from 3 and 7
declares those positions unusable, whatever ys happens to hold at them, so
there is no need to blank them first. Pass the zero [ErasureSet] when nothing
is missing.

Building the set is the expensive half of an erasure decode and depends on the
positions alone. Words that lost the same positions share one:

	lost, err := code.Erasures(3, 7)
	for _, word := range words {
		msg, err := code.Decode(word, lost)
	}

A set suits any code built with the same modulus, n, k and evaluation strategy.

Decode returns a message of exactly length k, zero-padded when the recovered
message has high-order zero symbols.

# Bytes

[Code.Bytes] is a view of the same code that works in bytes rather than
symbols. [ByteCode.Encode] packs the payload into symbols and serialises the
codeword, and [ByteCode.Decode] takes those bytes back:

	bc := code.Bytes()

	raw, err := bc.Encode([]byte("attack"))
	got, err := bc.Decode(raw, gao.ErasureSet{})

A symbol carries whole bytes of payload and occupies whole bytes on the wire,
both sized against the modulus. Over p=65537 that is 2 payload bytes in a
3-byte symbol, so the code above takes at most 8 bytes ([ByteCode.MaxBytes],
k times 2) and produces 48 (n times 3).

Decode returns MaxBytes bytes whatever was encoded, zero-padded past the
payload, and the padding is indistinguishable from payload afterwards: carry
the original length and slice the result. A codeword that is not the length
Encode produces is rejected with [ErrMismatchedLengths].

Lost byte ranges are named through [ByteCode.Erasures], which builds the same
[ErasureSet] the symbol interface takes, so a batch sharing a loss pattern
reuses one set either way. A symbol any range touches is erased whole, and
ranges may overlap, repeat, or fall partly outside the codeword:

	lost, err := bc.Erasures(gao.ByteRange{Off: 6, Len: 9})
	got, err := bc.Decode(raw, lost)

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

# Input mutation and concurrency

[Code.Decode] does not modify its input: it works on a copy of the slice it is
given.

A [Code] is immutable after construction and safe for concurrent use by
multiple goroutines, as are [ErasureSet] and [ByteCode].
*/
package gao
