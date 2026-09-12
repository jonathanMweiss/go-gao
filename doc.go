// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

/*
Package gao implements Reed-Solomon error correction using Gao's decoder.

Erasure coding repairs symbols only when you already know which ones are
missing. This package also repairs symbols that are silently wrong at
positions you do not know — the decoder finds the error locations itself. It
operates over arbitrary prime fields up to 63 bits rather than a fixed
GF(2^8), so the symbol alphabet is the caller's choice.

# Codes

A code is a field, a codeword length n, and a message length k. Encoding
treats the k data symbols as the coefficients of a polynomial and evaluates it
at n points; decoding recovers those coefficients from the n possibly
corrupted values.

	f, err := field.NewPrimeField(65537)
	if err != nil {
		return err
	}

	code, err := gao.NewCode(f, 16, 4)
	if err != nil {
		return err
	}

	codeword, err := code.EncodeToSlice([]uint64{10, 20, 30, 40})

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

Erasures are named by index: Decode(ys, 3, 7) declares positions 3 and 7
unusable, whatever ys happens to hold at them, so there is no need to blank
them first.

Decode returns a message of exactly length k, zero-padded when the recovered
message has high-order zero symbols.

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

To get the NTT, choose a prime with a large power of two dividing p-1. 65537
admits any n up to 2^16; 929, the PDF417 field, has p-1 = 2^5 * 29 and so
reaches only n = 32.

# Input mutation and concurrency

[Code.Decode] does not modify its input: it works on a copy of the slice it is
given.

A [Code] is immutable after construction and safe for concurrent use by
multiple goroutines.
*/
package gao
