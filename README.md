# Gao Decoder for Reed-Solomon Codes in Go

[![CI](https://github.com/jonathanMweiss/go-gao/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/jonathanMweiss/go-gao/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/jonathanMweiss/go-gao/branch/main/graph/badge.svg)](https://codecov.io/gh/jonathanMweiss/go-gao)
[![Go Reference](https://pkg.go.dev/badge/github.com/jonathanmweiss/go-gao.svg)](https://pkg.go.dev/github.com/jonathanmweiss/go-gao)
[![Go](https://img.shields.io/github/go-mod/go-version/jonathanMweiss/go-gao)](go.mod)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

## Overview

Gao's decoder for Reed-Solomon codes, over prime fields up to 63 bits.

It repairs two kinds of damage. An **erasure** is a symbol you know is missing. An
**error** is a symbol that is silently wrong, at a position you do not know. The
decoder finds those positions itself.

Additionally, this repo expose the [`field`](./field) package; it is the underlying arithmetic used throughout this repo: prime fields, dense polynomial rings, NTT (including fast **polynomial Division** based on Netwon-Raphson iterations), Lagrange interpolation, and a **half-GCD** extended Euclidean algorithm.

**No dependencies.** The library imports nothing outside the Go standard library;
`testify` is used by the tests only.

**Not a cryptographic library.** The arithmetic is not constant-time, and a decoded
message is not an authenticated one.

> **Status: v0.1.** The API may still change. Pin a version.

## Where this is useful

Wherever the damage is *wrong values at positions you do not know*.

**Byzantine agreement and reliable broadcast.** Reliable broadcast ensures that as long as fewer than a third of nodes are malicious, all honest nodes deliver the same message to their underlying apps. Classically it costs $\approx n^2 \cdot \ell$ bandwidth across the system, for an $\ell$-bit message: every node echoes the whole thing to every other.

Newer protocols, in the [COOL](https://drops.dagstuhl.de/storage/00lipics/lipics-vol209-disc2021/LIPIcs.DISC.2021.17/LIPIcs.DISC.2021.17.pdf) family, offer the same guarantee for $O(\max\{n\ell,\ n^2 \log q\})$ bandwidth: linear in $n$ per bit of message rather than quadratic, plus a term that does not grow with $\ell$. They do it by encoding the message with Reed-Solomon and dispersing one symbol to each peer, so no node has to relay the whole message. That works only because the code corrects corruptions: a malicious node returns a plausible wrong symbol rather than nothing, and its position is not known in advance.


**Information dispersal.** Spread a value across `n` nodes, recover it from a sufficient subset. If a node can return corrupted data rather than simply failing to answer, the reconstruction needs error correction. Nodes known to be down are cheaper: declare them as erasures and they cost one budget unit each instead of two.

**Secret sharing and threshold protocols.** A Shamir share is a polynomial evaluated at a point over a prime field, which is the same object this library encodes. Reconstructing from shares when some may be forged is error correction over that field.

**Silent corruption.** Bit rot, a misbehaving cache, a lossy link that delivers damaged frames rather than dropping them. Anything that hands you data without telling you which part is wrong.

**Verifiable secret sharing over asynchronous networks.** See Ben-Or, Canetti and
Goldreich, ["Asynchronous secure computation"](https://dl.acm.org/doi/10.1145/167088.167109)

When the damaged positions *are* known, a disk failed or a packet never arrived, an
erasure code is the better tool (see [below](#comparison-with-erasure-coding-libraries)).

## Installation

```sh
go get github.com/jonathanmweiss/go-gao
```

Requires Go 1.25 or later.

## Usage

The library works at two levels:

- **`ByteCode`**, the higher one, when you care about bytes only.
- **`Code`**, the lower one, when you want the symbols (field elements) directly.

If you are unsure, start with `ByteCode`. It is a view of a `Code`, so you build one of
those first either way.

Both fix two lengths: `n` symbols to a codeword, `k` of them data. The difference `n-k`
is the redundancy, so a larger `n` buys budget and a larger `k` spends it on payload.

### Symbols

Decoding succeeds while

$$2e+s\le n-k$$

where $e$ is the number of errors and $s$ the number of erasures. An erasure costs half
what an error does, because the decoder does not have to spend budget locating it.
`MaxErrors()` reports $(n-k)/2$, the all-errors corner of that budget.

```go
package main

import (
	"fmt"

	"github.com/jonathanmweiss/go-gao"
	"github.com/jonathanmweiss/go-gao/field"
)

func main() {
	// can use some other prime.
	f, _ := field.NewPrimeField(field.NTTFriendlyPrime)

	const n, k = 16, 4
	code, _ := gao.NewCode(f, n, k)

	data := []uint64{10, 20, 30, 40}
	codeword, _ := code.Encode(data)

	// Corrupt 6 symbols; the maximum this specific code can repair.
	for _, i := range []int{0, 3, 5, 9, 11, 14} {
		codeword[i] = 12345
	}

	decoded, _ := code.Decode(codeword, gao.ErasureSet{})
	fmt.Println(decoded) // [10 20 30 40]
}
```

A symbol is a field element, so every value must be below the modulus.
If you are not sure what size to pick; [`field.NTTFriendlyPrime`](./field) is a default sized for both things a codeword needs:
57 bits, so seven whole bytes fit in a symbol, and $2^{32}$ divides $p-1$, so transforms
run to $2^{32}$ points.

Codewords are positional, and reordering one makes it invalid: `Encode` returns `n`
values where index `i` is the evaluation at `EvaluationPoints()[i]`, and `Decode` expects
them back in that order. Store or split them as you like, one symbol per disk, per peer
or per packet, and reassemble them in the same order before decoding. A codeword of any
other length is rejected with `ErrMismatchedLengths`, since no amount of correction fixes
framing.

### Declaring erasures

An erasure is a position you know is unusable: a missing symbol, or a range of bytes
missing from a `ByteCode` codeword.

```go
codeword, _ := code.Encode(data)

codeword[0] = 0     // value irrelevant; index 0 is declared
codeword[1] = 999   // an error: not declared

lost, err := code.Erasures(0)         // 1 of the n-k budget, vs 2 for the error
decoded, err := code.Decode(codeword, lost)
```

Whatever the slice holds at a declared index is ignored, so there is no need to blank
those entries first. Pass the zero `ErasureSet` when nothing is missing.

`Erasures` builds the locator and evaluates it: the expensive half of an erasure decode,
and it depends on the positions alone. Words that lost the same positions share one set:

```go
lost, err := code.Erasures(3, 17, 42)

for _, word := range words {
	msg, err := code.Decode(word, lost)
}
```

That is the shape of a node or a disk going down: every codeword striped across it loses
the same index. Sharing the set saves compute time.

A set suits any code built with the same modulus, `n`, `k` and strategy, so the two ends
of a link can each build their own code and still share one.

If you received only some of the symbols (a k-of-n fetch) place what you have
and name the rest:

```go
ys := make([]uint64, code.N())
seen := make([]bool, code.N())

for _, sh := range shares {
	ys[sh.Index], seen[sh.Index] = sh.Value, true
}

var erased []int
for i, ok := range seen {
	if !ok {
		erased = append(erased, i)
	}
}

lost, err := code.Erasures(erased...)
decoded, err := code.Decode(ys, lost)
```

### Bytes

`code.Bytes()` is a view of the same code that works in bytes: `Encode` returns the
codeword already packed, ready to send or store, and `Decode` takes those bytes back.

```go
code, _ := gao.NewCode(f, 16, 4)
bc := code.Bytes()

bc.MaxBytes() // 28

raw, err := bc.Encode([]byte("attack at dawn"))
// len(raw) == 128

got, err := bc.Decode(raw, gao.ErasureSet{})
// got[:14] == "attack at dawn"
```

A symbol carries whole bytes of payload and occupies whole bytes on the wire, both sized
against the modulus. Over the 57-bit `NTTFriendlyPrime` that is 7 payload bytes in an
8-byte symbol, so this code takes at most 28 bytes (`k` times 7) and produces 128 (`n`
times 8).

`Decode` returns `MaxBytes()` bytes, zero-padded past whatever was encoded. The padding
is indistinguishable from payload afterwards, so keep the original length and slice the
result.

Byte ranges known to be lost (a dropped packet, a bad sector) are named as erasures. 
A symbol any range touches is erased whole, and ranges may overlap or repeat:

```go
lost, err := bc.Erasures(gao.ByteRange{Off: 24, Len: 4}
						gao.ByteRange{Off: 23, Len: 5})
got, err := bc.Decode(raw, lost)
```

### Too many errors

Past the budget `Decode` usually returns `ErrDecoding`, but it cannot always tell: with
enough errors a received word lands closer to a *different* valid codeword, and the
decoder returns a confidently wrong message. That is inherent to Reed-Solomon codes, not
to this implementation.

### Choosing parameters

`NewCode` picks the evaluation strategy:

| Strategy | Cost | Constraint on `n` |
|---|---|---|
| NTT (roots of unity) | quasi-linear | power of two dividing `p-1` |
| Pointwise (`1..n`) | quadratic in `n` | any `0 < n < p` |

The NTT is used whenever the field and `n` permit, and otherwise `NewCode` falls back to
pointwise evaluation. The fallback is silent and the difference at large `n` is
substantial, so state it when the fast path is a requirement:

```go
code, err := gao.NewCode(f, n, k, gao.RequireNTT()) // error instead of falling back
code, err := gao.NewCode(f, n, k, gao.Pointwise())  // force the classical path
code.UsesNTT()                                      // or check afterwards
```

```go
small, _ := field.NewPrimeField(65537) // p-1 = 2^16, so transforms stop at 65536

_, err := gao.NewCode(small, 65536, 32768, gao.RequireNTT())
// ErrUnsupportedSize: decoding needs a 2n-point transform,
// and 2n=131072 does not divide p-1=65536
```

Two properties of a prime matter, and they are independent. Its **size** sets how much
payload a symbol carries: a field operation costs the same whatever the modulus, so a
small prime spends a full 64-bit multiply to move very few bits. Its **2-adicity** (the
largest power of two dividing `p-1`) bounds the transform length, and so the codeword
length. Evaluating needs an `n`-point transform and decoding needs a `2n`-point one, so
size the prime against `2n`.

If you are not sure, pick `field.NTTFriendlyPrime`.

Invalid parameters are reported at construction: `NewCode` returns `ErrUnsupportedSize`,
`ErrNSmallerThanK` or `ErrNonPositiveK` rather than failing later.

### Notes

- A `*Code` is immutable after construction and safe for concurrent use, as are an
  `ErasureSet` and a `ByteCode`.
- `Decode` never modifies its input.
- `Decode` returns a message of exactly length `k`, zero-padded when the recovered
  message has high-order zero symbols. `[]uint64{10, 20, 30, 0}` decodes back to four
  symbols, not three.
- `EvaluationPoints()` is available for interoperating with another implementation,
  but neither `Encode` nor `Decode` requires it.

See [`example_test.go`](./example_test.go) and the unit tests for further examples.

## Comparison with erasure-coding libraries

Most Go projects reach for an erasure code such as
[`klauspost/reedsolomon`](https://github.com/klauspost/reedsolomon). The two solve
different problems, and the difference is what to check before choosing either.

An erasure code takes the damaged positions as an **input**. Given a complete but
corrupted shard set it has nothing to repair: it returns the corruption unchanged, with
no error. Locating the damage is the caller's job.

This decoder treats those positions as an **output**. It finds them, at the cost of two
budget units per error against one per erasure. It also accepts erasures, so a caller
who does know some positions pays the lower price for them.

Where an erasure code is the better fit: whole shards lost at known positions, large
payloads, byte-oriented transport; it is also considerably faster.

## Benchmarks

Apple M2 Pro, plugged in rather than on battery, which is worth saying because on battery
every number below drops by about 30%. `field.NTTFriendlyPrime` with NTT evaluation.

`per call` is one operation on one core. The MB/s columns are aggregate payload
throughput with that many independent codewords in flight: encode and decode are
single-threaded, so more cores mean more operations at once, not a faster one.

```
GOMAXPROCS=1 go test -run '^$' -bench 'BenchmarkEncode$|BenchmarkErasureDecode$|BenchmarkErrorDecode$'
GOMAXPROCS=4 go test -run '^$' -bench 'Parallel'
```

### Encoding

| `n`, `k` | per call | MB/s<br>1 core | MB/s<br>2 cores | MB/s<br>4 cores |
|---|---:|---:|---:|---:|
| 256, 128 | 3.3 µs | 279 | 481 | 552 |
| 1024, 512 | 14.6 µs | 260 | 446 | 547 |
| 4096, 2048 | 60.6 µs | 245 | 437 | 619 |

Encoding is one forward transform, an order of magnitude cheaper than any decode below.
It stops scaling around four cores: a call allocates a fresh codeword and runs for only a
few microseconds, so the allocation rate reaches the collector's limit before the cores
run out. With `GOGC=off` the same benchmark keeps climbing.

### Decoding erasures — positions known

The set is built once and reused, which is how a caller with a dead node uses one.
Building it is listed separately.

| `n`, `k` | erasures | per call | MB/s<br>1 core | MB/s<br>2 cores | MB/s<br>4 cores |
|---|---:|---:|---:|---:|---:|
| 256, 128 | 32 | 38 µs | 23.5 | 42.8 | 75.1 |
| | 64 | 42 µs | 21.5 | 39.4 | 70.7 |
| | 128 *(max)* | 46 µs | 19.9 | 36.6 | 66.0 |
| 1024, 512 | 128 | 173 µs | 20.7 | 39.3 | 74.5 |
| | 256 | 187 µs | 19.1 | 36.5 | 69.3 |
| | 512 *(max)* | 202 µs | 17.6 | 33.5 | 63.8 |
| 4096, 2048 | 512 | 762 µs | 18.7 | 35.7 | 67.4 |
| | 1024 | 828 µs | 17.2 | 33.3 | 62.6 |
| | 2048 *(max)* | 894 µs | 15.8 | 30.6 | 57.6 |

### Decoding errors — positions unknown

| `n`, `k` | errors | per call | MB/s<br>1 core | MB/s<br>2 cores | MB/s<br>4 cores |
|---|---:|---:|---:|---:|---:|
| 256, 128 | 8 | 115 µs | 7.9 | 14.4 | 25.9 |
| | 32 | 221 µs | 4.1 | 7.5 | 13.7 |
| | 64 *(max)* | 263 µs | 3.5 | 6.3 | 11.2 |
| 1024, 512 | 32 | 796 µs | 4.6 | 8.6 | 15.9 |
| | 128 | 1.18 ms | 3.1 | 5.7 | 10.1 |
| | 256 *(max)* | 1.59 ms | 2.3 | 4.2 | 6.8 |
| 4096, 2048 | 128 | 3.85 ms | 3.8 | 7.1 | 12.3 |
| | 512 | 6.14 ms | 2.4 | 4.4 | 7.1 |
| | 1024 *(max)* | 8.28 ms | 1.8 | 3.2 | 4.7 |

At the same shape and the same damage count, an error costs six to ten times an erasure.
That gap is the search: an erasure decode divides by a locator the caller's indices
already determine, while an error decode runs the partial GCD to find where the damage
is. It is also why an error spends two of the `n-k` budget and an erasure one.

Decoding scales further than encoding, because each call does much more work per
allocation: erasure recovery reaches about 3.6x on four cores, error correction about
3.3x.

### Building an ErasureSet

| `n`, `k` | erasures | once |
|---|---:|---:|
| 256, 128 | 64 | 12 µs |
| 1024, 512 | 256 | 81 µs |
| 4096, 2048 | 1024 | 470 µs |

Roughly the cost of a few decodes, paid once per erasure pattern. A caller decoding a
single word pays it in full; one decoding a stripe does not.

## Explanation about the decoding logic
GAO used a strong assumption:
If there is a fix to the polynomial, it'll look like Berlekamp-Welch equation:
$E(\omega_i)*f(\omega_i) = Q(\omega_i)$ for $i\in[n]$.
$f(\omega_i)=y_i$ is the original polynomial, without corruptions.
$E(x)$ will be the ERROR LOCATOR polynomial. $E(x)$ has roots where there are corruptions,
Assume $Q(\omega_i)=E(\omega_i)*y_i$
The above equation is true for any evaluation point $\omega_i$:
$E(\omega_i)*f(\omega_i)= E(\omega_i) * y_i$, for errors equals $0$ on both sides.
for non errors we get $E(\omega_i) *y_i$, a corrupt point on both sides of the equation

Berlekamp-Welch proved that there if there is a solution, then there is only 1 unique $Q$ and $E$ in existence,
also $Q$ has a specific degree and $E$ has a specific degree.
Berlekamp-Welch then solve an equation system to find $Q$ and $E$ in O(n^3), and returns $f=Q/E$.

Gao capitalize on their proof and statement, and go extracts $Q$ and $E$ via a partial GCD, returning $g$ of the specific degree of $Q$.
$$GCD(g0,g1)=g=g0*u+g1*v$$.
That is, since $g0=(x-\omega_1)(x-\omega_2)...(x-\omega_n)$, $g$ as GCD must have roots on the evaluation points too!
meaning $g(\omega_i)=0$ for $\omega_i$ (otherwise it isn't a GCD), and $g(\omega_i)=0$ for $\omega_i$ with errors too.
Since it has the same properties of $Q$, and it has the same degree, it must be $Q$, and thus we can get $f(x)$ by dividing $g$ by $v$

If there is a remainder, then GAO's strong assumption is violated, meaning there is no solution to $Q=Ey_i$, and thus we return an error.

(if we don't have some points, we can just fill with random points as errors)

### Decode with erasures (known missing points):
an erasure filled with an arbitrary value is simply an error, and plain GAO decodes it
as-is: the error locator picks up $\omega_j$ as one of its roots like any other corruption, and everything
goes through as long as $s+e\le\frac{n-k}{2}$. But that way each erasure costs a *whole* error, because
the decoder spends budget discovering a location we already knew.

Since we do know those locations, we can hand the decoder their roots instead of making it find them.
Split the locator into the part we know and the part we don't: build an erasure locator over the missing
points only (same shape as $g_0$, with fewer roots), $S(x) = \prod_j (x-\omega_j)$ for each erased index
$j$, and let $E(x)$ cover only the genuinely unknown corruptions. With $\tilde{E}(x)=S(x)E(x)$ and
$\tilde{Q}(x)=S(x)Q(x)$, Berlekamp-Welch reads:

$$\tilde{E}(\omega_i)\cdot y_i = \tilde{E}(\omega_i)\cdot f(\omega_i) = \tilde{Q}(\omega_i) \qquad \forall i\in[n]$$

Every point is covered: at an erased $i$ we have $S(\omega_i)=0$, at a corrupted $i$ we have $E(\omega_i)=0$,
and everywhere else $y_i=f(\omega_i)$. The first case is why the filler value is free — whatever we put at
an erased point is annihilated by $S(\omega_i)=0$, so we use zeros, but any value decodes the same.

So we can solve for $\tilde{E}(x)$ and $\tilde{Q}(x)$ with GAO, given two corrections:

First, the stop degree determined by $Q$ is $e<\frac{n-k}{2}$; now we want to find $\tilde{Q}=S\cdot E \cdot f$. This polynomial has a particular degree, too; from RS theorem the degree $< \frac{s}{2}+e< \frac{n-k+1}{2}$.
So we bump the Stop degree by $s/2$.

Second, we need $g_1(x)\cdot S(x)$. This is a bit troublesome, since a true product would need more points.
Instead, we compute $\widetilde{g_1(x)}=g_1(x)S(x) \bmod g_0(x)$. This is okay, since our wanted result is
always $\bmod\ g_0(x)$. In code this is done by scaling each received value by $S(\omega_i)$ before
interpolating: the degree-$<n$ interpolant of $S(\omega_i)y_i$ *is* $g_1S \bmod g_0$.

Because $\widetilde{g_1}$ already carries $S$, the Bézout coefficient that the partial GCD returns is the
plain error locator: $\tilde{Q}\equiv E\cdot\widetilde{g_1} \pmod{g_0}$, so we get $g=\underbrace{u\cdot g_0}_{0} +v \cdot \tilde{E}$.

In the no-erasure case, partial GCD gives $g/v = f$ directly.
With erasures, partial GCD gives $g(x)/v=\frac{S\cdot E\cdot f}{E}=S(x)\cdot f(x)=\widetilde{f(x)}$,
thus the return value is one more division, $\frac{\widetilde{f(x)}}{S(x)}$.

If either division leaves a remainder, the decoding assumptions are inconsistent and we return an error.

## Contributing
Contributions are welcome! If you’d like to contribute, please open an issue or submit a pull request.

## References
- Gao Shuhong. "A new algorithm for decoding Reed-Solomon codes", in Communications, Information and Network Security,  2003.
- Reed, I. S., & Solomon, G. (1960). "Polynomial codes over certain finite fields."
- [Reed-Solomon Codes - Wikipedia](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)
- Chen, Jinyuan. "Optimal Error-Free Multi-Valued Byzantine Agreement", DISC 2021 —
  the COOL protocol, which uses Reed-Solomon error correction to bound the communication
  cost of agreement.  (for a more digestable read: https://decentralizedthoughts.github.io/2025-08-01-graded-dispersal/)

## Author
Jonathan Weiss ([@jonathanmweiss](https://github.com/jonathanmweiss))

## License
Copyright 2025-2026 Jonathan Weiss.

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the
full text and [NOTICE](NOTICE) for attribution and third-party dependencies.

This code was originally developed as part of the
[Cohort](https://github.com/jonathanmweiss/Cohort) project and extracted into a
standalone module.
