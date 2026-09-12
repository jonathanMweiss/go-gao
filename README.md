# Gao Decoder for Reed-Solomon Codes in Go

## Overview

This repository implements Gao's decoder for Reed-Solomon codes in Go, over
**arbitrary prime fields** up to 63 bits. Can be extended to large prime fields using CRT/RNS.

Two things distinguish it from the erasure-coding libraries most Go projects reach
for:

- **It corrects errors, not just erasures.** An erasure is a symbol you know is
  missing; an error is a symbol that is silently wrong, at a position you do not
  know. Libraries such as `klauspost/reedsolomon` repair only the former — their
  README is explicit that "the encoder does not know which parts are invalid". This
  decoder repairs both, and a mixture of the two, as long as
  `2*errors + erasures <= n-k`.
- It works over prime fields.

The polynomial and finite-field arithmetic underneath are exported as a reusable
[`field`](./field) package: prime fields, dense polynomial rings, NTT, Lagrange
interpolation, and a half-GCD extended Euclidean algorithm.

> **Status: v0.x.** The API may still change. Pin a version.

## Installation

```sh
go get github.com/jonathanmweiss/go-gao
```

Requires Go 1.25 or later.

## Usage

Encode `k` data symbols into an `n`-symbol codeword, then recover the data after
up to `(n-k)/2` symbols have been corrupted:

```go
package main

import (
	"fmt"

	"github.com/jonathanmweiss/go-gao"
	"github.com/jonathanmweiss/go-gao/field"
)

func main() {
	f, _ := field.NewPrimeField(65537)

	const n, k = 16, 4
	code, _ := gao.NewCode(f, n, k)

	data := []uint64{10, 20, 30, 40}
	codeword, _ := code.Encode(data)

	// Corrupt 6 symbols; the maximum this 
	// specific code can repair.
	for _, i := range []int{0, 3, 5, 9, 11, 14} {
		codeword[i] = 12345
	}

	decoded, _ := code.Decode(codeword)
	fmt.Println(decoded) // [10 20 30 40]
}
```

### Declaring erasures

```go
codeword, _ := code.Encode(data)

codeword[0] = 0                           // value irrelevant; index 0 is declared
codeword[1] = 999                         // an error: not declared

decoded, err := code.Decode(codeword, 0)  // costs 1 of the n-k budget, vs 2 for the error
```

Whatever the slice holds at a declared index is ignored, so there is no need to blank
those entries first.

If you received only some of the symbols — a k-of-n fetch, say — place what you have
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

decoded, err := code.Decode(ys, erased...)
```

Past the budget `Decode` usually returns `ErrDecoding`, but it cannot always tell:
with enough errors a received word can land closer to a *different* valid codeword,
and you get a confidently wrong message. That is inherent to Reed-Solomon codes, not to this implementation.

### Choosing parameters

`NewCode` picks the evaluation strategy for you:

| Strategy | Cost | Constraint on `n` |
|---|---|---|
| NTT (roots of unity) | quasi-linear | power of two dividing `p-1` |
| Pointwise (`1..n`) | quadratic in `n` | any `0 < n < p` |

The NTT is used whenever the field and `n` permit, and otherwise it falls back to
pointwise evaluation. **That fallback is silent**, and at large `n` the difference is
substantial, so if the fast path is a requirement rather than a preference, say so:

```go
code, err := gao.NewCode(f, n, k, gao.RequireNTT()) // if not enough roots of unity: error instead of falling back
code, err := gao.NewCode(f, n, k, gao.Pointwise())  // force the classical path
code.UsesNTT()                                      // or just check afterwards
```

To get the NTT, pick a prime with a large power of two dividing `p-1` — and size it
against `2n`. Evaluating a codeword needs an `n`-point transform, and decoding
multiplies polynomials of degree up to `n` inside the partial GCD and needs a `2n`-point
one, so the NTT strategy requires both.

Since the fallback is silent, say so when you need the fast path:

```go
_, err := gao.NewCode(f, 65536, 32768, gao.RequireNTT())
// ErrUnsupportedSize: ... n and 2n must both be powers of two dividing p-1 ...
```

Invalid parameters are reported at construction — `NewCode` returns
`ErrUnsupportedSize`, `ErrNSmallerThanK` or `ErrNonPositiveK` rather than failing
later.

### Notes

- A `*Code` is immutable after construction and safe for concurrent use.
- `Decode` never modifies its input.
- `Decode` returns a message of exactly length `k`, zero-padded when the recovered
  message has high-order zero symbols. `[]uint64{10, 20, 30, 0}` decodes back to four
  symbols, not three.
- Codewords are positional throughout. `EvaluationPoints()` is available for
  interoperating with another implementation, but neither `Encode` nor `Decode`
  requires it.

See [`example_test.go`](./example_test.go) and the unit tests for further examples.

## Planned Improvements

- Remove the Lattigo dependency by implementing primitive-root search directly.
  It is currently pulled in for a single call.
- Benchmark numbers in this README, rather than only in `go test -bench`.

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

First, the stop degree determined by $Q$ is $e<\frac{n-k}{2}$; now we want to find $\tilde{Q}=S\cdot E \cdot f$. This polynomial has a particular degree, too; from RS theorem the degree $< \frac{s}{2}+e< \frac{n-k+1}/2$.
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


## Author
Jonathan Weiss ([@jonathanmweiss](https://github.com/jonathanmweiss))

## License
Copyright 2025-2026 Jonathan Weiss.

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the
full text and [NOTICE](NOTICE) for attribution and third-party dependencies.

This code was originally developed as part of the
[Cohort](https://github.com/jonathanmweiss/Cohort) project and extracted into a
standalone module.
