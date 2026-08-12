# Gao Decoder for Reed-Solomon Codes in Go

## Overview
This repository implements Gao's decoder for Reed-Solomon error-correcting codes in Go. 
The decoder can perform robust interpolation (fix corruptions), given a list of the original evaluation points.


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
	params, _ := gao.NewCodeParameters(gao.NewNttEvaluator(f), n, k)
	code := gao.NewCodeGao(params)

	data := []uint64{10, 20, 30, 40}
	codeword, _ := code.EncodeToSlice(data)

	// Corrupt 6 symbols -- the maximum this code can repair.
	for _, i := range []int{0, 3, 5, 9, 11, 14} {
		codeword[i] = 12345
	}

	decoded, _ := code.DecodeFromSlice(codeword)
	fmt.Println(decoded) // [10 20 30 40]
}
```

Two evaluation-point strategies are available. `NewNttEvaluator` uses roots of
unity and NTT-based polynomial arithmetic (fast; needs a field with suitable
roots of unity). `NewSlowEvaluator` uses successive powers of a generator and
classical arithmetic, and works for any prime field.

`Encode`/`Decode` operate on `map[uint64]uint64` of evaluation point to value,
which lets you express *erasures* by simply omitting entries.
`EncodeToSlice`/`DecodeFromSlice` are the positional equivalents.

See the unit tests for further examples.

## Planned Improvements:

- Remove the Lattigo import by implementing a prime factorization algorithm.

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



### Decode with erasures (known missing points):
Adding erasures uses similar ideas, in recursion!
When we have points that we know are missing, we create a specialized erasure locator (similar to g_0, but just for erasures): $S(x) = \prod (x-\omega_j)$ for each erased $j$  (same shape as g_0, with less points).
Now, when we multiply $S(x)$ with the berlekamp-welch equation from both sides we get:

$$\underbrace{ S(\omega_i)E(\omega_i)}_{\tilde{E}(x)}\cdot f(\omega_i) = \underbrace{S(\omega_i)Q(\omega_i)}_{\tilde{Q}(x)} $$

That is, we get $\tilde{E}(X)$ and $\tilde{Q}(X)$ and we can solve for them using GAO, we just need to do some corrections:
First, bump the StopDegree by $deg(S(x))$, which is of-course the number of erasures.

Second, we need to get $g_1(x)\cdot S(x)$. This is a bit troublesome, since it will mean we need more points. Instead, we compute $\widetilde{g_1(x)}=g_1(x)S(x) \mod g_0(x)$. This is okay, since our wanted result is always $\bmod g_0(x)$.

In the no-erasure case, partial GCD gives $g/v = f$ directly.
With erasures, partial GCD gives $g(x)/v=S(x)\cdot f(x)=\widetilde{f(x)}$, thus the return value is the division $\frac{\widetilde{f(x)}}{S(x)}$.

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
