# Gao Decoder for Reed-Solomon Codes in Go

## Overview
This repository implements Gao's decoder for Reed-Solomon error-correcting codes in Go. 
The decoder can perform robust interpolation (fix corruptions), given a list of the original evaluation points.


## Installation
To use the Gao decoder in your Go project, install it using:

```sh
 go get github.com/jonathanmweiss/go-gao
```

Then, import it in your code:

```go
 import "github.com/jonathanmweiss/go-gao"
```

## Usage
See unit tests for example.


## Planned Improvements:

- Implementing positive-wrapped and negative-wrapped NTT forward evaluation for fast encoding.
- Roots of unity fast decoding (Reducing the algorithms runtime from O(n^2 * logn) to O(n^2)).
- Optimised decoding for erasure only (erasure only faults decrease redundancies).
- Code optimisations to reduce allocations (and GC).
- Remove the Lattigo import by implementing a prime factorization algorithm (and switch to the MIT license).

## Contributing
Contributions are welcome! If you’d like to contribute, please open an issue or submit a pull request.

## References
- Gao Shuhong. "A new algorithm for decoding Reed-Solomon codes", in Communications, Information and Network Security,  2003.
- Reed, I. S., & Solomon, G. (1960). "Polynomial codes over certain finite fields."
- [Reed-Solomon Codes - Wikipedia](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)


