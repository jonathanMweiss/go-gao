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

- Optimised decoding for erasure only (erasure only faults decrease redundancies).
- Remove the Lattigo import by implementing a prime factorization algorithm (and switch to the MIT license).

## Disclosure
This project used ChatGPT to implement/modify/improve fast 
polynomial algorithms via Number Theoretic Transform (NTT). 
NTT-based methods were either optimized from naive versions (e.g., turning
a recursive NTT into an iterative one with cached roots of unity)
or generated directly and tested against classical (and written by me) 
existing implementations (e.g., checking LongDivNTT results against my own implementation of LongDiv).

## Contributing
Contributions are welcome! If you’d like to contribute, please open an issue or submit a pull request.

## References
- Gao Shuhong. "A new algorithm for decoding Reed-Solomon codes", in Communications, Information and Network Security,  2003.
- Reed, I. S., & Solomon, G. (1960). "Polynomial codes over certain finite fields."
- [Reed-Solomon Codes - Wikipedia](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)


