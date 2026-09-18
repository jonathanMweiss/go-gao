package field

import "math/bits"

// This file holds the overlap-add convolution used for lopsided products.
//
// A convolution by transform costs three transforms of the whole product length, and
// that length is set by the operands together: 163 coefficients times 16384 needs 16546
// of them, which rounds up to 32768 evaluation points. The short operand fills 163 of
// those and the transform spends the rest on padding.
//
// Splitting the long operand into blocks of blockLen coefficients, written d here,
//
//	B = B_0 + x^d*B_1 + ... + x^{d*(t-1)}*B_{t-1}
//	A*B = sum_j x^(d*j) * (A * B_j)
//
// turns that into one product per block, each needing only enough points for
// d + len(A) - 1 coefficients. Consecutive terms overlap -- A*B_j spans d + len(A) - 1
// coefficients while the next block starts only d later -- so the terms are summed into
// the output rather than laid end to end.
//
// A is the same in every term, so its transform is computed once and reused, and a block
// costs one forward and one inverse transform instead of three.

// blockedConvPlan returns the block size an overlap-add convolution of operands of
// lengths ls and ll would use, and reports whether it costs less than transforming the
// whole product. Callers pass ls <= ll.
func blockedConvPlan(ls, ll int) (blockLen int, ok bool) {
	if ls <= 0 || ll <= 0 {
		return 0, false
	}

	// A transform of size n is O(n*log2(n)).
	cost := func(n int) int { return n * bits.Len(uint(n-1)) }

	// blockLen >= ls, so a block's product fits in 2*blockLen coefficients.
	blockLen = nextPow2(ls)
	blocks := (ll + blockLen - 1) / blockLen

	// One forward and one inverse per block, plus the short operand's single forward.
	blockedCost := (2*blocks + 1) * cost(2*blockLen)
	// the cost is intt(ntt(a)*ntt(b)), so 3 ops.
	standardCost := 3 * cost(nextPow2(ls+ll-1))

	// Both counts are NTT costs only: they ignore the addition loop and the fixed cost
	// per transform, which flatters the blocked side.
	// Measurement puts blocked as better when
	// standardCost * 7/8 > blockedCost, so we require that much of a margin.
	//  TestBlockedConvPlan pins it.
	return blockLen, blockedCost < standardCost-standardCost/8
}

// mulBlockedInto writes short*long into c as an overlap-add over blocks of `long`, with
// blockLen from [blockedConvPlan]. Both operands must be in the coefficient domain, and
// c may alias either of them.
func (r *PolyRing) mulBlockedInto(c, short, long *Polynomial, blockLen int) {
	f := r.f
	ls, ll := len(short.inner), len(long.inner)
	total := ls + ll - 1
	n := 2 * blockLen

	// The short operand is a factor of every block's product, so it is transformed once
	// for the whole call. Its tail has to read as zero, so this one is borrowed zeroed.
	sh := r.borrowPolyZeroed(n)
	defer r.returnPoly(sh)

	copy(sh.inner, short.inner)

	if err := r.NttForward(sh); err != nil {
		panic(err)
	}

	// Each of these evaluations multiplies one point of every block, so its Shoup factor
	// is amortized over the blocks and pays for its division after the first one.
	shoup := r.borrowPoly(n) // the shoup factors of the short operand's transform.
	defer r.returnPoly(shoup)

	// sh, shoup and chunk are all n coefficients long and none of them is resliced below,
	// so cutting them to one length lets the loops here index all three without a bounds
	// check.
	shInner := sh.inner[:n]
	shoupInner := shoup.inner[:n]

	for i, w := range shInner {
		shoupInner[i] = f.shoupFactor(w)
	}

	// The blocks are read from long and the result is summed into out, so writing
	// straight into c is only safe where c is neither operand.
	var out []uint64
	if c != short && c != long {
		out = resizeZeroed(c.inner, total)
	} else {
		out = make([]uint64, total)
	}

	chunk := r.borrowPoly(n)
	defer r.returnPoly(chunk)

	chunkInner := chunk.inner[:n]

	for off := 0; off < ll; off += blockLen {
		m := min(blockLen, ll-off)

		// The copy covers at most the lower half, and the rest still holds the previous
		// block's inverse transform, which has to read as zero here.
		copy(chunk.inner[:m], long.inner[off:off+m])
		clear(chunk.inner[m:])
		chunk.isNTT = false

		if err := r.NttForward(chunk); err != nil {
			panic(err)
		}

		// not using multpointwise because the Shoup factors are already
		// computed and amortized over the blocks.
		for i, x := range chunkInner {
			chunkInner[i] = f.mulShoup(shInner[i], shoupInner[i], x)
		}

		if err := r.nttBackwardNoTrim(chunk); err != nil {
			panic(err)
		}

		// since we multiplied the two polynomials their degree got bigger, so we can't just copy
		// the result into the output; instead, we need to add the entire polynomial to the
		// output (with the original shift).
		//
		// This term spans m + ls - 1 coefficients from off, and the next block starts
		// only blockLen later, so its tail is added onto ground the next term also covers.
		seg := out[off : off+m+ls-1]
		tmp := chunk.inner[:len(seg)]
		for i, v := range tmp {
			seg[i] = f.Add(seg[i], v)
		}
	}

	c.inner, c.f, c.isNTT = out, f, false
}
