// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"bytes"
	"fmt"
	"math/bits"
	"math/rand"
	"testing"

	"github.com/jonathanmweiss/go-gao/field"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// corruptSymbol flips the low byte of symbol i in an encoded codeword. The value moves
// by at most 255, which is far below the modulus, so the symbol is always a different
// field element after reduction -- never an accidental no-op.
func corruptSymbol(code *Code, raw []byte, i int) {
	raw[i*code.encodedSymbolSize()] ^= 0xFF
}

// TestSymbolWidthsStayInTheField checks both. maxPayloadPerSymbol must be the widest
// payload that stays below the modulus -- one byte more must not fit. encodedSymbolSize
// must represent p-1, the largest a codeword symbol can be.
func TestSymbolWidthsStayInTheField(t *testing.T) {
	for _, prime := range []uint64{65537, field.NTTFriendlyPrime} {
		code, err := NewCode(newfield(t, prime), 16, 4)
		require.NoError(t, err)

		pay, wire := code.maxPayloadPerSymbol(), code.encodedSymbolSize()
		require.Positive(t, pay, "p=%d", prime)
		require.Greater(t, wire, pay, "p=%d: the wire needs a byte the payload does not", prime)

		widest := bytes.Repeat([]byte{0xFF}, pay)
		require.Less(t, bytes2symbols(widest, pay, 1)[0], prime,
			"p=%d: %d bytes of 0xFF must stay below the modulus", prime, pay)

		wider := bytes.Repeat([]byte{0xFF}, pay+1)
		require.GreaterOrEqual(t, bytes2symbols(wider, pay+1, 1)[0], prime,
			"p=%d: %d bytes would also fit, so maxPayloadPerSymbol is too small", prime, pay+1)

		require.LessOrEqual(t, bits.Len64(prime-1), 8*wire, "p=%d: p-1 must fit in %d bytes", prime, wire)
	}
}

// TestEncodeBytesRoundTrip covers every payload length from empty to full capacity,
// including the lengths that land mid-symbol.
func TestEncodeBytesRoundTrip(t *testing.T) {
	for _, prime := range []uint64{65537, field.NTTFriendlyPrime} {
		code, err := NewCode(newfield(t, prime), 32, 8)
		require.NoError(t, err)

		rng := rand.New(rand.NewSource(int64(prime)))

		for length := 0; length <= code.MaxBytes(); length++ {
			data := make([]byte, length)
			rng.Read(data)

			raw, err := code.EncodeBytes(data)
			require.NoError(t, err, "p=%d len=%d", prime, length)
			require.Len(t, raw, code.N()*code.encodedSymbolSize(), "p=%d len=%d", prime, length)

			got, err := code.DecodeBytes(raw)
			require.NoError(t, err, "p=%d len=%d", prime, length)

			require.Len(t, got, code.MaxBytes(), "p=%d len=%d", prime, length)
			require.True(t, bytes.Equal(data, got[:length]), "p=%d len=%d: payload differs", prime, length)
			require.True(t, bytes.Equal(make([]byte, code.MaxBytes()-length), got[length:]),
				"p=%d len=%d: padding must be zero", prime, length)
		}
	}
}

// TestEncodeBytesCorrectsErrors checks the byte path carries error correction through,
// at the full undeclared-error budget.
func TestEncodeBytesCorrectsErrors(t *testing.T) {
	const n, k = 64, 16

	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	rng := rand.New(rand.NewSource(11))

	data := make([]byte, code.MaxBytes())
	rng.Read(data)

	raw, err := code.EncodeBytes(data)
	require.NoError(t, err)

	for _, i := range rng.Perm(n)[:code.MaxErrors()] {
		corruptSymbol(code, raw, i)
	}

	got, err := code.DecodeBytes(raw)
	require.NoError(t, err)
	require.True(t, bytes.Equal(data, got))
}

// TestDecodeBytesWithLostRanges loses a contiguous byte run, as a dropped packet would,
// and names it rather than paying the undeclared-error price for it.
func TestDecodeBytesWithLostRanges(t *testing.T) {
	const n, k = 64, 16

	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), n, k, RequireNTT())
	require.NoError(t, err)

	rng := rand.New(rand.NewSource(8))

	data := make([]byte, code.MaxBytes())
	rng.Read(data)

	raw, err := code.EncodeBytes(data)
	require.NoError(t, err)

	lost := ByteRange{Off: 77, Len: 130}
	for i := lost.Off; i < lost.Off+lost.Len; i++ {
		raw[i] = 0xAA
	}

	got, err := code.DecodeBytes(raw, lost)
	require.NoError(t, err)
	require.True(t, bytes.Equal(data, got))
}

// TestDecodeBytesToleratesRedundantRanges: a caller reporting the same loss twice, or as
// several overlapping pieces, must not trip ErrDuplicateErasure.
func TestDecodeBytesToleratesRedundantRanges(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 64, 16, RequireNTT())
	require.NoError(t, err)

	rng := rand.New(rand.NewSource(12))

	data := make([]byte, code.MaxBytes())
	rng.Read(data)

	raw, err := code.EncodeBytes(data)
	require.NoError(t, err)

	for i := 80; i < 160; i++ {
		raw[i] = 0xAA
	}

	got, err := code.DecodeBytes(raw,
		ByteRange{Off: 80, Len: 80},
		ByteRange{Off: 80, Len: 80}, // the same loss again
		ByteRange{Off: 96, Len: 16}, // and a piece of it
	)
	require.NoError(t, err)
	require.True(t, bytes.Equal(data, got))
}

// TestErasedSymbols checks the mapping from lost byte ranges to erasure indices,
// including partial symbols, clamping, out-of-range input and overlap.
func TestErasedSymbols(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 16, 4)
	require.NoError(t, err)

	require.Equal(t, 8, code.encodedSymbolSize())

	for _, tc := range []struct {
		name string
		lost []ByteRange
		want []int
	}{
		{"one whole symbol", []ByteRange{{0, 8}}, []int{0}},
		{"one byte inside a symbol erases all of it", []ByteRange{{3, 1}}, []int{0}},
		{"straddling two symbols", []ByteRange{{7, 2}}, []int{0, 1}},
		{"a run of three", []ByteRange{{8, 24}}, []int{1, 2, 3}},
		{"clamped at the end", []ByteRange{{8*15 + 4, 99}}, []int{15}},
		{"negative offset is clipped", []ByteRange{{-4, 12}}, []int{0}},
		{"past the end", []ByteRange{{8 * 16, 4}}, nil},
		{"empty", []ByteRange{{0, 0}}, nil},
		{"no ranges", nil, nil},
		{"overlapping ranges dedupe", []ByteRange{{0, 3}, {4, 2}}, []int{0}},
		{"adjacent ranges merge", []ByteRange{{0, 8}, {8, 8}}, []int{0, 1}},
		{"unordered input comes back sorted", []ByteRange{{24, 8}, {0, 8}}, []int{0, 3}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			assert.Equal(t, tc.want, code.erasedSymbols(tc.lost))
		})
	}
}

// TestDecodeBytesReducesRatherThanRejects checks that symbols corrupted above the
// modulus are reduced and still decode, rather than being refused at the door.
func TestDecodeBytesReducesRatherThanRejects(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 64, 16, RequireNTT())
	require.NoError(t, err)

	rng := rand.New(rand.NewSource(3))

	data := make([]byte, code.MaxBytes())
	rng.Read(data)

	raw, err := code.EncodeBytes(data)
	require.NoError(t, err)

	// Smash each symbol's high wire byte, pushing it above the modulus.
	w := code.encodedSymbolSize()
	for _, i := range rng.Perm(code.N())[:code.MaxErrors()] {
		raw[i*w+w-1] = 0xFF
	}

	got, err := code.DecodeBytes(raw)
	require.NoError(t, err, "out-of-range symbols must not be rejected")
	require.True(t, bytes.Equal(data, got))
}

// TestEncodeBytesRejectsOversizedData checks the capacity boundary: MaxBytes is
// accepted, one more is not.
func TestEncodeBytesRejectsOversizedData(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 16, 4)
	require.NoError(t, err)

	require.Equal(t, 28, code.MaxBytes(), "4 symbols * 7 bytes")

	_, err = code.EncodeBytes(make([]byte, code.MaxBytes()))
	assert.NoError(t, err)

	_, err = code.EncodeBytes(make([]byte, code.MaxBytes()+1))
	assert.ErrorIs(t, err, ErrDataTooLarge)
}

// TestDecodeBytesRejectsBadLength checks a wrong length fails loudly, including a buffer
// sized at the payload width rather than the wire width.
func TestDecodeBytesRejectsBadLength(t *testing.T) {
	code, err := NewCode(newfield(t, field.NTTFriendlyPrime), 16, 4)
	require.NoError(t, err)

	_, err = code.DecodeBytes(make([]byte, 3))
	assert.ErrorIs(t, err, ErrMismatchedLengths)

	_, err = code.DecodeBytes(make([]byte, code.N()*code.maxPayloadPerSymbol()))
	assert.ErrorIs(t, err, ErrMismatchedLengths)
}

// TestPackUnpackIsLossless checks the packing alone, independent of coding.
func TestPackUnpackIsLossless(t *testing.T) {
	rng := rand.New(rand.NewSource(5))

	for _, bps := range []int{1, 2, 3, 7} {
		for _, k := range []int{1, 4, 9} {
			t.Run(fmt.Sprintf("bps=%d/k=%d", bps, k), func(t *testing.T) {
				data := make([]byte, k*bps)
				rng.Read(data)

				require.True(t, bytes.Equal(data, symbols2bytes(bytes2symbols(data, bps, k), bps)))
			})
		}
	}
}

// FuzzEncodeBytesRoundTrip drives arbitrary payloads and damage through the byte path.
// The seeds cover the packing boundaries: a length landing mid-symbol, a payload of
// exactly capacity, an empty one, and one ending in a zero byte.
func FuzzEncodeBytesRoundTrip(fz *testing.F) {
	fz.Add([]byte(nil), uint16(0), uint16(0))
	fz.Add([]byte{0}, uint16(1), uint16(0))
	fz.Add(bytes.Repeat([]byte{0xFF}, 112), uint16(0), uint16(48))
	fz.Add([]byte("payload ending in a zero\x00"), uint16(3), uint16(2))

	const n, k = 64, 16

	code, err := NewCode(newfield(fz, field.NTTFriendlyPrime), n, k, RequireNTT())
	if err != nil {
		fz.Fatal(err)
	}

	fz.Fuzz(func(t *testing.T, data []byte, rawE, rawS uint16) {
		if len(data) > code.MaxBytes() {
			data = data[:code.MaxBytes()]
		}

		raw, err := code.EncodeBytes(data)
		require.NoError(t, err)

		// Split the budget: 2e + s <= n-k, counted in symbols.
		budget := n - k
		e := int(rawE) % (budget/2 + 1)
		s := int(rawS) % max(1, budget-2*e+1)

		rng := rand.New(rand.NewSource(int64(len(data))*7919 + int64(rawE)*31 + int64(rawS)))

		perm := rng.Perm(n)
		for _, i := range perm[:e] {
			corruptSymbol(code, raw, i)
		}

		w := code.encodedSymbolSize()

		lost := make([]ByteRange, 0, s)
		for _, i := range perm[e : e+s] {
			for j := range w {
				raw[i*w+j] = 0xAA
			}

			lost = append(lost, ByteRange{Off: i * w, Len: w})
		}

		got, err := code.DecodeBytes(raw, lost...)
		require.NoError(t, err, "len=%d e=%d s=%d", len(data), e, s)

		require.Len(t, got, code.MaxBytes())
		require.True(t, bytes.Equal(data, got[:len(data)]), "len=%d e=%d s=%d: payload differs", len(data), e, s)
	})
}

// TestEncodedSizeIsOneMoreThanPayload pins the identity encodedSymbolSize relies on:
// ceil(bitlen(p)/8) == (bitlen(p)-1)/8 + 1, for every prime the field accepts.
func TestEncodedSizeIsOneMoreThanPayload(t *testing.T) {
	for _, prime := range []uint64{257, 65537, 4294967311, 144115158011084801, 9223372006790004737} {
		code, err := NewCode(newfield(t, prime), 16, 4)
		require.NoError(t, err)

		require.Equal(t, (bits.Len64(prime)+7)/8, code.encodedSymbolSize(),
			"p=%d: must match the direct ceil(bitlen/8)", prime)
		require.Equal(t, code.maxPayloadPerSymbol()+1, code.encodedSymbolSize(), "p=%d", prime)
	}
}
