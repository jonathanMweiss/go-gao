// Copyright 2025-2026 Jonathan Weiss
// SPDX-License-Identifier: Apache-2.0

package gao

import (
	"encoding/binary"
	"fmt"
	"math/bits"
	"slices"
)

// A ByteRange marks bytes [Off, Off+Len) of an encoded codeword as lost.
type ByteRange struct {
	Off, Len int
}

// A ByteCode encodes and decodes byte payloads over the code it wraps.
//
// It carries no state beyond that code, so it is safe for concurrent use and two
// ByteCodes over the same Code are interchangeable.
type ByteCode struct {
	code *Code
}

// Bytes returns a view of the code that speaks bytes instead of symbols.
func (gao *Code) Bytes() *ByteCode {
	return &ByteCode{code: gao}
}

// Code returns the code the byte view was built from.
func (bc *ByteCode) Code() *Code {
	return bc.code
}

// MaxBytes is the largest payload [ByteCode.Encode] accepts.
func (bc *ByteCode) MaxBytes() int {
	return bc.code.K() * bc.maxPayloadPerSymbol()
}

// Encode encodes data and returns the codeword as bytes, ready to send or store.
//
// Short data is zero-padded, and the padding is indistinguishable from payload
// afterwards, so the caller has to carry the original length: [ByteCode.Decode] returns
// MaxBytes bytes whatever was encoded.
//
// It returns ErrDataTooLarge if data exceeds MaxBytes.
func (bc *ByteCode) Encode(data []byte) ([]byte, error) {
	if len(data) > bc.MaxBytes() {
		return nil, fmt.Errorf("%w: %d bytes exceeds %d", ErrDataTooLarge, len(data), bc.MaxBytes())
	}

	word, err := bc.code.Encode(bytes2symbols(data, bc.maxPayloadPerSymbol(), bc.code.K()))
	if err != nil {
		return nil, err
	}

	return bc.marshal(word), nil
}

// Decode decodes a codeword produced by [ByteCode.Encode], repairing it, and returns
// MaxBytes payload bytes.
//
// erasures names the byte ranges that aren't known, and [ByteCode.Erasures] builds it.
// Pass the zero ErasureSet when nothing is missing.
//
// It returns ErrMismatchedLengths if raw is not the length Encode produces.
func (bc *ByteCode) Decode(raw []byte, erasures ErasureSet) ([]byte, error) {
	word, err := bc.unmarshal(raw)
	if err != nil {
		return nil, err
	}

	msg, err := bc.code.Decode(word, erasures)
	if err != nil {
		return nil, err
	}

	return symbols2bytes(msg, bc.maxPayloadPerSymbol()), nil
}

// maxPayloadPerSymbol is the most payload one symbol can carry
func (bc *ByteCode) maxPayloadPerSymbol() int {
	return (bits.Len64(bc.code.PrimeField().Modulus()) - 1) / 8
}

// encodedSymbolSize is the space one symbol occupies once encoded, which is the smallest number of bytes that can hold the modulus.
func (bc *ByteCode) encodedSymbolSize() int {
	return (bits.Len64(bc.code.PrimeField().Modulus()) + 7) / 8
}

// bytes2symbols splits data into k symbols of payloadBytes bytes each, little-endian. The
// last partial group and any symbols past the data are zero.
func bytes2symbols(data []byte, payloadBytes, k int) []uint64 {
	out := make([]uint64, k)

	for i := range out {
		start := i * payloadBytes
		if start >= len(data) {
			break
		}

		var buf [8]byte
		copy(buf[:], data[start:min(start+payloadBytes, len(data))])

		out[i] = binary.LittleEndian.Uint64(buf[:])
	}

	return out
}

// symbols2bytes is bytes2symbols reversed: each symbol contributes its low payloadBytes
// bytes, little-endian.
func symbols2bytes(msg []uint64, payloadBytes int) []byte {
	out := make([]byte, len(msg)*payloadBytes)

	for i, sym := range msg {
		var buf [8]byte
		binary.LittleEndian.PutUint64(buf[:], sym)

		copy(out[i*payloadBytes:], buf[:payloadBytes])
	}

	return out
}

// marshal serialises a codeword little-endian at encodedSymbolSize bytes per symbol.
func (bc *ByteCode) marshal(c Codeword) []byte {
	w := bc.encodedSymbolSize()
	out := make([]byte, len(c)*w)

	for i, sym := range c {
		var buf [8]byte
		binary.LittleEndian.PutUint64(buf[:], sym)

		copy(out[i*w:], buf[:w])
	}

	return out
}

// unmarshal is marshal reversed.
//
// Symbols at or above the modulus are reduced to the modulus; items over the modulus
// are corrupted (not marked as erasures).
// A wrong length is rejected, since no amount of correction
// fixes framing.
func (bc *ByteCode) unmarshal(b []byte) (Codeword, error) {
	w := bc.encodedSymbolSize()

	want := bc.code.N() * w
	if len(b) != want {
		return nil, fmt.Errorf("%w: got %d bytes, want %d", ErrMismatchedLengths, len(b), want)
	}

	f := bc.code.PrimeField()
	out := make(Codeword, bc.code.N())

	for i := range out {
		var buf [8]byte
		copy(buf[:], b[i*w:(i+1)*w])

		out[i] = f.Reduce(binary.LittleEndian.Uint64(buf[:]))
	}

	return out, nil
}

// erasedSymbols maps lost byte ranges onto the symbol indices they cover, sorted and
// without repeats. Malformed ranges are skipped rather than reported as erasures.
func (bc *ByteCode) erasedSymbols(lost []ByteRange) []int {
	symbolRanges := bc.symbolRangesOf(lost)
	mergedRanges := mergeRanges(symbolRanges)
	return symbolRanges2Indices(mergedRanges)
}

// A symbolRange is an inclusive range of symbol indices.
type symbolRange struct {
	first, last int
}

// symbolRangesOf converts byte ranges to the symbol ranges they cover, dropping those
// that fall outside the codeword.
func (bc *ByteCode) symbolRangesOf(lost []ByteRange) []symbolRange {
	out := make([]symbolRange, 0, len(lost))

	for _, r := range lost {
		if sr, ok := bc.byteRange2symbolRange(r); ok {
			out = append(out, sr)
		}
	}

	return out
}

// mergeRanges folds overlapping and adjacent ranges together.
//
// Adjacent ones merge too:
// [0,1] and [2,3] expand to the same indices as [0,3].
func mergeRanges(rs []symbolRange) []symbolRange {
	if len(rs) == 0 {
		return nil
	}

	slices.SortFunc(rs, func(a, b symbolRange) int { return a.first - b.first })

	merged := rs[:1]

	for _, r := range rs[1:] {
		tail := len(merged) - 1

		if r.first <= merged[tail].last+1 {
			// take the larger last so that [0,2] and [2,3] merge to [0,3]
			merged[tail].last = max(merged[tail].last, r.last)

			continue
		}

		merged = append(merged, r)
	}

	return merged
}

// symbolRanges2Indices lists every index the ranges cover,
// sized up front so it allocates once.
func symbolRanges2Indices(rs []symbolRange) []int {
	total := 0
	for _, r := range rs {
		total += r.last - r.first + 1
	}

	if total == 0 {
		return nil
	}

	out := make([]int, 0, total)

	for _, r := range rs {
		for i := r.first; i <= r.last; i++ {
			out = append(out, i)
		}
	}

	return out
}

// byteRange2symbolRange returns the symbols a lost byte range covers.
// Ranges reaching outside the codeword are clipped, and ones falling wholly outside report ok false.
func (bc *ByteCode) byteRange2symbolRange(r ByteRange) (symbolRange, bool) {
	w := bc.encodedSymbolSize()
	wire := bc.code.N() * w

	off, length := r.Off, r.Len
	if off < 0 {
		length += off
		off = 0
	}

	if length <= 0 || off >= wire {
		return symbolRange{}, false
	}

	return symbolRange{first: off / w, last: min(off+length-1, wire-1) / w}, true
}
