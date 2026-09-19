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

// MaxBytes is the largest payload [Code.EncodeBytes] accepts.
func (gao *Code) MaxBytes() int {
	return gao.K() * gao.maxPayloadPerSymbol()
}

// EncodeBytes encodes data and returns the codeword as bytes, ready to send or store.
//
// Short data is zero-padded, and the padding is indistinguishable from payload
// afterwards, so the caller has to carry the original length: [Code.DecodeBytes] returns
// MaxBytes bytes whatever was encoded.
//
// It returns ErrDataTooLarge if data exceeds MaxBytes.
func (gao *Code) EncodeBytes(data []byte) ([]byte, error) {
	if len(data) > gao.MaxBytes() {
		return nil, fmt.Errorf("%w: %d bytes exceeds %d", ErrDataTooLarge, len(data), gao.MaxBytes())
	}

	word, err := gao.Encode(bytes2symbols(data, gao.maxPayloadPerSymbol(), gao.K()))
	if err != nil {
		return nil, err
	}

	return gao.marshal(word), nil
}

// DecodeBytes decodes a codeword produced by [Code.EncodeBytes], repairing it, and
// returns MaxBytes payload bytes.
//
// `lost` states byte ranges that aren't known (erasures).
// A symbol any range touches is erased whole.
// ranges may overlap, repeat, or fall partly outside the codeword.
//
// It returns ErrMismatchedLengths if raw is not the length EncodeBytes produces.
func (gao *Code) DecodeBytes(raw []byte, lost ...ByteRange) ([]byte, error) {
	word, err := gao.unmarshal(raw)
	if err != nil {
		return nil, err
	}

	msg, err := gao.Decode(word, gao.erasedSymbols(lost)...)
	if err != nil {
		return nil, err
	}

	return symbols2bytes(msg, gao.maxPayloadPerSymbol()), nil
}

// maxPayloadPerSymbol is the most payload one symbol can carry
func (gao *Code) maxPayloadPerSymbol() int {
	return (bits.Len64(gao.PrimeField().Modulus()) - 1) / 8
}

// encodedSymbolSize is the space one symbol occupies once encoded, which is the smallest number of bytes that can hold the modulus.
func (gao *Code) encodedSymbolSize() int {
	return (bits.Len64(gao.PrimeField().Modulus()) + 7) / 8
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
func (gao *Code) marshal(c Codeword) []byte {
	w := gao.encodedSymbolSize()
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
func (gao *Code) unmarshal(b []byte) (Codeword, error) {
	w := gao.encodedSymbolSize()

	want := gao.N() * w
	if len(b) != want {
		return nil, fmt.Errorf("%w: got %d bytes, want %d", ErrMismatchedLengths, len(b), want)
	}

	f := gao.PrimeField()
	out := make(Codeword, gao.N())

	for i := range out {
		var buf [8]byte
		copy(buf[:], b[i*w:(i+1)*w])

		out[i] = f.Reduce(binary.LittleEndian.Uint64(buf[:]))
	}

	return out, nil
}

// erasedSymbols maps lost byte ranges onto the symbol indices they cover, sorted and
// without repeats. Malformed ranges are skipped rather than reported as erasures.
func (gao *Code) erasedSymbols(lost []ByteRange) []int {
	symbolRanges := gao.symbolRangesOf(lost)
	mergedRanges := mergeRanges(symbolRanges)
	return symbolRanges2Indices(mergedRanges)
}

// A symbolRange is an inclusive range of symbol indices.
type symbolRange struct {
	first, last int
}

// symbolRangesOf converts byte ranges to the symbol ranges they cover, dropping those
// that fall outside the codeword.
func (gao *Code) symbolRangesOf(lost []ByteRange) []symbolRange {
	out := make([]symbolRange, 0, len(lost))

	for _, r := range lost {
		if sr, ok := gao.byteRange2symbolRange(r); ok {
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
func (gao *Code) byteRange2symbolRange(r ByteRange) (symbolRange, bool) {
	w := gao.encodedSymbolSize()
	wire := gao.N() * w

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
