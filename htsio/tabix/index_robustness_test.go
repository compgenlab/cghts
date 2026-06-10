package tabix

import (
	"bytes"
	"encoding/binary"
	"testing"
)

// le32 encodes int32 values little-endian, the wire format of TBI/CSI counts.
func le32(vals ...int32) []byte {
	buf := &bytes.Buffer{}
	for _, v := range vals {
		_ = binary.Write(buf, binary.LittleEndian, v)
	}
	return buf.Bytes()
}

// TestReadBinIndexMalformedNoOOM feeds index bytes with implausible or negative
// reference counts and asserts the parser returns an error promptly rather than
// pre-allocating from an attacker-controlled count or panicking.
func TestReadBinIndexMalformedNoOOM(t *testing.T) {
	cases := map[string][]byte{
		"huge n_ref, no data": le32(0x7FFFFFFF),
		"negative n_ref":      le32(-1),
	}
	for name, data := range cases {
		t.Run(name, func(t *testing.T) {
			if _, err := readBinIndex(bytes.NewReader(data)); err == nil {
				t.Fatalf("expected error for %q, got nil", name)
			}
		})
	}
}

// TestReadRefBinsMalformedNoOOM exercises the per-reference bin/chunk/interval
// counts, all of which were previously used to size slices directly.
func TestReadRefBinsMalformedNoOOM(t *testing.T) {
	cases := map[string][]byte{
		"huge n_bins":          le32(0x7FFFFFFF),
		"negative n_bins":      le32(-1),
		"huge n_chunks":        append(le32(1, 0), le32(0x7FFFFFFF)...), // nBins=1, binNum=0, nChunks=huge
		"negative n_chunks":    append(le32(1, 0), le32(-1)...),
		"huge n_intervals":     append(le32(0), le32(0x7FFFFFFF)...), // nBins=0, nIntervals=huge
		"negative n_intervals": append(le32(0), le32(-1)...),
	}
	for name, data := range cases {
		t.Run(name, func(t *testing.T) {
			var ref binIndex
			if err := readRefBins(bytes.NewReader(data), &ref); err == nil {
				t.Fatalf("expected error for %q, got nil", name)
			}
		})
	}
}

// TestReadRefBinsValidEmpty confirms the hardening does not reject a legitimate
// empty reference (no bins, no intervals).
func TestReadRefBinsValidEmpty(t *testing.T) {
	var ref binIndex
	if err := readRefBins(bytes.NewReader(le32(0, 0)), &ref); err != nil {
		t.Fatalf("valid empty ref rejected: %v", err)
	}
}

// TestReadCSIRefMalformedNoOOM covers the CSI variant of the per-reference path.
func TestReadCSIRefMalformedNoOOM(t *testing.T) {
	idx := &CSIIndex{}
	var ref csiRefIndex
	if err := idx.readCSIRef(bytes.NewReader(le32(0x7FFFFFFF)), &ref); err == nil {
		t.Fatal("expected error for huge CSI n_bins, got nil")
	}
}

// TestParseTabixAuxMalformedNoPanic feeds the in-memory CSI aux block with a
// huge/short names length and asserts it never panics on the slice bounds.
func TestParseTabixAuxMalformedNoPanic(t *testing.T) {
	hugeNames := make([]byte, 28)
	binary.LittleEndian.PutUint32(hugeNames[24:28], 0xFFFFFFFF)

	cases := map[string][]byte{
		"huge namesLen": hugeNames,
		"too short":     make([]byte, 10),
		"exactly 28":    make([]byte, 28),
	}
	for name, data := range cases {
		t.Run(name, func(t *testing.T) {
			idx := &CSIIndex{}
			_ = idx.parseTabixAux(data) // must return without panicking
		})
	}
}
