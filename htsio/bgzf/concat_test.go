package bgzf

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// concatStreams compresses each part as a complete BGZF stream -- data blocks
// and its own EOF marker -- and joins them, the way a separately compressed
// header and body are often concatenated. Every join leaves an empty block in
// mid-file. It returns the file and the uncompressed whole.
func concatStreams(t *testing.T, parts ...string) ([]byte, string) {
	t.Helper()
	var file bytes.Buffer
	for _, p := range parts {
		var b bytes.Buffer
		w := NewWriter(&b)
		if _, err := w.Write([]byte(p)); err != nil {
			t.Fatal(err)
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
		file.Write(b.Bytes())
	}
	return file.Bytes(), strings.Join(parts, "")
}

// Parts big enough that the body spans several blocks.
func concatParts() []string {
	header := "#header line one\n#header line two\n"
	var body strings.Builder
	for i := 0; body.Len() < 200_000; i++ {
		body.WriteString("chr1\t")
		body.WriteString(strings.Repeat("x", i%50))
		body.WriteString("\n")
	}
	return []string{header, body.String()}
}

func TestReaderReadsThroughMidFileEOF(t *testing.T) {
	file, want := concatStreams(t, concatParts()...)
	got, err := io.ReadAll(NewReader(bytes.NewReader(file)))
	if err != nil || string(got) != want {
		t.Fatalf("Read: got %d bytes (err %v), want %d", len(got), err, len(want))
	}
}

// ReadByte panicked here: nextBlock consumed the mid-file marker, returned nil
// with an empty buffer, and ReadByte indexed buf[0].
func TestReaderReadByteThroughMidFileEOF(t *testing.T) {
	file, want := concatStreams(t, concatParts()...)
	r := NewReader(bytes.NewReader(file))
	var got bytes.Buffer
	for {
		b, err := r.ReadByte()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		got.WriteByte(b)
	}
	if got.String() != want {
		t.Fatalf("ReadByte: got %d bytes, want %d", got.Len(), len(want))
	}
}

// IndexedReader.advance took the mid-file marker for the end, so a read that
// crossed the join -- a tabix chunk spanning header and body -- ended early.
func TestIndexedReaderReadsThroughMidFileEOF(t *testing.T) {
	file, want := concatStreams(t, concatParts()...)
	r := NewIndexedReader(bytes.NewReader(file))
	if err := r.SeekToVirtualOffset(0); err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(r)
	if err != nil || string(got) != want {
		t.Fatalf("Read: got %d bytes (err %v), want %d", len(got), err, len(want))
	}

	if err := r.SeekToVirtualOffset(0); err != nil {
		t.Fatal(err)
	}
	var byByte bytes.Buffer
	for {
		b, err := r.ReadByte()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		byByte.WriteByte(b)
	}
	if byByte.String() != want {
		t.Fatalf("ReadByte: got %d bytes, want %d", byByte.Len(), len(want))
	}
}

// Several joins in a row, and a stream that is nothing but its EOF marker.
func TestReadersThroughConsecutiveEmptyBlocks(t *testing.T) {
	file, want := concatStreams(t, "a\n", "", "", "b\n", "")
	if got, err := io.ReadAll(NewReader(bytes.NewReader(file))); err != nil || string(got) != want {
		t.Errorf("Reader: got %q (err %v), want %q", got, err, want)
	}
	r := NewIndexedReader(bytes.NewReader(file))
	r.SeekToVirtualOffset(0)
	if got, err := io.ReadAll(r); err != nil || string(got) != want {
		t.Errorf("IndexedReader: got %q (err %v), want %q", got, err, want)
	}
	ir := NewIndexingReader(bytes.NewReader(file))
	if got, err := io.ReadAll(ir); err != nil || string(got) != want {
		t.Errorf("IndexingReader: got %q (err %v), want %q", got, err, want)
	}
}

// IndexingReader stopped at the first empty block: it read only the header
// stream, and a .gzi built from it covered nothing past the join.
func TestIndexingReaderReadsThroughMidFileEOF(t *testing.T) {
	file, want := concatStreams(t, concatParts()...)
	ir := NewIndexingReader(bytes.NewReader(file))
	got, err := io.ReadAll(ir)
	if err != nil || string(got) != want {
		t.Fatalf("got %d bytes (err %v), want %d", len(got), err, len(want))
	}

	// Every recorded block must start with a real, non-empty block, at the
	// uncompressed offset recorded for it.
	blocks := ir.Blocks()
	if len(blocks) < 3 {
		t.Fatalf("want the body's blocks recorded too, got %d", len(blocks))
	}
	for _, b := range blocks {
		r := NewIndexedReader(bytes.NewReader(file))
		if err := r.SeekToVirtualOffset(NewVirtualOffset(b.Compressed, 0)); err != nil {
			t.Fatalf("block at %d: %v", b.Compressed, err)
		}
		buf := make([]byte, 16)
		n, _ := io.ReadFull(r, buf)
		if string(buf[:n]) != want[b.Uncompressed:b.Uncompressed+int64(n)] {
			t.Errorf("block at %d does not hold uncompressed offset %d", b.Compressed, b.Uncompressed)
		}
	}

	// And the .gzi made from those blocks seeks into the body.
	dir := t.TempDir()
	path := filepath.Join(dir, "concat.gz")
	if err := os.WriteFile(path, file, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := WriteGZIndex(path+".gzi", blocks); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	gzi, err := LoadGZIndex(path + ".gzi")
	if err != nil {
		t.Fatal(err)
	}
	r := NewIndexedReader(f)
	r.setGZIndex(gzi)
	off := int64(len(want) - 100)
	if _, err := r.Seek(off, io.SeekStart); err != nil {
		t.Fatal(err)
	}
	tail, err := io.ReadAll(r)
	if err != nil || string(tail) != want[off:] {
		t.Errorf("Seek(%d) via .gzi: got %q (err %v), want %q", off, tail, err, want[off:])
	}
}
