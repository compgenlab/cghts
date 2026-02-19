package seq_test

// import (
// 	"io"
// 	"strings"
// 	"testing"

// 	"github.com/compgen-io/cgkmer/seq"
// )

// func TestFastqReaderChunking(t *testing.T) {
// 	input := strings.Join([]string{
// 		"@read1",
// 		"ACGTACGT",
// 		"+",
// 		"!!!!!!!!",
// 		"@read2",
// 		"TT",
// 		"+",
// 		"##",
// 		"",
// 	}, "\n")
// 	reader := seq.NewFastqReader(strings.NewReader(input), seq.NewSeqReaderOptions().SetBufSize(3))

// 	type chunk struct {
// 		seq  string
// 		more bool
// 	}
// 	var got []chunk
// 	for {
// 		seq, more, err := reader.ReadSeq()
// 		if err == io.EOF {
// 			break
// 		}
// 		if err != nil {
// 			t.Fatalf("ReadSeq error: %v", err)
// 		}
// 		got = append(got, chunk{seq: seq, more: more})
// 	}

// 	want := []chunk{
// 		{seq: "ACG", more: true},
// 		{seq: "TAC", more: true},
// 		{seq: "GT", more: false},
// 		{seq: "TT", more: false},
// 	}
// 	if len(got) != len(want) {
// 		t.Fatalf("expected %d chunks, got %d", len(want), len(got))
// 	}
// 	for i := range want {
// 		if got[i] != want[i] {
// 			t.Fatalf("chunk %d: expected %+v, got %+v", i, want[i], got[i])
// 		}
// 	}
// }

// func TestFastqReaderInvalidHeader(t *testing.T) {
// 	input := strings.Join([]string{
// 		"read1",
// 		"ACGT",
// 		"+",
// 		"!!!!",
// 		"",
// 	}, "\n")
// 	reader := seq.NewFastqReader(strings.NewReader(input), seq.NewSeqReaderOptions().SetBufSize(4))

// 	_, _, err := reader.ReadSeq()
// 	if err == nil {
// 		t.Fatalf("expected error for invalid FASTQ header")
// 	}
// }

// func TestFastqReaderInvalidQualityLength(t *testing.T) {
// 	input := strings.Join([]string{
// 		"@read1",
// 		"ACGT",
// 		"+",
// 		"!!",
// 		"",
// 	}, "\n")
// 	reader := seq.NewFastqReader(strings.NewReader(input), seq.NewSeqReaderOptions().SetBufSize(4))

// 	_, _, err := reader.ReadSeq()
// 	if err == nil {
// 		t.Fatalf("expected error for invalid FASTQ quality length")
// 	}
// }
