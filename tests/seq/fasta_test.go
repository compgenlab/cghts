package seq_test

// import (
// 	"io"
// 	"strings"
// 	"testing"

// 	"github.com/compgen-io/cgltk/seqio"
// )

// func TestFastaReaderChunking(t *testing.T) {
// 	input := strings.Join([]string{
// 		">seq1",
// 		"ACGTACGT",
// 		"",
// 		">seq2",
// 		"TTTT",
// 		"",
// 	}, "\n")
// 	reader := seqio.NewFastaReader(strings.NewReader(input), seqio.NewSeqReaderOptions().SetBufSize(4))

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
// 		{seq: "ACGT", more: true},
// 		{seq: "ACGT", more: false},
// 		{seq: "TTTT", more: false},
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
