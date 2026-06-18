package seqanalysis

import "github.com/compgenlab/hts/seqio"

// CalcGC returns the GC content of a sequence record as a fraction in the range
// [0, 1]. It is the count of G and C bases (case-insensitive) divided by the
// total number of bases. The sequence is streamed in chunks, so records of any
// length can be processed without loading the whole sequence into memory. An
// empty sequence returns 0.
func CalcGC(s seqio.SeqRecord) float64 {
	gcCount := 0
	total := 0
	for chunk := range s.Chunks(1024) {
		for _, base := range chunk.Seq() {
			switch base {
			case 'G', 'C', 'g', 'c':
				gcCount++
			}
			total++
		}
	}

	if total > 0 {
		return (float64(gcCount) / float64(total))
	}
	return 0.0
}
