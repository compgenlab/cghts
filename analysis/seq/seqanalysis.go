package seqanalysis

import "github.com/compgen-io/cgkit/seqio"

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
