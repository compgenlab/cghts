package tabix

import "github.com/compgenlab/cghts/internal/vcfspan"

// Preset values for the tabix index header's `format` field, matching htslib's
// TBX_* constants. The low 16 bits hold the preset; bit 0x10000 is the zero-based
// (BED) coordinate flag.
const (
	PresetGeneric int32 = 0
	PresetSAM     int32 = 1
	PresetVCF     int32 = 2

	presetMask   int32 = 0xffff
	zeroBasedBit int32 = 0x10000
)

// vcfSpanEnd returns the 0-based half-open end of a VCF data line already split on
// tabs, whose 0-based start is beg.
//
// The rules live in internal/vcfspan so that this and the parsed-record path in the
// vcf package cannot drift: an index claiming a record spans 2 kb while a reader
// claims one base is a disagreement nothing would report. tabix cannot call the vcf
// package directly -- vcf imports tabix, not the reverse -- which is why the shared
// copy is internal rather than a method on VcfRecord.
func vcfSpanEnd(fields []string, beg int) int {
	return vcfspan.FieldsEnd(fields, beg)
}
