// Package seqanalysis provides sequence analysis routines for FASTA/FASTQ
// records.
//
// It operates on the seqio.SeqRecord interface and streams sequence data in
// chunks so that arbitrarily large records can be analyzed without loading the
// whole sequence into memory. Currently it provides [CalcGC] for computing GC
// content.
package seqanalysis
