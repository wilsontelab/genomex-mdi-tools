//! Support functions for manipulating BAM CIGAR strings.

// dependencies
use rust_htslib::bam::{Record as BamRecord};

/* -------------------------------------------------------------------------
left and right alignment clips
------------------------------------------------------------------------- */
/// Get the BAM left clip length from a CIGAR string.
pub fn get_clip_left(aln: &BamRecord) -> u32 {
    aln.cigar().leading_softclips() as u32
}
/// Get the BAM right clip length from a CIGAR string.
pub fn get_clip_right(aln: &BamRecord) -> u32 {
    aln.cigar().trailing_softclips() as u32
}
/* -------------------------------------------------------------------------
query start and end
------------------------------------------------------------------------- */
/// Get the PAF-like start of a SamRecord's alignment on the query sequence (0-based).
pub fn get_query_start0(aln: &BamRecord) -> u32 {
    if aln.is_unmapped() { return 0; }
    if aln.is_reverse() {
        get_clip_right(aln)
    } else {
        get_clip_left(aln)
    }
}
/// Get the PAF-like end of a SamRecord's alignment on the query sequence (1-based).
pub fn get_query_end1(aln: &BamRecord, qlen: u32) -> u32 {
    if aln.is_unmapped() { return qlen; }
    if aln.is_reverse() {
        qlen - get_clip_left(aln)
    } else {
        qlen - get_clip_right(aln)
    }
}
