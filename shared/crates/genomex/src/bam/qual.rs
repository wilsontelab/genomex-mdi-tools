//! Support functions for calculating against BAM Record QUAL fields.

// dependencies
use rust_htslib::bam::{Record as BamRecord};
use crate::bam::cigar::{get_clip_left, get_clip_right};

/// Determine whether QUAL is present, i.e., not '*'.
pub fn get_qual(aln: &BamRecord) -> Option<&[u8]> {
    let qual = aln.qual();
    if qual.len() > 0 && qual[0] != b'*' {
        Some(qual)
    } else {
        None
    }
}

/// Calculate the median base quality of the entire read in a BAM record. 
/// The resulting value does NOT have PHRED_OFFSET 33 added, it is the raw 
/// value from the BAM file. 
pub fn median_qual_all(aln: &BamRecord) -> u8 {
    if let Some(qual) = get_qual(aln) {
        median_qual_vec(qual.to_vec())
    } else {
        0
    }
}

/// Calculate the median base quality of the aligned portion of a BAM record. 
/// The resulting value does NOT have PHRED_OFFSET 33 added, it is the raw 
/// value from the BAM file. 
pub fn median_qual_aln(aln: &BamRecord) -> u8 {
    if let Some(qual) = get_qual(aln) {
        let start0 = get_clip_left(aln) as usize;
        let end1 = qual.len() - get_clip_right(aln) as usize;
        median_qual_vec(qual[start0..end1].to_vec())
    } else {
        0
    }
}

/// Calculate the median base quality of a QUAL as mut Vec<u8>.
fn median_qual_vec(mut qual: Vec<u8>) -> u8 {
    let mid = qual.len() / 2;
    let (_, median, _) = qual.select_nth_unstable(mid);
    *median
}

/// Calculate the average base quality of the entire read in a BAM record. 
/// The resulting value does NOT have PHRED_OFFSET 33 added, it is the raw 
/// value from the BAM file. 
pub fn avg_qual_all(aln: &BamRecord) -> f64 {
    if let Some(qual) = get_qual(aln) {
        avg_qual_vec(qual)
    } else {
        0.0
    }
}

/// Calculate the average base quality of the aligned portion of a BAM record. 
/// The resulting value does NOT have PHRED_OFFSET 33 added, it is the raw 
/// value from the BAM file. 
pub fn avg_qual_aln(aln: &BamRecord) -> f64 {
    if let Some(qual) = get_qual(aln) {
        let start0 = get_clip_left(aln) as usize;
        let end1 = qual.len() - get_clip_right(aln) as usize;
        avg_qual_vec(&qual[start0..end1])
    } else {
        0.0
    }
}

/// Calculate the average base quality of a QUAL as mut Vec<u8>.
fn avg_qual_vec(qual: &[u8]) -> f64 {
    if qual.is_empty() { return 0.0; }
    let sum: usize = qual.iter().map(|q| *q as usize).sum();
    sum as f64 / qual.len() as f64
}
