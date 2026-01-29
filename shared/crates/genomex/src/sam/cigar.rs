//! Module with CIGAR string regular expressions and utilities.
//! These methods mostly depend on just the CIGAR string itself.

// dependencies
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;
use regex::{Regex, Captures};

// regular expressions to process read clips
static CLIP_LEFT_LOCK: OnceLock<Regex> = OnceLock::new();
fn clip_left_regex() -> &'static Regex {
    CLIP_LEFT_LOCK.get_or_init(|| Regex::new(r"^(\d+)S").unwrap())
}
static CLIP_RIGHT_LOCK: OnceLock<Regex> = OnceLock::new();
fn clip_right_regex() -> &'static Regex {
    CLIP_RIGHT_LOCK.get_or_init(|| Regex::new(r"(\d+)S$").unwrap())
}
static CIGAR_OP_LOCK: OnceLock<Regex> = OnceLock::new();
fn cigar_op_regex() -> &'static Regex {
    CIGAR_OP_LOCK.get_or_init(|| Regex::new(r"(\d+)(\w)").unwrap())
}

/// SamFlag struct for working with SAM flag values.
#[derive(Serialize, Deserialize)]
pub struct CigarString{
    pub cigar: String,
}
impl CigarString {
    /// Create a new CigarString instance from a CIGAR string.
    pub fn new(cigar: &str) -> Self {
        CigarString {
            cigar: cigar.to_string(),
        }
    }
    /* -------------------------------------------------------------------------
    left and right alignment clips
    ------------------------------------------------------------------------- */
    /// Get the SAM left clip length from a CigarString.
    pub fn get_clip_left(&self) -> u32 {
        clip_left_regex()
            .captures(&self.cigar)
            .and_then(|caps| caps.get(1))
            .map_or(0, |m| m.as_str().parse::<u32>().unwrap_or(0))
    }
    /// Get the SAM right clip length from a CigarString.
    pub fn get_clip_right(&self) -> u32 {
        clip_right_regex()
            .captures(&self.cigar)
            .and_then(|caps| caps.get(1))
            .map_or(0, |m| m.as_str().parse::<u32>().unwrap_or(0))
    }
    /* -------------------------------------------------------------------------
    CIGAR string operations
        M: Alignment match (includes both sequence matches and mismatches).​
            =: Exact sequence match (distinguishes from mismatches in extended format).​
            X: Sequence mismatch (distinguishes from matches in extended format).
        I: Insertion to the reference (gap in reference, bases from read).​
        S: Soft-clipped sequence in read (present but unaligned).
        H: Hard-clipped sequence in read (removed from processing).​
        D: Deletion from the reference (gap in read).​​
        N: Skipped region on reference (e.g., introns in RNA-seq).​
    ------------------------------------------------------------------------- */
    /// Return a vector of CigarString operations as regex captures.
    /// Values are: `cap[1]` = length, `cap[2]` = operation.
    pub fn get_operations(&self) -> Vec<Captures<'_>> {
        cigar_op_regex().captures_iter(&self.cigar).collect()
    }
    /* ---------------------------------------------------------------------- */
    /// Check if a CIGAR operation consumes query bases as found in the original insert.
    /// This value includes soft clipped bases that are not part of the alignment
    /// as well as hard clipped bases removed from SEQ by the aligner.
    pub fn is_query_op_insert(op: &str) -> bool {
        op == "M" || op == "=" || op == "X" || op == "I" || op == "S" || op == "H"
    }
    /// Check if a CIGAR operation consumes query bases as found in SEQ.
    /// This value includes soft clipped bases that are not part of the alignment
    /// but excludes hard clipped bases not found in SEQ.
    pub fn is_query_op_seq(op: &str) -> bool {
        op == "M" || op == "=" || op == "X" || op == "I" || op == "S"
    }
    /// Check if a CIGAR operation consumes query bases in an alignment.
    /// This value excludes clipped bases that are not part of the alignment.
    pub fn is_query_op_aln(op: &str) -> bool {
        op == "M" || op == "=" || op == "X" || op == "I" 
    }
    /* ---------------------------------------------------------------------- */
    /// Check if a CIGAR operation consumes reference sequence.
    /// This value excludes query insertions and clipped bases.
    pub fn is_reference_op(op: &str) -> bool {
        op == "M" || op == "=" || op == "X" || op == "D" || op == "N"
    }
    /* -------------------------------------------------------------------------
    reference and query spans
    ------------------------------------------------------------------------- */
    /// Get the aligned size from a CigarString instance, i.e., the number of aligned bases.
    pub fn get_aligned_size(&self) -> u32 {
        let mut size = 0;
        for cap in self.get_operations() {
            if Self::is_query_op_aln(&cap[2]) {
                size += cap[1].parse().unwrap_or(0);
            }
        }
        size
    }
    /// Get the reference span from a CigarString instance, i.e., the number of reference bases spanned.
    pub fn get_ref_span(&self) -> u32 {
        let mut span = 0;
        for cap in self.get_operations() {
            if Self::is_reference_op(&cap[2]) {
                span += cap[1].parse().unwrap_or(0);
            }
        }
        span
    }
}
