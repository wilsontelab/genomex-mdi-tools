//! Module for parsing and handling SAM alignment records.
//! Provides access to methods specific to flags, CIGAR strings, and tags,
//! as well as methods that act on combinations of values.

// modules
pub mod flag;
pub mod cigar;
pub mod qual;
pub mod tags;
pub mod junction;

// dependencies
use std::str::FromStr;
use std::fmt;
use serde::{Deserialize, Serialize};
use crate::sequence::rc_acgtn_str;
use flag::SamFlag;
use cigar::CigarString;
pub use qual::SamQual;
use tags::SamTags;

// structures
#[derive(Serialize, Deserialize)]
pub struct SamRecord { // a single SAM alignment record
    pub qname: String,
    pub flag:  SamFlag,
    pub rname: String,
    pub pos1:  usize, 
    pub mapq:  u8,
    pub cigar: CigarString,
    pub rnext: String,
    pub pnext: usize,
    pub tlen:  i32,
    pub seq:   String, // when streamed from STDIN, cannot used borrowed fields, must be DeserializeOwned
    pub qual:  SamQual,
    pub tags:  SamTags,
}
impl SamRecord {
    /* -------------------------------------------------------------------------
    SamFlag methods
    ------------------------------------------------------------------------- */
    /// Check a SamRecord's flag against one or more flag bits in a mask.
    /// Return true if any bits in the mask are set in the flag.
    pub fn check_flag_any(&self, mask: u16) -> bool {
        self.flag.check_any(mask)
    }
    /// Check a SamRecord's flag against one or more flag bits in a mask.
    /// Return true if all bits in the mask are set in the flag.
    pub fn check_flag_all(&self, mask: u16) -> bool {
        self.flag.check_all(mask)
    }
    /// Check a SamRecord's flag against one or more flag bits in a mask.
    /// Return true if none of the bits in the mask are set in the flag.
    pub fn check_flag_none(&self, mask: u16) -> bool {
        self.flag.check_none(mask)
    }    
    /* -------------------------------------------------------------------------
    CigarString methods
    ------------------------------------------------------------------------- */
    /// Get the left clip length from a SamRecord's CIGAR string.
    pub fn get_clip_left(&self) -> usize {
        self.cigar.get_clip_left()
    }
    /// Get the right clip length from a SamRecord's CIGAR string.
    pub fn get_clip_right(&self) -> usize {
        self.cigar.get_clip_right()
    }
    /* ---------------------------------------------------------------------- */
    /// Get the PAF-like start of a SamRecord's alignment on the query sequence (0-based).
    pub fn get_query_start0(&self) -> usize {
        if self.check_flag_any(flag::UNMAPPED) { return 0; }
        if self.check_flag_any(flag::REVERSE) {
            self.get_clip_right()
        } else {
            self.get_clip_left()
        }
    }
    /// Get the PAF-like end of a SamRecord's alignment on the query sequence (1-based).
    pub fn get_query_end1(&self) -> usize {
        let qlen = self.seq.len();
        if self.check_flag_any(flag::UNMAPPED) { return qlen; }
        if self.check_flag_any(flag::REVERSE) {
            qlen - self.get_clip_left()
        } else {
            qlen - self.get_clip_right()
        }
    }
    /* ---------------------------------------------------------------------- */
    /// Get the rightmost mapped read position in the reference genome from
    /// a SamRecord's alignment (1-based).
    pub fn get_end1(&self) -> usize {
        let mut end1 = self.pos1 - 1;
        for cap in  self.cigar.get_operations() {
            if CigarString::is_reference_op(&cap[2]) {
                end1 += cap[1].parse().unwrap_or(0);
            }
        }
        end1
    }
    /* ---------------------------------------------------------------------- */
    /// Get the aligned size from a SamRecord's CIGAR string, i.e., the number of aligned bases.
    pub fn get_aligned_size(&self) -> usize {
        self.cigar.get_aligned_size()
    }
    /// Get the reference span from a SamRecord's CIGAR string, i.e., the number of reference bases spanned.
    pub fn get_ref_span(&self) -> usize {
        self.cigar.get_ref_span()
    }
    /* -------------------------------------------------------------------------
    SEQ methods
    ----||||||||||--          cL=4, cR=2, alen=10, qlen=16, rS1=5, rE1=14, qS0=4, qE1=14
    =========================
      --||||||||||----        cL=2, cR=4, alen=10, qlen=16, rS1=5, rE1=14, qS0=4, qE1=14
    ------------------------------------------------------------------------- */
    /// Get a SamRecord's SEQ for the entire read in the original sequenced orientation. 
    /// The reference-oriented SEQ is simply the value of SamRecord.seq.
    pub fn get_seq_read(&self) -> String {
         if self.check_flag_any(flag::REVERSE){
            rc_acgtn_str(&self.seq)
        } else {
            self.seq.clone()
        }
    }
    /// Get a SamRecord's SEQ for the aligned portion of the read, if any, as Option<String>.
    /// If as_sequenced==true, SEQ is returned in the original sequenced orientation,
    /// otherwise it is reference-oriented and may have been reverse-complemented by the aligner.
    pub fn get_seq_aln(&self, as_sequenced: bool) -> Option<String> {
        if self.check_flag_any(flag::UNMAPPED) { return None; }
        let qstart0 = self.get_query_start0();
        let qend1 = self.get_query_end1();
        if self.check_flag_any(flag::REVERSE){
            let qlen = self.seq.len();
            let seq = &self.seq[qlen - qend1..qlen - qstart0];
            if as_sequenced {
                Some(rc_acgtn_str(&seq))
            } else {
                Some(seq.to_string())
            }
        } else {
            Some(self.seq[qstart0..qend1].to_string())
        }
    }
    /* -------------------------------------------------------------------------
    SamQual methods
    ------------------------------------------------------------------------- */
    /// Get the average per-base Phred QUAL for the entire read in a SamRecord
    /// (not just the aligned portion). 
    pub fn get_avg_qual_read(&self) -> f64 {
        self.qual.get_avg_qual()
    }
    /// Get the average per-base Phred QUAL for the aligned portion of a SamRecord.
    pub fn get_avg_qual_aln(&self) -> f64 {
        if self.check_flag_any(flag::UNMAPPED) { return 0.0; }
        let qstart0 = self.get_query_start0();
        let qend1 = self.get_query_end1();
        if self.check_flag_any(flag::REVERSE){
            let qlen = self.seq.len();
            return SamQual::get_avg_qual_str(&self.qual.qual[qlen - qend1..qlen - qstart0]);
        } else {
            return SamQual::get_avg_qual_str(&self.qual.qual[qstart0..qend1]);
        }
    }
    /* -------------------------------------------------------------------------
    SamTags methods
    ------------------------------------------------------------------------- */
    /// Extract a parsed tag value from a SamRecord's tags as Some(T)
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., `"AS:"` or `"AS:i:"`. Return None if the tag is not present.
    pub fn get_tag_value_parsed<T>(&self, prefix: &str) -> Option<T>
    where
        T: FromStr,
        <T as FromStr>::Err: fmt::Debug,
    {
        self.tags.get_tag_value_parsed(prefix)
    }
    /// Extract a parsed tag value from a SamRecord's tags as Some(String),
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., `"AS:"` or `"AS:i:"`. Return None if the tag is not present.
    pub fn get_tag_value(&self, prefix: &str) -> Option<String>{
        self.get_tag_value_parsed::<String>(prefix)
    }
    /// Retain only those tags whose prefix is in the provided &[&str],
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., either `"AS:"` or `"AS:i:"`.
    /// 
    /// Tags on the retention list not found in the SamTags are simply ignored.
    pub fn retain_tags(&mut self, prefixes: &[&str]) {
        self.tags.retain(prefixes);
    }
}
