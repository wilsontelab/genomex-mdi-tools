//! Module with SAM QUAL utilities.
//! These methods depend on just the QUAL values themselves.

// dependencies
use serde::{Deserialize, Serialize};

// constants
const PHRED_OFFSET: u8 = 33;

/// SamQual struct for working with SAM QUAL values.
#[derive(Serialize, Deserialize)]
pub struct SamQual{
    pub qual: String,
}
impl SamQual {
    /// Create a new SamQual instance from a QUAL string.
    pub fn new(qual: &str) -> Self {
        SamQual { qual: qual.to_string() }
    }
    /* -------------------------------------------------------------------------
    quality score averaging
    ------------------------------------------------------------------------- */
    /// Get the average per-base Phred QUAL for a SamQual string.
    pub fn get_avg_qual(&self) -> f64 {
        Self::get_avg_qual_str(&self.qual)
    }
    /// Get the average per-base Phred QUAL from a QUAL string.
    /// May be just a portion of a read as determined by the caller.
    pub fn get_avg_qual_str(qual: &str) -> f64 {
        if qual.is_empty() { return 0.0; }
        let sum: usize = qual.bytes().map(|b| (b as usize).saturating_sub(33)).sum();
        sum as f64 / qual.len() as f64
    }
    /* -------------------------------------------------------------------------
    quality score quanization to reduce bam file size
    ------------------------------------------------------------------------- */
    /// Quantize quality scores to 8 levels to reduce BAM file size, using levels
    /// as defined in: 
    ///      https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac054/6661371
    ///
    /// This is an unsafe function because it modifies the input string in place,
    /// which creates undefined behavior if the input string is not valid ASCII.
    /// Please ensure that the input string is a valid SAM QUAL string.
    pub unsafe fn quantize_qual_scores(qual_str: &mut str) {
        unsafe { for byte in qual_str.as_bytes_mut().iter_mut() {
            let x: u8 = match byte.saturating_sub(PHRED_OFFSET) {
                0..=6 => 5,
                7..=11 => 10,
                12..=16 => 15,
                17..=21 => 20,
                22..=26 => 25,
                27..=31 => 30,
                32..=36 => 35,
                _ => 40,
            };
            *byte = x + PHRED_OFFSET; // safe because we are converting from ASCII to ASCII
        } }
    }
}
