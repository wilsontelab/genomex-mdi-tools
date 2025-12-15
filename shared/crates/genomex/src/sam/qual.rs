//! Module with SAM QUAL utilities.
//! These methods depend on just the QUAL values themselves.

// dependencies
use serde::{Deserialize, Serialize};

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
}
