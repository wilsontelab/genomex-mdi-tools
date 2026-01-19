use crate::sam::flag::SamFlag;

/// Support SAM/BAM file size control and other data clearing
/// by allowing certain fields to be set to null values.

// dependencies
use super::{SamRecord, SamQual, SamTags};

// constants
pub const FLAG: &str  = "flag"; // nullable fields
pub const RNEXT: &str = "rnext";
pub const PNEXT: &str = "pnext";
pub const TLEN: &str  = "tlen";
pub const SEQ: &str   = "seq";
pub const QUAL: &str  = "qual";
pub const TAGS: &str  = "tags";

// implemenation
impl SamRecord {
    /* -------------------------------------------------------------------------
    record size and data content management
    ------------------------------------------------------------------------- */
    /// Simplify the requested fields of a SamRecord to null values.
    pub fn set_to_null(&mut self, fields: &[&str]) {
        for field in fields {
            match *field {
                FLAG  => self.flag  = SamFlag::new(0),
                RNEXT => self.rnext = "*".to_string(),
                PNEXT => self.pnext = 0,
                TLEN  => self.tlen  = 0,
                SEQ   => self.seq   = "*".to_string(),
                QUAL  => self.qual  = SamQual::new("*"),
                TAGS  => self.tags  = SamTags::new(vec![]),
                _ => panic!("SamRecord::set_to_null unsupported field: {}", field),
            }
        }
    }
}
