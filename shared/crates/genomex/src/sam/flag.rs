//! Module with SAM flag constants and utilities.
//! These methods depend on just the flag values themselves.

// dependencies
use serde::{Deserialize, Serialize};

// pub constants
pub const IS_PAIRED: u16      = 0b0000_0001; // 1
pub const PROPER_PAIR: u16    = 0b0000_0010; // 2
pub const UNMAPPED: u16       = 0b0000_0100; // 4
pub const MATE_UNMAPPED: u16  = 0b0000_1000; // 8
pub const REVERSE: u16        = 0b0001_0000; // 16
pub const MATE_REVERSE: u16   = 0b0010_0000; // 32
pub const FIRST_IN_PAIR: u16  = 0b0100_0000; // 64
pub const SECOND_IN_PAIR: u16 = 0b1000_0000; // 128
pub const SECONDARY: u16      = 0b0001_0000_0000; // 256
pub const FAILED_QC: u16      = 0b0010_0000_0000; // 512
pub const DUPLICATE: u16      = 0b0100_0000_0000; // 1024
pub const SUPPLEMENTAL: u16   = 0b1000_0000_0000; // 2048

/// SamFlag struct for working with SAM flag values.
#[derive(Serialize, Deserialize)]
pub struct SamFlag{
    pub flag: u16,
}
impl SamFlag {
    /// Create a new SamFlag instance from a u16 flag value.
    pub fn new(flag: u16) -> Self {
        SamFlag { flag }
    }
    /* -------------------------------------------------------------------------
    mask check methods
    ------------------------------------------------------------------------- */
    /// Check a SamFlag against one or more flag bits in a mask.
    /// Return true if any bits in the mask are set in the flag.
    pub fn check_any(&self, mask: u16) -> bool {
        (self.flag & mask) != 0
    }
    /// Check a SamFlag against one or more flag bits in a mask.
    /// Return true if all bits in the mask are set in the flag.
    pub fn check_all(&self, mask: u16) -> bool {
        (self.flag & mask) == mask
    }
    /// Check a SamFlag against one or more flag bits in a mask.
    /// Return true if none of the bits in the mask are set in the flag.
    pub fn check_none(&self, mask: u16) -> bool {
        (self.flag & mask) == 0
    }
}
