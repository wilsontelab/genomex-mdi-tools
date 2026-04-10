//! Module for parsing and handling PAF alignment records.

// modules

// dependencies
use serde::{Serialize, Deserialize};
pub use crate::sam::tags::SamTags;

// structures
#[derive(Serialize, Deserialize)]
pub struct PafRecord { // a single PAF alignment record
    // read-level fields (query)
    pub qname:   String,
    pub qlen:    u32,
    // query alignment span  
    pub qstart0: u32,
    pub qend1:   u32,
    // alignment metadata
    pub strand:  char,
    // reference-level fields (target)   
    pub tname:   String,
    pub tlen:    u32,
    // query alignment span  
    pub tstart0: u32,
    pub tend1:   u32,
    // alignment metadata
    pub n_match: u32,
    pub n_bases: u32,
    pub mapq:    u8,
    // mixed read and alignment-level fields
    pub tags:  SamTags,
}
impl PafRecord {

}
