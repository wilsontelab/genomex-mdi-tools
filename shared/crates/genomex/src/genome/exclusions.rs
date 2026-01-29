//! Suport manipulations and assessments related to excluded 
//! genome regions, e.g., the ENCODE blacklist regions, or
//! any other set of regions to be excluded from analysis.

// dependencies
use std::ops::Deref;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use super::regions::*;

// constants
pub_key_constants!(
    // environment variables
    GENOME_EXCLUSIONS_BED
);

/// Exclusions structure for tracking multiple genome exclusion regions.
pub struct Exclusions(pub GenomeRegions); // a struct with one unnamed field
impl Deref for Exclusions { // call GenomeRegions methods directly on Exclusions
    type Target = GenomeRegions;
    fn deref(&self) -> &Self::Target { &self.0 }
}
impl Exclusions {
    /// Create a new Exclusions collection
    /// from a BED file provided as environment variable GENOME_EXCLUSIONS_BED.
    /// 
    /// See GenomeRegions::from_env() for details.
    pub fn from_env(w: &mut Workflow, has_header: bool) -> Self {
        Exclusions(
            GenomeRegions::from_env(
                GENOME_EXCLUSIONS_BED, 
                &mut w.cfg,
                has_header, 
                0, // no padding for exclusions, actual regions are expected to be precise
            )
        )
    }

    /// Create a new Exclusions collection
    /// from a BED file path provided as an argument.
    /// 
    /// See GenomeRegions::new() for details.
    pub fn new(
        bed_file:         &str, 
        has_header:       bool,
    ) -> Self {
        Exclusions(
            GenomeRegions::new(
                bed_file, 
                has_header, 
                0,
                true,
            )
        )
    }
}
