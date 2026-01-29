//! Suport manipulations and assessments related to gene annotations.

// dependencies
use std::ops::Deref;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use super::regions::*;

// constants
pub_key_constants!(
    // environment variables
    GENES_BED
);

/// Genes structure for tracking multiple genome target regions.
pub struct Genes(pub GenomeRegions); // a struct with one unnamed field
impl Deref for Genes { // call GenomeRegions methods directly on Genes
    type Target = GenomeRegions;
    fn deref(&self) -> &Self::Target { &self.0 }
}
impl Genes {
    /// Create a new Genes collection
    /// from a BED file provided as environment variable GENES_BED.
    /// 
    /// See GenomeRegions::from_env() for details.
    pub fn from_env(w: &mut Workflow, has_header: bool, gene_padding: u32) -> Self {
        Genes(
            GenomeRegions::from_env(
                GENES_BED, 
                &mut w.cfg,
                has_header, 
                gene_padding,
            )
        )
    }

    /// Create a new Genes collection
    /// from a BED file path provided as an argument.
    /// 
    /// See GenomeRegions::new() for details.
    pub fn new(
        bed_file:         &str, 
        has_header:       bool,
        gene_padding:     u32, 
    ) -> Self {
        Genes(
            GenomeRegions::new(
                bed_file, 
                has_header, 
                gene_padding,
                true,
            )
        )
    }
}