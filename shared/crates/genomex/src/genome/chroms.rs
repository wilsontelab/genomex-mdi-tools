//! Metadata about chromosomes as read from genome FASTA .fai index file.

// dependencies
use rustc_hash::FxHashMap;
use std::fs::read_to_string;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use mdi::OutputFile;

// constants
pub_key_constants!(
    GENOME
    GENOME_FASTA // from user options
    USE_ALL_CHROMS
    IS_COMPOSITE_GENOME
    GENOME_CHROMS // from shared/modules/set_genome_vars.sh
);

/// Chroms structure holding canonical chromosomes and their indices.
pub struct Chroms {
    genome:            String,
    pub is_composite_genome: bool,
    pub canonical:     Vec<String>,
    pub nuclear:       Vec<String>,
    pub canonical_indices: Vec<u8>, // in number order
    pub nuclear_indices:   Vec<u8>,
    pub index:         FxHashMap<String, u8>,
    pub rev_index:     FxHashMap<u8,     String>,
    pub nuclear_index: FxHashMap<String, u8>,
    pub sizes:         FxHashMap<String, u32>,
    pub index_sizes:   FxHashMap<u8,     u32>,
}
impl Chroms {
    /// Create a new Chroms instance from a genome FASTA .fai index file
    /// as identified from environment variables.
    pub fn new(cfg: &mut Config) -> Self {

        // load config variables
        cfg.set_string_env(&[GENOME, GENOME_FASTA, GENOME_CHROMS]);
        cfg.set_bool_env(  &[USE_ALL_CHROMS, IS_COMPOSITE_GENOME]);
        let is_composite_genome = *cfg.get_bool(IS_COMPOSITE_GENOME);

        // all placed chromosome sequences including chrM and chrEBV if present (but not chrXX_XX)
        // set by shared/modules/set_genome_vars.sh or similar script
        let chroms: Vec<String> = cfg.get_string(GENOME_CHROMS)
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();

        // parser for restricting work to properly ordered canonical chromosomes
        let canonical = if *cfg.get_bool(USE_ALL_CHROMS) || is_composite_genome {
            chroms // chroms used in fai indexed order
        } else {
            let is_roman = cfg.get_string(GENOME_CHROMS).to_uppercase().contains("CHRI");
            let candidates: Vec<String> = if is_roman { // handle sacCer3, etc.
                vec![
                    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                    "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
                    "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII", "XXVIII", "XXIX", "XXX",
                ].iter().map(|s| s.to_string()).collect()
            } else {
                let mut init: Vec<String> = (1..=90).map(|n| n.to_string()).collect();
                init.extend_from_slice(&["X".to_string(), "Y".to_string(), "M".to_string(), "EBV".to_string()]);
                let mut all: Vec<String> = vec![];
                for x in init { // handle Drosophila etc. L/R arms
                    all.push(x.clone());
                    all.push(format!("{}L", x));
                    all.push(format!("{}R", x));
                }
                all
            };
            let candidates: Vec<String> = candidates.into_iter().map(|s| format!("chr{}", s)).collect();
            let mut matched: Vec<String> = vec![];
            for candidate in candidates {
                if chroms.contains(&candidate) {
                    matched.push(candidate);
                }
            }
            matched
        };

        // nuclear chromosomes only
        let nuclear: Vec<String> = canonical.iter().filter(|chrom| {
            let uchrom = chrom.to_uppercase();
            !(uchrom.contains("CHRM") || uchrom.contains("CHREBV"))
        }).cloned().collect();

        // chromosome forward and reverse indices
        let mut index:         FxHashMap<String, u8> = FxHashMap::default();
        let mut rev_index:     FxHashMap<u8, String> = FxHashMap::default();
        let mut nuclear_index: FxHashMap<String, u8> = FxHashMap::default();
        let mut canonical_indices: Vec<u8> = vec![];
        let mut nuclear_indices:   Vec<u8> = vec![];
        for (i, chrom) in canonical.iter().enumerate() {
            let i1 = i as u8 + 1; // 1-referenced chrom indices, i.e., chr3 => 3
            index.insert(chrom.clone(), i1);
            rev_index.insert(i1, chrom.clone());
            canonical_indices.push(i1);
            if nuclear.contains(chrom) {
                nuclear_index.insert(chrom.clone(), i1);
                nuclear_indices.push(i1);
            }
        }
        index.insert("*".to_string(), 99); // special handling of unmapped reads
        rev_index.insert(99, "*".to_string());


        // chromosome sizes from .fai index
        let fai_file = format!("{}.fai", cfg.get_string(GENOME_FASTA));
        let mut sizes:       FxHashMap<String, u32> = FxHashMap::default();
        let mut index_sizes: FxHashMap<u8,     u32> = FxHashMap::default();
        if let Ok(fai_content) = read_to_string(&fai_file) {
            for line in fai_content.lines() {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    let chrom = parts[0].to_string();
                    if let Ok(size) = parts[1].parse() {
                        sizes.insert(chrom.clone(), size);
                        if let Some(&index) = index.get(&chrom) {
                            index_sizes.insert(index, size);
                        }
                    }
                }
            }
        }

        // return the Chroms instance
        Chroms {
            genome: cfg.get_string(GENOME).to_string(),
            is_composite_genome,
            canonical,
            nuclear,
            canonical_indices,
            nuclear_indices,
            index,
            rev_index,
            nuclear_index,
            sizes,
            index_sizes,
        }
    }

    /// Report if a chromosome name is a canonical chromosome.
    pub fn is_canonical(&self, chrom: &str) -> bool {
        self.index.contains_key( &chrom.to_string() )
    }

    /// Report if a chromosome index is a canonical chromosome.
    pub fn is_canonical_index(&self, chrom_index: u8) -> bool {
        self.rev_index.contains_key(&chrom_index)
    }

    /// Report if a chromosome name is a nuclear chromosome.
    pub fn is_nuclear(&self, chrom: &str) -> bool {
        self.nuclear_index.contains_key( &chrom.to_string() )
    }

    /// Report if a chromosome index is a nuclear chromosome.
    pub fn is_nuclear_index(&self, chrom_index: u8) -> bool {
        self.nuclear_index.contains_key(
            self.rev_index.get(&chrom_index).unwrap_or(&"*".to_string() )
        )
    }

    /// Get the genome suffix for composite genomes, e.g., "_hg38" or "_mm10".
    /// 
    /// For unaligned reads, return "*".
    pub fn get_genome_suffix(&self, chrom: &str) -> String {
        if self.is_composite_genome {
            chrom.split('_').rev().next().unwrap_or("*").to_string()
        } else {
            return self.genome.clone();
        }
    }
    /// Check if two alignments share the same genome suffix.
    pub fn is_same_genome_suffix(&self, chrom1: &str, chrom2: &str) -> bool {
        self.get_genome_suffix(chrom1) == self.get_genome_suffix(chrom2)
    }

    /// Write a chroms file with essential chromosome metadata.
    pub fn write_chroms_file(&self, filepath: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = OutputFile::open_file(
            filepath, 
            b'\t', 
            Some(&["chrom", "chrom_index1", "chrom_size"]),
        );
        for chrom_index1 in &self.canonical_indices {
            let chrom = self.rev_index.get(chrom_index1).unwrap();
            let chrom_size = self.index_sizes.get(chrom_index1).unwrap_or(&0);
            writer.write_record(vec![
                chrom,
                &chrom_index1.to_string(),
                &chrom_size.to_string(),
            ]);
        }
        Ok(())
    }
}
