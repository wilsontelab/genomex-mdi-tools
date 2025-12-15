//! Metadata about chromosomes as read from genome FASTA .fai index file.

// dependencies
use std::collections::HashMap;
use std::fs::read_to_string;
use mdi::pub_key_constants;
use mdi::workflow::Config;

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
    pub index:         HashMap<String, usize>,
    pub rev_index:     HashMap<usize, String>,
    pub nuclear_index: HashMap<String, usize>,
    pub sizes:         HashMap<String, usize>,
    pub index_sizes:   HashMap<usize, usize>,
}
impl Chroms {
    /// Create a new Chroms instance from a genome FASTA .fai index file
    /// as identified from environment variables.
    pub fn new(cfg: &mut Config) -> Self {

        // load config variables
        cfg.set_string_env(&[GENOME, GENOME_FASTA, GENOME_CHROMS]);
        cfg.set_bool_env(&[USE_ALL_CHROMS, IS_COMPOSITE_GENOME]);
        let is_composite_genome = *cfg.get_bool(IS_COMPOSITE_GENOME);

        // all placed chromosome sequences including chrM and chrEBV if present (but not chrXX_XX)
        // set by shared/modules/set_genome_vars.sh or similar script
        let chroms: Vec<String> = cfg.get_string(GENOME_CHROMS).split(" ").map(|s| s.to_string()).collect();

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
        let mut index:         HashMap<String, usize> = HashMap::new();
        let mut rev_index:     HashMap<usize, String> = HashMap::new();
        let mut nuclear_index: HashMap<String, usize> = HashMap::new();
        for (i, chrom) in canonical.iter().enumerate() {
            index.insert(chrom.clone(), i + 1); // 1-referenced chrom indices, i.e., chr3 => 3
            rev_index.insert(i + 1, chrom.clone());
            if nuclear.contains(chrom) {
                nuclear_index.insert(chrom.clone(), i + 1);
            }
        }
        index.insert("*".to_string(), 99); // special handling of unmapped reads
        rev_index.insert(99, "*".to_string());


        // chromosome sizes from .fai index
        let fai_file = format!("{}.fai", cfg.get_string(GENOME_FASTA));
        let mut sizes:       HashMap<String, usize> = HashMap::new();
        let mut index_sizes: HashMap<usize, usize>  = HashMap::new();
        if let Ok(fai_content) = read_to_string(&fai_file) {
            for line in fai_content.lines() {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    let chrom = parts[0].to_string();
                    if let Ok(size) = parts[1].parse::<usize>() {
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
    pub fn is_canonical_index(&self, chrom_index: usize) -> bool {
        self.rev_index.contains_key(&chrom_index)
    }

    /// Report if a chromosome name is a nuclear chromosome.
    pub fn is_nuclear(&self, chrom: &str) -> bool {
        self.nuclear_index.contains_key( &chrom.to_string() )
    }

    /// Report if a chromosome index is a nuclear chromosome.
    pub fn is_nuclear_index(&self, chrom_index: usize) -> bool {
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
}

// sub writeChromsFile {
//     my ($chromsFile, $genomeFasta) = @_;
//     @canonicalChroms or setCanonicalChroms();
//     my @chromSizes = getChromIndexSizes("$genomeFasta.fai");
//     open my $chrH, ">", $chromsFile or die "could not open: $chromsFile: $!\n";
//     foreach my $chrom(keys %chromIndex){
//         my $chromIndex1 = $chromIndex{$chrom};
//         print $chrH join("\t", $chrom, $chromIndex1, $chromSizes[$chromIndex1] || "NA"), "\n";
//     }
//     close $chrH;
// }
