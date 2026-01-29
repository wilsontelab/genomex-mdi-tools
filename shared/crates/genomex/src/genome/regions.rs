//! Suport manipulations and assessment related to 
//! genome query regions, i.e., spans.
//! 
//! Specific NewTypes are provided for Genes, Exclusions, 
//! and TargetRegions
//! 
//! At present, regions and comparisons are unstranded.

// dependencies
use rustc_hash::FxHashMap;
use serde::{Serialize, Deserialize};
use mdi::workflow::Config;
use mdi::InputFile;
use crate::genome::chroms::Chroms;

// constants
// region-matching states for single positions
pub const INSIDE_REGION: &str    = "I"; // "Inside region"
pub const NEAR_REGION: &str      = "A"; // "Adjacent to region"
pub const OUTSIDE_REGION: &str   = "O"; // "Outside region"
// region-matching states for pairs of positions
pub const IN_IN_REGION: &str     = "II";
pub const OUT_OUT_REGION: &str   = "OO";
// null region constants
pub const NULL_REGION_I: u32     = 0;
pub const NULL_REGION_NAME: &str = "*";
pub const NULL_DISTANCE: i32     = 0;

/// A GenomeRegion describes one span in a reference genome.
/// 
/// Name and center fields are Vec to be able to aggregate
/// regions together, e.g., when splitting overlapping regions.
#[derive(Clone)]
pub struct GenomeRegion {
    pub index1:  u32,    // 1-referenced region index; with 3-bit shift allows up to 536,870,911 regions
    pub padded:  bool,
    pub chrom:   String,
    pub start0:  u32,    // 0-based inclusive
    pub end1:    u32,    // 1-based inclusive (BED convention)
    pub names:   Vec<String>, // the gene or other text identifier of a region; "*" if not specified in BED
    pub centers: Vec<u32>,
}

/// GenomeRegions collects multiple GenomeRegion instances.
pub struct GenomeRegions {
    pub has_data:               bool, // flag whether any regions were loaded into this GenomeRegions collection
    pub bed_file:               String,
    pub has_header:             bool,
    pub region_padding:         u32,
    pub n_regions:              u32,
    pub sum_region_lens:        u32, // u32 exceeds the total length of most genomes
    pub sum_padded_region_lens: u32,
    pub region_chroms:          Vec<String>, // lists of the properties of all regions
    pub region_start0s:         Vec<u32>,
    pub region_end1s:           Vec<u32>,
    pub region_names:           Vec<Vec<String>>,
    pub region_centers:         Vec<Vec<u32>>,
    pub unique_chroms:          Vec<String>, // unique chromosomes included in regions
    pub regions:                FxHashMap<String, Vec<GenomeRegion>>, // unique_chroms -> sorted regions
    pub padded_regions:         FxHashMap<String, Vec<GenomeRegion>>,
    pub paired_states:          FxHashMap<String, u8>,
    pub paired_state_states:    Vec<&'static str>,
    pub null_region:            GenomeRegion,
}
impl GenomeRegions {

    /// Create a new GenomeRegions collection by loading regions 
    /// from a BED file provided as an environment variable key.
    /// 
    /// Use 0 for `region_padding` to load only actual regions;
    /// use >0 to also load padded regions.
    /// 
    /// If bed_file_env_key does not match an environment variable
    /// with a valid BED file path, no regions are loaded and
    /// self.has_data is set to false.
    pub fn from_env(
        bed_file_env_key: &str, 
        cfg:              &mut Config,
        has_header:       bool,
        region_padding:   u32, 
    ) -> Self {
        cfg.set_bool_env(  &[bed_file_env_key]);
        cfg.set_string_env(&[bed_file_env_key]);
        let has_data = *cfg.get_bool(bed_file_env_key);
        Self::new(
            cfg.get_string(bed_file_env_key), 
            has_header, 
            region_padding,
            has_data,
        )
    }

    /// Create a new GenomeRegions collection by loading regions 
    /// from a BED file path provided as an argument.
    /// 
    /// Use 0 for `region_padding` to load only actual regions;
    /// use >0 to also load padded regions.
    /// 
    /// Typically `has_data` is set to true if this function
    /// is called directly. See also `GenomeRegions::from_env()`.
    pub fn new(
        bed_file:         &str, 
        has_header:       bool,
        region_padding:   u32, 
        has_data:         bool,
    ) -> Self {

        // load config variables
        // cfg.set_bool_env(&[bed_file_env_key]);
        // cfg.set_string_env(&[bed_file_env_key]);
        // let has_data = *cfg.get_bool(bed_file_env_key);
        let paired_states: FxHashMap<String, u8> = [
            ("OO".to_string(), 0), // outside-outside region pair, etc.
            ("II".to_string(), 1),
            ("IA".to_string(), 2),
            ("IO".to_string(), 3),
            ("AA".to_string(), 4),
            ("AO".to_string(), 5)  // max paired_state is 5, fits in 3 bits
        ].into_iter().collect();
        let paired_state_states = vec![
            OUTSIDE_REGION, // for mapping paired_states to single region matching states
            INSIDE_REGION,
            INSIDE_REGION,
            INSIDE_REGION,  // e.g., IO -> IN_REGION, etc.
            NEAR_REGION,
            NEAR_REGION
        ];

        // initialize GenomeRegions instance
        let mut regions = GenomeRegions {
            has_data,
            bed_file: bed_file.to_string(),
            has_header,
            region_padding,
            n_regions: 0,
            sum_region_lens: 0,
            sum_padded_region_lens: 0,
            region_chroms:  vec![],
            region_start0s: vec![],
            region_end1s:   vec![],
            region_names:   vec![],
            region_centers: vec![],
            unique_chroms:  vec![],
            regions: FxHashMap::default(),
            padded_regions: FxHashMap::default(),
            paired_states,
            paired_state_states,
            null_region: GenomeRegion{
                index1:  NULL_REGION_I,
                padded:  false,
                chrom:   NULL_REGION_NAME.to_string(),
                start0:  0,
                end1:    0,
                names:   vec![NULL_REGION_NAME.to_string()],
                centers: vec![0],
            }
        };  

        // load genome regions, if any
        if has_data {

            // first pass to load unpadded regions
            regions.load_regions(0);

            // second pass to load padded regions to support the A=adjacent/near state
            if region_padding > 0 {
                regions.load_regions(region_padding);
            } else {
                regions.sum_padded_region_lens = regions.sum_region_lens;
            }
        }

        // return the GenomeRegions instance
        regions
    }

    // load the genome regions from the requested BED file
    fn load_regions(&mut self, padding: u32) {
        let mut index1 = 0;
        let mut sum_region_lens: u32 = 0;
        for line in InputFile::get_lines(&self.bed_file, self.has_header) {
            index1 += 1; // 1-referenced, so the index1 == 0 means no region
            let line = line.replace("\r", "");
            let parts: Vec<&str> = line.split('\t').collect();
            let chrom = parts[0];
            let mut start0: u32 = parts[1].parse().unwrap_or(0);
            let mut end1:   u32 = parts[2].parse().unwrap_or(0);
            let name = if parts.len() >= 4 { parts[3] } else { NULL_REGION_NAME };
            sum_region_lens += end1 - start0;
            let center = ((start0 as f64 + 1.0 + end1 as f64) / 2.0) as u32;
            if padding == 0 {
                self.n_regions += 1;
                self.region_chroms .push(chrom.to_string());
                self.region_start0s.push(start0);
                self.region_end1s  .push(end1);
                self.region_names  .push(vec![name.to_string()]);
                self.region_centers.push(vec![center]);
                if !self.unique_chroms.contains(&chrom.to_string()) {
                    self.unique_chroms.push(chrom.to_string());
                }
            } else {
                start0 = start0.saturating_sub(padding);
                end1 += padding;
            }
            let region = GenomeRegion {
                index1,
                padded:  padding > 0,
                chrom:   chrom.to_string(),
                start0,
                end1,
                names:   vec![name.to_string()],
                centers: vec![center],
            };
            if padding == 0 {
                self.regions.entry(chrom.to_string()).or_insert_with(Vec::new).push(region);
            } else {
                self.padded_regions.entry(chrom.to_string()).or_insert_with(Vec::new).push(region);
            }
        }

        // sort regions by start position within each chromosome
        if padding == 0 {
            for regions in self.regions.values_mut() {
                regions.sort_by_key(|r| r.start0);
            }
            self.sum_region_lens = sum_region_lens;
        } else {
            for regions in self.padded_regions.values_mut() {
                regions.sort_by_key(|r| r.start0);
            }
            self.sum_padded_region_lens = sum_region_lens;
        }
    }

    /// Create a new, empty GenomeRegions collection using
    /// a parent GenomeRegions collection as a template.
    pub fn from_template(&self) -> Self {
        GenomeRegions {
            has_data:               self.has_data, // some values won't change in any derived collection
            bed_file:               self.bed_file.clone(),
            has_header:             self.has_header,
            region_padding:         self.region_padding,
            n_regions:              0, // other values reset to null pending region modification or aggregations
            sum_region_lens:        0,
            sum_padded_region_lens: 0,
            region_chroms:          vec![],
            region_start0s:         vec![],
            region_end1s:           vec![],
            region_names:           vec![],
            region_centers:         vec![],
            unique_chroms:          self.unique_chroms.clone(),
            regions:                FxHashMap::default(),
            padded_regions:         FxHashMap::default(),
            paired_states:          self.paired_states.clone(),
            paired_state_states:    self.paired_state_states.clone(),
            null_region:            self.null_region.clone(),
        }
    }
    
    /// Split regions into distinct overlap regions that aggregate all 
    /// input regions that overlap each output region. Takes ownership
    /// of the input GenomeRegions to return a new GenomeRegions collection.
    /// 
    /// Currently does not handle padding of split regions; padding
    /// is always set to 0 with unpadded regions as output.
    pub fn split_overlaps(self) -> GenomeRegions {
        let mut genome_regions = self.from_template();
        genome_regions.region_padding = 0; // split regions are unpadded
    
        // process each chromosome separately
        for (chrom, regions_in) in self.regions.iter() {
            let mut index1_out: u32 = 0;
            genome_regions.regions.insert(chrom.clone(), vec![]);
            let regions_out = genome_regions.regions.get_mut(chrom).unwrap();

            // group regions into chains of overlaps
            let mut i = 0;
            while i < regions_in.len() {

                // initialize the first/next group after an overlap break
                let mut overlap_group: Vec<&GenomeRegion> = vec![&regions_in[i]];
                let mut group_end1 = regions_in[i].end1;

                // extend the group while subsequent regions overlap
                let mut j = i + 1;
                while j < regions_in.len() && regions_in[j].start0 < group_end1 {
                    overlap_group.push(&regions_in[j]);
                    if regions_in[j].end1 > group_end1 { group_end1 = regions_in[j].end1; }
                    j += 1;
                }

                // no overlaps, add the single region as-is
                if overlap_group.len() == 1 {
                    index1_out += 1;
                    let region_in = overlap_group[0];
                    let size = region_in.end1 - region_in.start0;
                    let region_out = GenomeRegion {
                        index1:  index1_out,
                        padded:  false,
                        chrom:   chrom.clone(),
                        start0:  region_in.start0,
                        end1:    region_in.end1,
                        names:   region_in.names.clone(),
                        centers: region_in.centers.clone(),
                    };
                    genome_regions.n_regions += 1;
                    genome_regions.sum_region_lens += size;
                    genome_regions.sum_padded_region_lens += size;
                    genome_regions.region_chroms .push(chrom.clone());
                    genome_regions.region_start0s.push(region_out.start0);
                    genome_regions.region_end1s  .push(region_out.end1);
                    genome_regions.region_names  .push(region_out.names.clone());
                    genome_regions.region_centers.push(region_out.centers.clone());
                    regions_out.push(region_out);

                // overlaps exist, collect change points and split into output regions
                } else {
                    let mut change_points0: Vec<u32> = vec![];
                    for region in overlap_group.iter() {
                        change_points0.push(region.start0); // first 0-based index where we step into   an input region
                        change_points0.push(region.end1);   // first 0-based index where we step out of an input region
                    }
                    change_points0.sort_unstable();
                    change_points0.dedup(); // 0-based indices where region membership changes

                    // create one split region for each change point interval with overlapping input regions
                    for k in 0..change_points0.len() - 1 {
                        let start0_out = change_points0[k];
                        let end1_out   = change_points0[k + 1];
                        
                        // find which regions from the group overlap this interval
                        let mut region_out_regions = vec![];
                        for region in overlap_group.iter() {
                            if region.start0 < end1_out && region.end1 > start0_out {
                                region_out_regions.push(*region);
                            }
                        }
                        index1_out += 1;
                        let size = end1_out - start0_out;
                        let region_out = GenomeRegion {
                            index1:  index1_out,
                            padded:  false,
                            chrom:   chrom.clone(),
                            start0:  start0_out,
                            end1:    end1_out,
                            names:   region_out_regions.iter().flat_map(|r| &r.names)  .cloned().collect(),
                            centers: region_out_regions.iter().flat_map(|r| &r.centers).copied().collect(),
                        };
                        genome_regions.n_regions += 1;
                        genome_regions.sum_region_lens += size;
                        genome_regions.sum_padded_region_lens += size;
                        genome_regions.region_chroms .push(chrom.clone());
                        genome_regions.region_start0s.push(region_out.start0);
                        genome_regions.region_end1s  .push(region_out.end1);
                        genome_regions.region_names  .push(region_out.names.clone());
                        genome_regions.region_centers.push(region_out.centers.clone());
                        regions_out.push(region_out);
                    }
                }

                // move to next non-overlapping region
                i = j; 
            }
        }

        // return the new GenomeRegions collection
        genome_regions
    }

    /// Get a list of unique region chromosomes as Vec<(chrom, chrom_index)>,
    /// sorted by chrom_index.
    /// 
    /// If no regions, i.e., !self.has_data, return all nuclear chromosomes.
    pub fn get_region_chroms(&self, chroms: &Chroms) -> Vec<(String, u8)> {
        let unique_chroms = if self.has_data {
            &self.unique_chroms
        } else {
            &chroms.nuclear
        };
        let mut unique_chroms = unique_chroms.iter()
            .map(|chrom| (chrom.to_string(), *chroms.index.get(chrom).unwrap()))
            .collect::<Vec<(String, u8)>>();
        unique_chroms.sort_by_key(|(_, chrom_index)| *chrom_index);
        unique_chroms
    }

    /// Find a region containing a position, if any, returning 
    /// (&GenomeRegion, Vec<pos1 - region.center>, region_match_state)
    /// for the first matching region.
    /// 
    /// Includes an assessment of padded regions if no unpadded regions match.
    /// 
    /// Returns (null_region, [0], "I") if not regions.has_data, to indicate "inside the genome".
    /// 
    /// Returns (null_region, [0], "O") if pos1 is outside any (padded) region.
    pub fn get_pos_region(&self, chrom: &str, pos1: u32) -> (&GenomeRegion, Vec<i32>, &'static str) {
        if !self.has_data {
            (&self.null_region, vec![NULL_DISTANCE], INSIDE_REGION) // always "on target/in genome" if no regions
        } else if let Some((region, distances)) = Self::get_pos_region_(chrom, pos1, &self.regions) {
            (region, distances, INSIDE_REGION)
        } else if self.region_padding == 0 {
            (&self.null_region, vec![NULL_DISTANCE], OUTSIDE_REGION)
        } else if let Some((region, distances)) = Self::get_pos_region_(chrom, pos1, &self.padded_regions) {
            (region, distances, NEAR_REGION)
        } else {
            (&self.null_region, vec![NULL_DISTANCE], OUTSIDE_REGION)
        }
    }

    /// Find a region containing a position, if any, returning
    /// (&GenomeRegion, Vec<pos1 - region.center>)
    /// for the first matching region.
    /// 
    /// Does not consider region padding; pos1 must be within an actual region.
    /// 
    /// Returns (null_region, [0]) if not regions.has_data or if pos1 is outside any region.
    pub fn get_pos_region_unpadded(&self, chrom: &str, pos1: u32) -> (&GenomeRegion, Vec<i32>) {
        if !self.has_data {
            (&self.null_region, vec![NULL_DISTANCE])
        } else if let Some((region, distances)) = Self::get_pos_region_(chrom, pos1, &self.regions) {
            (region, distances)
        } else  {
            (&self.null_region, vec![NULL_DISTANCE])
        }
    }

    // internal function for position matching
    fn get_pos_region_<'a>(
        chrom:   &str, 
        pos1:    u32, 
        regions: &'a FxHashMap<String, Vec<GenomeRegion>>
    ) -> Option<(&'a GenomeRegion, Vec<i32>)> {
        let chrom_regions = regions.get(chrom)?;
        let region = chrom_regions.iter().find(|region| { 
            pos1 >= region.start0 + 1 && pos1 <= region.end1
        })?;
        let distances = region.centers.iter()
            .map(|center| pos1 as i32 - *center as i32)
            .collect::<Vec<i32>>(); // in case region was aggregated from multiple initial regions
        Some((region, distances))
    }

    /// Report whether a position falls within any  unpadded region.
    /// 
    /// Returns `false` if not regions.has_data.
    pub fn pos_in_region(&self, chrom: &str, pos1: u32) -> bool {
        if !self.has_data {
            false
        } else if let Some(chrom_regions) = self.regions.get(chrom) {
            chrom_regions.iter().any(|region| { 
                pos1 >= region.start0 + 1 && pos1 <= region.end1
            })
        } else {
            false
        }
    }
}


/// BedGraphU16 structure for output regions
/// with integer or other coverage values.
#[derive(Serialize, Deserialize)]
pub struct BedGraphU16 {
    pub chrom_index: u8,
    pub start0:      u32,
    pub end1:        u32,
    pub coverage:    u16,
}
