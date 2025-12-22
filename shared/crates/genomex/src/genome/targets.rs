//! Suport manipulations and assessment related to genome target regions.

// dependencies
use std::collections::HashMap;
use std::fs::read_to_string;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use crate::sam::SamRecord;

// constants
pub_key_constants!(
    // environment variables
    TARGETS_BED
    REGION_PADDING
    // derived variables set upstream
    IS_ONT
    // counters
    N_ALNS_BY_TARGET_CLASS
    N_BASES_BY_TARGET_CLASS
);
// target states for single positions and alignments
pub const ON_TARGET: &str        = "T";
pub const NEAR_TARGET: &str      = "A"; // for "adjacent"
pub const OFF_TARGET: &str       = "-";
// target states for alignment pairs
pub const ON_ON_TARGET: &str     = "TT";
pub const OFF_OFF_TARGET: &str   = "--";
// null target constants
pub const NULL_TARGET_I: usize   = 0;
pub const NULL_TARGET_NAME: &str = "*";
pub const NULL_DISTANCE: isize   = 0;

/// TargetRegion structure to describe one target region.
#[derive(Clone)]
pub struct TargetRegion {
    pub chrom:       String,
    pub start0:      usize, // 0-based inclusive
    pub end1:        usize, // 1-based inclusive (BED convention)
    pub name:        String,
    pub target_i1:   usize, // 1-referenced index in source BED file
    pub center:      usize,
    pub padded:      bool,
}

/// TargetRegions structure for tracking multiple genome target regions.
pub struct TargetRegions {
    is_targeted:            bool,
    targets_bed:            String,
    region_padding:         usize,
    n_regions:              usize,
    sum_target_lens:        usize,
    sum_padded_target_lens: usize,
    region_centers:         Vec<usize>,
    target_chroms:          Vec<String>,
    regions:                HashMap<String, Vec<TargetRegion>>, // chromosome -> sorted regions
    padded_regions:         HashMap<String, Vec<TargetRegion>>,
    target_classes:         HashMap<String, usize>,
    count_classes:          Vec<&'static str>,
    null_target_region:     TargetRegion,   
}
impl TargetRegions {
    /// Create a new TargetRegions instance.
    pub fn new(w: &mut Workflow) -> Self {

        // load config variables
        w.cfg.set_bool_env(&[TARGETS_BED]);
        w.cfg.set_string_env(&[TARGETS_BED]);
        w.cfg.set_usize_env(&[REGION_PADDING]);
        let is_targeted = *w.cfg.get_bool(TARGETS_BED);
        let region_padding = *w.cfg.get_usize(REGION_PADDING);
        let target_classes: HashMap<String, usize> = [
            ("--".to_string(), 0),
            ("TT".to_string(), 1),
            ("TA".to_string(), 2),
            ("T-".to_string(), 3),
            ("AA".to_string(), 4),
            ("A-".to_string(), 5)
        ].into_iter().collect();
        let count_classes = vec![
            OFF_TARGET,
            ON_TARGET,
            ON_TARGET,
            ON_TARGET,
            NEAR_TARGET,
            NEAR_TARGET
        ];

        // initialize counters
        w.ctrs.add_keyed_counters(&[
            (N_ALNS_BY_TARGET_CLASS,  "alignment counts by target class"),
            (N_BASES_BY_TARGET_CLASS, "aligned base counts by target class")
        ]);

        // initialize TargetRegions instance
        let mut targets = TargetRegions {
            is_targeted,
            targets_bed: w.cfg.get_string(TARGETS_BED).to_string(),
            region_padding,
            n_regions: 0,
            sum_target_lens: 0,
            sum_padded_target_lens: 0,
            region_centers: vec![],
            target_chroms: vec![],
            regions: HashMap::new(),
            padded_regions: HashMap::new(),
            target_classes,
            count_classes,
            null_target_region: TargetRegion{
                chrom:       NULL_TARGET_NAME.to_string(),
                start0:      0,
                end1:        0,
                name:        NULL_TARGET_NAME.to_string(),
                target_i1:   NULL_TARGET_I,
                center:      0,
                padded:      false,
            }
        };  

        // load target regions, if any
        if is_targeted {
            // first pass to load target regions (T=target type)
            targets.load_target_regions(0);

            // second pass to load padded regions (A=adjacent type)
            if region_padding > 0 {
                targets.load_target_regions(region_padding);
            } else {
                targets.sum_padded_target_lens = targets.sum_target_lens;
            }
        }

        // return the Targets instance
        targets
    }
    // load the target regions from the BED file
    fn load_target_regions(&mut self, padding: usize) {
        let bed_content = read_to_string(&self.targets_bed).expect(
            &format!("could not open {}: ", self.targets_bed)
        );  
        let mut target_i1 = 0;
        let mut sum_target_lens = 0;
        for line in bed_content.lines() {
            target_i1 += 1; // 1-referenced
            let line = line.replace("\r", "");
            let parts: Vec<&str> = line.split('\t').collect();
            let chrom = parts[0];
            let mut start0: usize = parts[1].parse().unwrap_or(0);
            let mut end1: usize = parts[2].parse().unwrap_or(0);
            let name = if parts.len() >= 4 { parts[3] } else { "unnamed" };
            sum_target_lens += end1 - start0;
            let center = ((start0 as f64 + 1.0 + end1 as f64) / 2.0) as usize;
            if padding == 0 {
                self.region_centers.push(center);
                self.n_regions += 1;
                if !self.target_chroms.contains(&chrom.to_string()) {
                    self.target_chroms.push(chrom.to_string());
                }
            } else {
                start0 = start0.saturating_sub(padding);
                end1 += padding;
            }
            let region = TargetRegion {
                chrom: chrom.to_string(),
                start0,
                end1,
                name: name.to_string(),
                target_i1,
                center,
                padded: padding > 0,
            };
            if padding == 0 {
                self.regions.entry(chrom.to_string()).or_insert_with(Vec::new).push(region);
            } else {
                self.padded_regions.entry(chrom.to_string()).or_insert_with(Vec::new).push(region);
            }
        }
        if padding == 0 {
            for regions in self.regions.values_mut() {
                regions.sort_by_key(|r| r.start0);
            }
            self.sum_target_lens = sum_target_lens;
        } else {
            for regions in self.padded_regions.values_mut() {
                regions.sort_by_key(|r| r.start0);
            }
            self.sum_padded_target_lens = sum_target_lens;
        }
    }


//     /// Print a summary of the target region count and sizes to the log.
//     pub fn print_targets_summary(&self, genome_size: usize) {
// //         printCount($nRegions,            'nRegions',            'target regions');
// //         printCount(commify($sumTargetLens),       'sumTargetLens',       'total bp covered by target regions');
// //         $REGION_PADDING and 
// //         printCount(commify($sumPaddedTargetLens), 'sumPaddedTargetLens', 'total bp covered by padded target regions');
//     }

    /// Find the target region containing a position, if any.
    pub fn get_pos_target(&self, chrom: &str, pos1: usize) -> (TargetRegion, isize, String) {
        if !self.is_targeted {
            (self.null_target_region.clone(), NULL_DISTANCE, ON_TARGET.to_string()) // always on target if not targeted
        } else if let Some((region, distance)) = self.get_pos_target_(chrom, pos1, &self.regions) {
            (region, distance, ON_TARGET.to_string())
        } else if self.region_padding == 0 {
            (self.null_target_region.clone(), NULL_DISTANCE, OFF_TARGET.to_string())
        } else if let Some((region, distance)) = self.get_pos_target_(chrom, pos1, &self.padded_regions) {
            (region, distance, NEAR_TARGET.to_string())
        } else {
            (self.null_target_region.clone(), NULL_DISTANCE, OFF_TARGET.to_string())
        }
    }
    // internal function for position matching
    fn get_pos_target_(&self, chrom: &str, pos1: usize, regions: &HashMap<String, Vec<TargetRegion>>) -> Option<(TargetRegion, isize)> {
        let regions = regions.get(chrom)?;
        let region = regions.iter().find(|region| { 
            pos1 >= region.start0 + 1 && pos1 <= region.end1
        })?;
        Some((region.clone(), pos1 as isize - region.center as isize))
    }

    /// Find the target region matching an alignment, if any.
    /// 
    /// There is a general presumption that alignments are short relative to target regions,
    /// such that the two ends of an alignment are unlikely to be in different targets.
    /// If they are, only the target region matching the leftmost position is reported
    pub fn get_aln_target(&self, aln: &SamRecord) -> (usize, TargetRegion, usize, usize, bool) {
        let left_pos1 = aln.pos1;
        let right_pos1 = aln.get_end1();
        if !self.is_targeted {
            return (self.target_classes[ON_ON_TARGET], self.null_target_region.clone(), left_pos1, right_pos1, true);
        }
        let (region_left, _distance_left, target_type_left) = self.get_pos_target(&aln.rname, left_pos1);
        let (region_right, _distance_right, target_type_right) = self.get_pos_target(&aln.rname, right_pos1);
        if region_left.chrom  != NULL_TARGET_NAME && 
           region_right.chrom != NULL_TARGET_NAME {
            let (key, is_on_target) = if target_type_left == ON_TARGET { // TT, TA, or AA, but never AT
                (&format!("{}{}", target_type_left, target_type_right), true)
            } else {
                (&format!("{}{}", target_type_right, target_type_left), target_type_right == ON_TARGET)
            };
            (self.target_classes[key], region_left, left_pos1, right_pos1, is_on_target)
        } else if region_left.chrom != NULL_TARGET_NAME {
            let key = format!("{}{}", target_type_left, OFF_TARGET); // T- or A-
            (self.target_classes[&key], region_left, left_pos1, right_pos1, target_type_left == ON_TARGET)
        } else if region_right.chrom != NULL_TARGET_NAME {
            let key = format!("{}{}", target_type_right, OFF_TARGET);
            (self.target_classes[&key], region_right, left_pos1, right_pos1, target_type_right == ON_TARGET)
        } else {
            (self.target_classes[OFF_OFF_TARGET], region_left, left_pos1, right_pos1, false) // --
        }
    }

    /// Tag a read's alignments with a target class that reflects the target status of 
    /// the outermost alignments alongside each alignments target status.
    /// 
    /// Return `true` if the 5' most (ONT adaptive) or any (other platforms) alignment 
    /// of the read was on target
    pub fn set_aln_target_classes (
        &self,
        alns: &mut [SamRecord],
        w:    &mut Workflow,
        class_prefix:      &'static str,
        off_target_prefix: &'static str,
    ) -> bool {
        if !self.is_targeted { return true; } // always on target if not targeted
        let (target_class_5, region_5, left_pos1_5, right_pos1_5, mut is_on_target) 
            = self.get_aln_target(&alns[0]);
        let n_alns = alns.len();
        if n_alns > 1 {
            let mut any_on_target = false;
            for aln in &mut *alns {
                let (target_class, region, left_pos1, right_pos1, aln_on_target) = self.get_aln_target(aln);
                let tag = region.target_i1 << 3 | target_class;
                aln.tags.tags.push(format!("{}{}", class_prefix, tag));
                let count_class = self.count_classes[target_class];
                w.ctrs.increment_keyed(N_ALNS_BY_TARGET_CLASS, count_class);
                w.ctrs.add_to_keyed(N_BASES_BY_TARGET_CLASS, count_class, right_pos1 - left_pos1 + 1);
                any_on_target |= aln_on_target;
            }
            // ONT adaptive requires 5' alignment to be on-target
            // other platforms allow any alignment to be on-target
            if !*w.cfg.get_bool(IS_ONT) {
                is_on_target = any_on_target;
            }
        } else {
            let tag: usize = region_5.target_i1 << 3 | target_class_5;
            alns[0].tags.tags.push(format!("{}{}", class_prefix, tag));
            let count_class = self.count_classes[target_class_5];
            w.ctrs.increment_keyed(N_ALNS_BY_TARGET_CLASS, count_class);
            w.ctrs.add_to_keyed(N_BASES_BY_TARGET_CLASS, count_class, right_pos1_5 - left_pos1_5 + 1);
        }
        if !is_on_target {
            for aln in alns {
                aln.tags.tags.push(format!("{}{}", off_target_prefix, 1_u8));
            }
        }
        is_on_target
    }

    // sub getTargetRegions {
    //     ($TARGETS_BED and $TARGETS_BED ne "null" and $TARGETS_BED ne "NA") or return [];
    //     my @regions;
    //     open my $inH, "<", $TARGETS_BED or die "could not open $TARGETS_BED: $!\n";
    //     while(my $line = <$inH>){
    //         chomp $line;
    //         $line =~ s/\r//g;
    //         my ($chr, $start, $end, $name) = split("\t", $line);
    //         push @regions, {
    //             name  => $name,
    //             chr   => $chr,
    //             start => $start,
    //             end   => $end,
    //             paddedStart  => $start - $REGION_PADDING, # all still half-open like the source BED
    //             paddedEnd    => $end   + $REGION_PADDING,
    //             paddedStart1 => $start - $REGION_PADDING + 1 # 1-referenced start, unlike above
    //         };
    //     }
    //     close $inH;
    //     return \@regions;
    // }

}
