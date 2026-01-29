//! Suport manipulations and assessment related to genome target regions.

// dependencies
use std::ops::Deref;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use crate::sam::SamRecord;
use super::regions::*;

// constants
pub_key_constants!(
    // environment variables
    TARGETS_BED
    REGION_PADDING
    // derived variables set upstream
    IS_ONT
);

/// TargetRegions structure for tracking multiple genome target regions.
pub struct TargetRegions(pub GenomeRegions); // a struct with one unnamed field
impl Deref for TargetRegions { // call GenomeRegions methods directly on TargetRegions
    type Target = GenomeRegions;
    fn deref(&self) -> &Self::Target { &self.0 }
}
impl TargetRegions {
    /// Create a new TargetRegions collection
    /// from a BED file provided as environment variable TARGETS_BED.
    /// 
    /// See GenomeRegions::from_env() for details.
    pub fn from_env(w: &mut Workflow, has_header: bool) -> Self {
        w.cfg.set_u32_env(&[REGION_PADDING]);
        let region_padding = *w.cfg.get_u32(REGION_PADDING);
        TargetRegions(
            GenomeRegions::from_env(
                TARGETS_BED, 
                &mut w.cfg,
                has_header, 
                region_padding,
            )
        )
    }

    /// Create a new TargetRegions collection
    /// from a BED file path provided as an argument.
    /// 
    /// See GenomeRegions::new() for details.
    pub fn new(
        bed_file:         &str, 
        has_header:       bool,
        region_padding:   u32, 
    ) -> Self {
        TargetRegions(
            GenomeRegions::new(
                bed_file, 
                has_header, 
                region_padding,
                true,
            )
        )
    }

    /// Find the target region matching an alignment, if any, as
    /// (paired_state, GenomeRegion, left_pos1, right_pos1, is_on_target).
    /// 
    /// There is a general presumption that alignments are short relative to target regions,
    /// such that the two ends of an alignment are unlikely to be in different targets.
    /// If they are, only the target region matching the leftmost position is reported.
    pub fn get_aln_target(&self, aln: &SamRecord) -> (u8, &GenomeRegion, u32, u32, bool) {
        let left_pos1  = aln.pos1;
        let right_pos1 = aln.get_end1();
        if !self.has_data {
            return (self.paired_states[IN_IN_REGION], &self.null_region, left_pos1, right_pos1, true);
        }
        let (region_left,  _distance_left,  target_type_left)  = 
            self.get_pos_region(&aln.rname, left_pos1);
        let (region_right, _distance_right, target_type_right) = 
            self.get_pos_region(&aln.rname, right_pos1);
        if region_left.chrom  != NULL_REGION_NAME && 
           region_right.chrom != NULL_REGION_NAME {
            let (key, is_on_target) = if target_type_left == INSIDE_REGION { // TT, TA, or AA, but never AT
                (&format!("{}{}", target_type_left, target_type_right), true)
            } else {
                (&format!("{}{}", target_type_right, target_type_left), target_type_right == INSIDE_REGION)
            };
            (self.paired_states[key], region_left, left_pos1, right_pos1, is_on_target)
        } else if region_left.chrom != NULL_REGION_NAME {
            let key = format!("{}{}", target_type_left, OUTSIDE_REGION); // T- or A-
            (self.paired_states[&key], region_left, left_pos1, right_pos1, target_type_left == INSIDE_REGION)
        } else if region_right.chrom != NULL_REGION_NAME {
            let key = format!("{}{}", target_type_right, OUTSIDE_REGION);
            (self.paired_states[&key], region_right, left_pos1, right_pos1, target_type_right == INSIDE_REGION)
        } else {
            (self.paired_states[OUT_OUT_REGION], region_left, left_pos1, right_pos1, false) // --
        }
    }

    /// Tag each of a read's alignments with its target matching outcome.
    /// 
    /// Each tag is index1 << 3 | paired_state, where index1 is the 1-based index of the target 
    /// region in the source BED file and paired_state is an integer code (0-5) reflecting the alignment's  
    /// target match status at the two ends of the alignment (pairs of I (inside), A (adjacent), 
    /// and O (outside)).
    /// 
    /// Returns `true` if the 5' most (ONT adaptive) or any (other platforms) alignment 
    /// of the read was on target.
    pub fn set_aln_targets (
        &self,
        alns: &mut [SamRecord],
        w:    &mut Workflow,
        target_match_prefix: &'static str, // SAM tag prefixes
        off_target_prefix:   &'static str,
        n_alns_by_key:       &'static str,
        n_bases_by_key:      &'static str,
    ) -> bool {
        if !self.has_data { return true; } // always on target if not targeted
        let (paired_state_5, region_5, left_pos1_5, right_pos1_5, mut is_on_target) 
            = self.get_aln_target(&alns[0]);
        let n_alns = alns.len();
        if n_alns > 1 {
            let mut any_on_target = false;
            let is_ont = *w.cfg.get_bool(IS_ONT);
            for aln in &mut *alns {
                let (paired_state, region, left_pos1, right_pos1, aln_on_target) = self.get_aln_target(aln);
                let target_match_tag = region.index1 << 3 | paired_state as u32;
                aln.tags.tags.push(format!("{}{}", target_match_prefix, target_match_tag));
                let match_state = self.paired_state_states[paired_state as usize];
                w.ctrs.increment_keyed(n_alns_by_key,  match_state);
                w.ctrs.add_to_keyed(n_bases_by_key, match_state, (right_pos1 - left_pos1 + 1) as usize);
                any_on_target |= aln_on_target;
            }
            // ONT adaptive requires 5' alignment to be on-target
            // other platforms allow any alignment to be on-target
            if !is_ont { is_on_target = any_on_target; }
        } else {
            let target_match_tag = region_5.index1 << 3 | paired_state_5 as u32;
            alns[0].tags.tags.push(format!("{}{}", target_match_prefix, target_match_tag));
            let match_state = self.paired_state_states[paired_state_5 as usize];
            w.ctrs.increment_keyed(n_alns_by_key,  match_state);
            w.ctrs.add_to_keyed(n_bases_by_key, match_state, (right_pos1_5 - left_pos1_5 + 1) as usize);
        }
        if !is_on_target {
            for aln in alns {
                aln.tags.tags.push(format!("{}{}", off_target_prefix, 1_u8));
            }
        }
        is_on_target
    }

}
