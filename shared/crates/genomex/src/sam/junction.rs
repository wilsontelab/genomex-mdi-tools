//! SV junction methods that apply to pairs of SamRecords from the same read

// dependencies
use super::{SamRecord, flag};
use crate::sequence::rc_acgtn_str;

// enumerations
/// JunctionType describes the placement of two alignments in the genome.
pub enum JunctionType {
    Proper        = 0,
    Deletion      = 1,
    Duplication   = 2,
    Inversion     = 3,
    Translocation = 4,
}
/// JunctionMetadata decribes the properties of the junction between
/// two SamRecord alignments, including inserted, i.e., non-reference  
/// bases not present in either alignment.
/// 
/// Junction bases are in the orientation of the original read, so
/// are reverse complemented from aln1.seq when the REVERSE bit is 
/// set on aln1.flag.
pub struct JunctionMetadata {
    pub jxn_type:         JunctionType,
    pub alignment_offset: isize,
    pub jxn_bases:        String,
}

// implementation
impl SamRecord {
    /* -------------------------------------------------------------------------
    SV junction methods that apply to pairs of SamRecords from the same read
    ------------------------------------------------------------------------- */
    // Get the (jxn_type, alignment_offset, jxn_bases) from a pairs of SamRecords
    // from the same read. Read should be in 5' to 3' order on the read such that
    // aln1 is 5' of aln2.
    pub fn get_junction_metadata(aln1: &SamRecord, aln2: &SamRecord) -> JunctionMetadata {

        // determine the offset/overlap between the two alignments
        let read_len = aln1.seq.len();
        let is_reverse1 = aln1.check_flag_any(flag::REVERSE);
        let clip1 = if is_reverse1 {
            aln1.get_clip_left()
        } else {
            aln1.get_clip_right()
        };
        let clip2 = if aln2.check_flag_any(flag::REVERSE) {
            aln2.get_clip_right()
        } else {
            aln2.get_clip_left()
        };
        let alignment_offset = clip1 as isize + clip2 as isize - read_len as isize;

        // then determine the read bases corresponding to the insertion or microhomology
        // note that inserted bases must be provided by the read, not reference alignments
        // JXN_BASES leave in the orientation of the original read, thus,
        // will have two distinct rc'ed values when then same junction was sequenced in opposite orientations
        let mut jxn_bases = "*"; // a blunt joint
        if alignment_offset > 0 {
            jxn_bases = if is_reverse1 {
                &aln1.seq[clip1 - alignment_offset as usize..clip1]
            } else {
                let start = read_len - clip1;
                &aln1.seq[start..start + alignment_offset as usize]
            };
        } else if alignment_offset < 0 {
            let offset_abs = (-alignment_offset) as usize;
            jxn_bases = if is_reverse1 {
                &aln1.seq[clip1..clip1 + offset_abs]
            } else {
                let start = read_len - clip1 - offset_abs;
                &aln1.seq[start..start + offset_abs]
            };
        }
        if alignment_offset != 0 && is_reverse1 {
            JunctionMetadata{
                jxn_type: Self::get_junction_type(aln1, aln2), 
                alignment_offset,
                jxn_bases: rc_acgtn_str(jxn_bases),
            }
        } else {
            JunctionMetadata{
                jxn_type: Self::get_junction_type(aln1, aln2), 
                alignment_offset,
                jxn_bases: jxn_bases.to_string(),
            }
        }
    }

    /// Determine the JunctionType from two SamRecord alignments.
    pub fn get_junction_type(aln1: &SamRecord, aln2: &SamRecord) -> JunctionType {
        if aln1.rname != aln2.rname {
            JunctionType::Translocation
        } else if aln1.check_flag_any(flag::REVERSE) != aln2.check_flag_any(flag::REVERSE) {
            JunctionType::Inversion
        } else {
            if aln1.check_flag_any(flag::REVERSE) {
                if aln2.pos1 < aln1.pos1 {
                    return JunctionType::Deletion;
                }
            } else {
                if aln1.pos1 < aln2.pos1 {
                    return JunctionType::Deletion;
                }
            }
            JunctionType::Duplication
        }
    }
}
