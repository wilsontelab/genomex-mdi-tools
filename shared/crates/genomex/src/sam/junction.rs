//! SV junction methods that apply to pairs of SamRecords from the same read
//! to define a structural variant junction between their alignments.

// dependencies
use super::{SamRecord, flag};
use crate::sequence::rc_acgtn_str;
use crate::genome::Chroms;

// enumerations
/// The JunctionType enum describes the relative placement in the  
/// genome of two alignments ordered 5' to 3' along a sequenced read.
#[repr(u8)]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum JunctionType {
    Proper        = 0,
    Deletion      = 1,
    Duplication   = 2,
    Inversion     = 3,
    Translocation = 4,
}
/// The JunctionStrands enum describes the strand orientations 
/// of two alignments ordered 5' to 3' along a sequenced read.
#[repr(u8)]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum JunctionStrands {
    ForwardForward = 0,
    ForwardReverse = 1,
    ReverseForward = 2,
    ReverseReverse = 3,
}

/// The Junction struct describes the properties of the junction between
/// two SamRecord alignments, including inserted, i.e., non-reference  
/// bases not present in either alignment.
/// 
/// `seq` and `qual` are in the orientation of the original read, so
/// are reverse complemented from `aln5.seq` when the REVERSE bit is 
/// set on `aln5.flag`.
pub struct Junction {
    pub jxn_type:       JunctionType,
    pub strands:        JunctionStrands,
    pub node_aln5_end3: isize,  // +|- chrom_index << 9 + junction node pos1
    pub node_aln3_end5: isize,  // node_aln3_end5 = 5' end of the 3'-most alignment
    pub is_canonical:   bool,   // node5 < node3
    pub was_flipped:    bool,   // if this junction was flipped relative the original read
    pub sv_size:        usize,  // abs(abs(node3) - abs(node5)), 0 for translocations
    pub offset:         isize,  // + offset = insertion, - offset = microhomology overlap
    pub jxn_seq:        String, // in the orientation of the original read
    pub jxn_qual:       String, // qual must come last in case it has a comma, the separator used for serialization
}
impl Junction {
    /* -------------------------------------------------------------------------
    serialization methods for storing Junctions in SAM tags
    ------------------------------------------------------------------------- */
    /// Serialize a Junction into a string for storage in a SAM tag.
    pub fn serialize(&self) -> String {
        format!("{},{},{},{},{},{},{},{},{},{}",
            self.jxn_type as u8,
            self.strands as u8,
            self.node_aln5_end3,
            self.node_aln3_end5,
            self.is_canonical as u8,
            self.was_flipped as u8,
            self.sv_size,
            self.offset,
            self.jxn_seq,
            self.jxn_qual, // may have commas, so need to use splitn when deserializing
        )
    }
    /// Deserialize a Junction from a string stored in a SAM tag.
    pub fn deserialize(s: &str) -> Junction {
        let mut parts = s.splitn(10, ',');
        Junction{
            jxn_type: match parts.next().unwrap().parse::<u8>().unwrap() {
                0 => JunctionType::Proper,
                1 => JunctionType::Deletion,
                2 => JunctionType::Duplication,
                3 => JunctionType::Inversion,
                4 => JunctionType::Translocation,
                _ => panic!("SamRecord::deserialize_junction: invalid jxn_type"),
            },
            strands: match parts.next().unwrap().parse::<u8>().unwrap() {
                0 => JunctionStrands::ForwardForward,
                1 => JunctionStrands::ForwardReverse,
                2 => JunctionStrands::ReverseForward,
                3 => JunctionStrands::ReverseReverse,
                _ => panic!("SamRecord::deserialize_junction: invalid strands"),
            },
            node_aln5_end3: parts.next().unwrap().parse::<isize>().unwrap(),
            node_aln3_end5: parts.next().unwrap().parse::<isize>().unwrap(),
            is_canonical:   parts.next().unwrap().parse::<u8>().unwrap() != 0,
            was_flipped:    parts.next().unwrap().parse::<u8>().unwrap() != 0,
            sv_size:        parts.next().unwrap().parse::<usize>().unwrap(),
            offset:         parts.next().unwrap().parse::<isize>().unwrap(),
            jxn_seq:        parts.next().unwrap().to_string(),
            jxn_qual:       parts.next().unwrap().to_string(),
        }
    }
}

// implementation
impl SamRecord {
    /* -------------------------------------------------------------------------
    codifying junctions
    ------------------------------------------------------------------------- */
    // Describe a junction between two alignments from a pair of SamRecords from the
    // same read. The alignments should be in 5' to 3' order on the read such that
    // aln5 is 5' of aln3. The alignments are typically but not necessarily adjacent
    // in the query read.
    pub fn get_junction(aln5: &SamRecord, aln3: &SamRecord, chroms: &Chroms) -> Junction {

        // collect data on the read and alignments
        let read_len = aln5.seq.len();
        let is_reverse5 = aln5.check_flag_any(flag::REVERSE);
        let is_reverse3 = aln3.check_flag_any(flag::REVERSE);
        let strands = match (is_reverse5, is_reverse3) {
            (false, false) => JunctionStrands::ForwardForward,
            (false, true)  => JunctionStrands::ForwardReverse,
            (true, false)  => JunctionStrands::ReverseForward,
            (true, true)   => JunctionStrands::ReverseReverse,
        };

        // define the nodes at the innermost junction positions
        let node_aln5_end3 = Self::pack_signed_node_aln(aln5, 3, &chroms); // the 3' end of the 5' alignment
        let node_aln3_end5 = Self::pack_signed_node_aln(aln3, 5, &chroms); // the 5' end of the 3' alignment
        let node_aln5_end3_abs= node_aln5_end3.abs();
        let node_aln3_end5_abs= node_aln3_end5.abs();

        // determine the offset/overlap between the two alignments
        // + offset = insertion, - offset = microhomology overlap
        let jxn_clip5 = if is_reverse5 {
            aln5.get_clip_left()
        } else {
            aln5.get_clip_right()
        };
        let jxn_clip3 = if is_reverse3 {
            aln3.get_clip_right()
        } else {
            aln3.get_clip_left()
        };
        let offset = jxn_clip5 as isize + jxn_clip3 as isize - read_len as isize;

        // determine the read bases corresponding to the insertion or microhomology
        // inserted bases can only be provided by the read, not reference bases
        // seq leaves in the orientation of the original read and thus two reads that
        // sequenced the same junction will have seq values that are reverse complements
        let (mut jxn_seq, mut jxn_qual) = ("", ""); // a blunt joint
        if offset > 0 { // inserted bases
            let range = if is_reverse5 {
                jxn_clip5 - offset as usize .. jxn_clip5
            } else {
                let start = read_len - jxn_clip5;
                start .. start + offset as usize
            };
            (jxn_seq, jxn_qual) = (&aln5.seq[range.clone()], &aln5.qual.qual[range])
        } else if offset < 0 { // microhomology overlap
            let offset_abs = (-offset) as usize;
            let range = if is_reverse5 {
                jxn_clip5 .. jxn_clip5 + offset_abs
            } else {
                let start = read_len - jxn_clip5 - offset_abs;
                start .. start + offset_abs
            };
            (jxn_seq, jxn_qual) = (&aln5.seq[range.clone()], &aln5.qual.qual[range])
        }
        let (jxn_seq, jxn_qual) = if offset != 0 && is_reverse5 {
            (rc_acgtn_str(jxn_seq), jxn_qual.chars().rev().collect())
        } else {
            (jxn_seq.to_string(), jxn_qual.to_string())
        };

        // determine the junction type
        let jxn_type = if aln5.rname != aln3.rname {
            JunctionType::Translocation
        } else if is_reverse5 != is_reverse3 {
            JunctionType::Inversion
        } else {
            if is_reverse5 {
                if node_aln3_end5_abs < node_aln5_end3_abs {
                    JunctionType::Deletion
                } else {
                    JunctionType::Duplication
                }
            } else {
                if node_aln5_end3_abs < node_aln3_end5_abs {
                    JunctionType::Deletion
                } else {
                    JunctionType::Duplication
                }
            }
        };

        // determine the structural variant size of same-chromosome junctions
        let sv_size = match jxn_type {
            JunctionType::Translocation => 0,
            _ => (node_aln3_end5_abs - node_aln5_end3_abs).abs() as usize,
        };

        // return the junction
        Junction{
            jxn_type, 
            strands,
            node_aln5_end3,
            node_aln3_end5,
            is_canonical: node_aln5_end3 < node_aln3_end5,
            was_flipped:  false,
            sv_size,
            offset,
            jxn_seq,
            jxn_qual,
        }
    }
    /* -------------------------------------------------------------------------
    codifying the breakpoint nodes connected by junctions
    ------------------------------------------------------------------------- */
    /// Get a signed 64-bit (isize)-encoded node at one end of a SamRecord alignment.
    /// Value is +|- chrom_index << 9 + position (1-based).
    pub fn pack_signed_node_aln(&self, end: u8, chroms: &Chroms) -> isize {
        let is_reverse = self.check_flag_any(flag::REVERSE);
        let pos1 = match (end, is_reverse) {
            (5, true)  => self.get_end1(),
            (5, false) => self.pos1,
            (3, true)  => self.pos1,
            (3, false) => self.get_end1(),
            _ => panic!("SamRecord::get_signed_node: `end` must be 5 or 3"),
        };
        Self::pack_signed_node(&self.rname, pos1, is_reverse, chroms)
    }
    /// Get a signed 64-bit (isize)-encoded node based on a SamRecord alignment's
    /// chrom and strand and a caller-provided position (1-based) and strand inversion.
    /// 
    /// This function is useful for calculating a node at a position other than the
    /// alignment end itself, e.g., at a projected site position.
    pub fn pack_signed_node_at_pos(&self, pos1: usize, invert_strand: bool, chroms: &Chroms) -> isize {
        let rev_bool = if invert_strand { false } else { true };
        let is_reverse = if self.check_flag_any(flag::REVERSE) { rev_bool } else { !rev_bool };
        Self::pack_signed_node(&self.rname, pos1, is_reverse, chroms)
    }
    /// Encode a chrom, pos1, and strand into a signed 64-bit (isize) node.
    /// Value is +|- chrom_index << 9 + position (1-based).
    pub fn pack_signed_node(chrom: &str, pos1: usize, is_reverse: bool, chroms: &Chroms) -> isize {
        let chrom_index_shift9 = (*chroms.index.get(chrom).or(Some(&0)).unwrap() << 9) as isize;
        match is_reverse{
            true  => (chrom_index_shift9 + pos1 as isize) * -1,
            false =>  chrom_index_shift9 + pos1 as isize,
        }
    }
    /// Extract the chromosome, position and strand (is_reverse) from a signed node
    /// provided as isize.
    pub fn unpack_signed_node(node: isize, chroms: &Chroms) -> (String, usize, bool) {
        let chrom_index = (node.abs() >> 9) as usize;
        let pos1 = (node.abs() & 0x1FF) as usize;
        let is_reverse = node < 0;
        let chrom = chroms.rev_index.get(&chrom_index).unwrap().to_string();
        (chrom, pos1, is_reverse)
    }
    /// Extract the chromosome, position and strand from a signed node
    /// provided as &str.
    pub fn unpack_signed_node_str(node: &str, chroms: &Chroms) -> (String, usize, bool) {
        let node = node.parse::<isize>().unwrap();
        Self::unpack_signed_node(node, chroms)
    }
}
