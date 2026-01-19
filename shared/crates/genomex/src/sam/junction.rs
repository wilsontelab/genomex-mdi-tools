//! SV junction methods that apply to pairs of SamRecords from the same read
//! to define a structural variant junction between their alignments.
//
//  logic of junction node ordering such that canonical == node1.abs() <= node2.abs()
//    translocations sort low chromI to high chromI, i.e., by refPos1 along the composite genome
//    same-chrom junctions sort low refPos1 to high refPos1
//    as a consequence, the canonical ordering of junctions is:
//      deletions     sort to +/+, i.e., 1--->~~~2---> as canonical, <---2~~~<---1 gets flipped to it
//                                           *   *                       *   *
//      duplications  sort to -/-, i.e., <---1~~~<---2 as canonical, 2--->~~~1---> gets flipped to it
//                                       *           *               *           *
//      FF inversions sort to +/-, i.e., 1--->~~~<---2 as canonical, <---2~~~1---> gets flipped to it
//                                           *       *                   *       *
//      RR inversions sort to -/+, i.e., <---1~~~2---> as canonical, 2--->~~~<---1 gets flipped to it
//                                       *       *                   *       *
//    the above patterns apply the same to translocations, where ~~|~~ simply crosses (a) chromosome boundary(ies)
//
//  similiary, the logic of outer node ordering such that canonical == node1.abs() <= node2.abs()
//      nonSV/del     sort to +/+, i.e., 1--->~~~2---> as canonical, <---2~~~<---1 gets flipped to it
//                                       *           *               *           *
//      most dups     sort to -/-, i.e., <---1~~~<---2 as canonical, 2--->~~~1---> gets flipped to it
//                                           *   *                       *   *
//      FF inversions sort to +/-, i.e., 1--->~~~<---2 as canonical, <---2~~~1---> gets flipped to it
//                                       *       *                   *       *
//      RR inversions sort to -/+, i.e., <---1~~~2---> as canonical, 2--->~~~<---1 gets flipped to it
//                                           *       *                   *       *

// dependencies
use crate::sequence::rc_acgtn_str;
use crate::genome::Chroms;
use super::{SamRecord, flag};

/// JunctionType describes the relative placement in the  
/// genome of two alignments ordered 5' to 3' along a sequenced read.
/// 
/// Use bitwise encoding since a read may need to represent >1 junctions,
/// e.g., one at each end of an alignment.
#[repr(u8)]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum JunctionType {
    Proper        = 0,
    Deletion      = 1,
    Duplication   = 2,
    Inversion     = 4,
    Translocation = 8,
}

/// JunctionStrands describe the strand orientations 
/// of two alignments ordered 5' to 3' along a sequenced read.
/// 
/// JunctionStrands is required in addition to JunctionType
/// to fully describe inversion and translocation junctions.
#[repr(u8)]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum JunctionStrands {
    ForwardForward = 0,
    ForwardReverse = 1,
    ReverseForward = 2,
    ReverseReverse = 3,
}

/// OrderedNodes holds two nodes in a canonical order.
/// 
/// The two nodes might flank an SV junction or represent
/// the outer ends of a sequenced insert.
#[derive(Clone, Copy)]
pub struct OrderedNodes {
    pub node1:         isize,
    pub node2:         isize,
    pub was_reordered: bool,
}

/// Junction holds the properties of the junction between
/// two SamRecord alignments, including inserted, i.e., non-reference  
/// bases not present in either alignment.
/// 
/// `seq` and `qual` are in the orientation of the original read, so
/// are reverse complemented from `aln5.seq` when the REVERSE bit is 
/// set on `aln5.flag`.
pub struct Junction {
    /* ---------------------------------------------------------------------- */
    // two nodes from flanking alignments define a junction
    pub node_aln5_end3: isize,  // +|- chrom_index << 29 + junction node pos1
    pub node_aln3_end5: isize,  // node_aln3_end5 = 5' end of the 3'-most alignment
    /* ---------------------------------------------------------------------- */
    // node pairs come in two orientations, only one of which is the canonical comparison orientation
    pub is_canonical:   bool,   // node_aln5_end3.abs() <= node_aln3_end5.abs()
    pub was_flipped:    bool,   // if this junction was flipped relative the original read to force is_canonical=true
    /* ----------------------------------------------------------------------
    the following properties derive deterministically from the two defining nodes
    ------------------------------------------------------------------------ */
    // JunctionType and JunctionStrands together describe the junction type
    pub jxn_type:       JunctionType,
    pub strands:        JunctionStrands,
    // SV size applies only to same-chromosome junctions only
    pub sv_size:        u32,  // abs(abs(node3) - abs(node5)), 0 for translocations
    /* ----------------------------------------------------------------------
    the following properties are specific to one instance of the junction based in SEQ and QUAL
    ------------------------------------------------------------------------ */
    pub offset:         i16,    // + offset = insertion, - offset = microhomology overlap
    pub jxn_seq:        String, // in the orientation of the original read, or its reverse complement if was_flipped=true
    pub jxn_qual:       String, // qual must come last in case it has a comma, the separator used for serialization
}
impl Junction {
    /* -------------------------------------------------------------------------
    serialization methods for storing Junctions in SAM tags
    ------------------------------------------------------------------------- */
    /// Serialize a Junction into a string for storage in a SAM tag.
    pub fn serialize(&self) -> String {
        format!("{},{},{},{},{},{},{},{},{},{}",
            self.node_aln5_end3,
            self.node_aln3_end5,
            self.is_canonical as u8,
            self.was_flipped as u8,
            self.jxn_type as u8,
            self.strands as u8,
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
            node_aln5_end3: parts.next().unwrap().parse::<isize>().unwrap(),
            node_aln3_end5: parts.next().unwrap().parse::<isize>().unwrap(),
            is_canonical:   parts.next().unwrap().parse::<u8>().unwrap() != 0,
            was_flipped:    parts.next().unwrap().parse::<u8>().unwrap() != 0,
            jxn_type: match parts.next().unwrap().parse::<u8>().unwrap() {
                0 => JunctionType::Proper,
                1 => JunctionType::Deletion,
                2 => JunctionType::Duplication,
                4 => JunctionType::Inversion,
                8 => JunctionType::Translocation,
                _ => panic!("SamRecord::deserialize_junction: invalid jxn_type"),
            },
            strands: match parts.next().unwrap().parse::<u8>().unwrap() {
                0 => JunctionStrands::ForwardForward,
                1 => JunctionStrands::ForwardReverse,
                2 => JunctionStrands::ReverseForward,
                3 => JunctionStrands::ReverseReverse,
                _ => panic!("SamRecord::deserialize_junction: invalid strands"),
            },
            sv_size:        parts.next().unwrap().parse::<u32>().unwrap(),
            offset:         parts.next().unwrap().parse::<i16>().unwrap(),
            jxn_seq:        parts.next().unwrap().to_string(),
            jxn_qual:       parts.next().unwrap().to_string(),
        }
    }
    /* -------------------------------------------------------------------------
    junction manipulation
    ------------------------------------------------------------------------- */
    /// Flip a junction so that its nodes are in canonical order, to support
    /// comparison of junctions independently of their read orientation.
    pub fn orient_junction(mut self) -> Self {
        let ordered_nodes = SamRecord::order_paired_nodes(self.node_aln5_end3, self.node_aln3_end5);
        self.is_canonical = true;
        if ordered_nodes.was_reordered {
            self.node_aln5_end3 = ordered_nodes.node1;
            self.node_aln3_end5 = ordered_nodes.node2;
            self.was_flipped  = true;
            self.strands = match self.strands {
                JunctionStrands::ForwardForward => JunctionStrands::ReverseReverse,
                JunctionStrands::ForwardReverse => JunctionStrands::ReverseForward,
                JunctionStrands::ReverseForward => JunctionStrands::ForwardReverse,
                JunctionStrands::ReverseReverse => JunctionStrands::ForwardForward,
            };
            self.jxn_seq = rc_acgtn_str(&self.jxn_seq);
            self.jxn_qual = self.jxn_qual.chars().rev().collect();
        }
        self
    }
}

// implementation
impl SamRecord {
    /* -------------------------------------------------------------------------
    codifying junctions
    ------------------------------------------------------------------------- */
    /// Describe a junction between two alignments from a pair of SamRecords from the
    /// same read. The alignments should be in 5' to 3' order on the read such that
    /// aln5 is 5' of aln3. The alignments are typically but not necessarily adjacent
    /// in the query read. No re-orientation of the nodes is performed yet.
    pub fn get_junction(aln5: &SamRecord, aln3: &SamRecord, chroms: &Chroms) -> Junction {

        // collect data on the read and alignments
        let read_len     = aln5.seq.len();
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
            aln5.get_clip_left() as usize
        } else {
            aln5.get_clip_right() as usize
        };
        let jxn_clip3 = if is_reverse3 {
            aln3.get_clip_right() as usize
        } else {
            aln3.get_clip_left() as usize
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
            _ => (node_aln3_end5_abs - node_aln5_end3_abs).abs() as u32,
        };

        // return the junction
        Junction{
            node_aln5_end3,
            node_aln3_end5,
            is_canonical: node_aln5_end3_abs <= node_aln3_end5_abs,
            was_flipped:  false,
            jxn_type, 
            strands,
            sv_size,
            offset: offset as i16,
            jxn_seq,
            jxn_qual,
        }
    }
    /* -------------------------------------------------------------------------
    codifying the breakpoint nodes connected by junctions; also applies to read outer nodes
    ------------------------------------------------------------------------- */
    /// Get a signed 64-bit (isize)-encoded node at one end of a SamRecord alignment.
    /// Value is +|- chrom_index << 29 + position (1-based).
    pub fn pack_signed_node_aln(
        &self, 
        end:    u8, 
        chroms: &Chroms
    ) -> isize {
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
    pub fn pack_signed_node_at_pos(
        &self, 
        pos1:          u32, 
        invert_strand: bool, 
        chroms:        &Chroms
    ) -> isize {
        let rev_bool = if invert_strand { false } else { true };
        let is_reverse = if self.check_flag_any(flag::REVERSE) { rev_bool } else { !rev_bool };
        Self::pack_signed_node(&self.rname, pos1, is_reverse, chroms)
    }
    /// Encode a chrom, pos1, and strand into a signed 64-bit (isize) node.
    /// Value is +|- chrom_index << 29 + position (1-based).
    pub fn pack_signed_node(
        chrom:      &str, 
        pos1:       u32, 
        is_reverse: bool, 
        chroms:     &Chroms
    ) -> isize {
        // 29 bits for position allows up to 536,870,911 1-based positions per chromosome
        let chrom_index = *chroms.index.get(chrom).or(Some(&0)).unwrap();
        let chrom_index_shift29 = (chrom_index << 29) as isize;
        let pos1 = pos1 & 0x1FFFFFFF;
        match is_reverse{
            true  => (chrom_index_shift29 + pos1 as isize) * -1,
            false =>  chrom_index_shift29 + pos1 as isize,
        }
    }
    /// Extract the chromosome, chrom_index, position and strand (is_reverse) from a signed node
    /// provided as isize.
    pub fn unpack_signed_node(node: isize, chroms: &Chroms) -> (String, u8, u32, bool) {
        let chrom_index = (node.abs() >> 29) as u8;
        let pos1 = (node.abs() & 0x1FFFFFFF) as u32;
        let is_reverse = node < 0;
        let chrom = chroms.rev_index.get(&chrom_index).unwrap().to_string();
        (chrom, chrom_index, pos1, is_reverse)
    }
    /// Extract one end's chrom, chrom_index, position and strand (is_reverse) 
    /// from a SamRecord alignment.
    pub fn get_unpacked_node_aln(
        &self, 
        end:    u8, 
        chroms: &Chroms
    ) -> (String, u8, u32, bool) {
        let chrom_index = *chroms.index.get(&self.rname).or(Some(&0)).unwrap();
        let is_reverse = self.check_flag_any(flag::REVERSE);
        let pos1 = match (end, is_reverse) {
            (5, true)  => self.get_end1(),
            (5, false) => self.pos1,
            (3, true)  => self.pos1,
            (3, false) => self.get_end1(),
            _ => panic!("SamRecord::get_unpacked_node_aln: `end` must be 5 or 3"),
        };
        (self.rname.clone(), chrom_index, pos1, is_reverse)
    }
    /* -------------------------------------------------------------------------
    codifying node pairs as SAM tags, especially outer nodes
    ------------------------------------------------------------------------- */
    /// Encode two nodes as a SAM tag string, in the order provided.
    pub fn paired_nodes_to_sam_tag(node1: isize, node2: isize) -> String {
        format!("{},{}", node1, node2)
    }
    /// Decode two nodes from a SAM tag string value encoded by `paired_nodes_to_sam_tag()`.
    /// Optionally sort the nodes using `order_paired_nodes()`.
    pub fn sam_tag_to_paired_nodes(tag: &str, sort: bool) -> OrderedNodes {
        let parts: Vec<&str> = tag.split(',').collect();
        let node1: isize = parts[0].parse().unwrap();
        let node2: isize = parts[1].parse().unwrap();
        if sort {
            Self::order_paired_nodes(node1, node2)
        } else {
            OrderedNodes { node1, node2, was_reordered: false }
        }
    }
    //  1            9  => +1,+9 (1 < 9, no reorder) ==> 1,9
    //  ----->  ----->
    //  <-----  <-----
    // -1           -9  => -9,-1 (9 > 1, reorder)    ==> 1,9
    /// Order two nodes such that first.abs() is less than or equal to second.abs().
    /// If nodes are reordered, return was_reordered=true and invert their signs.
    pub fn order_paired_nodes(node1: isize, node2: isize) -> OrderedNodes {
        if node1.abs() <= node2.abs() {
            OrderedNodes { node1, node2, was_reordered: false }
        } else {
            OrderedNodes { node1: -node2, node2: -node1, was_reordered: true }
        }
    }
}
