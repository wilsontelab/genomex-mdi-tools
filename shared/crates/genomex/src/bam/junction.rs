//! SV junction methods that apply to pairs of BAM Records from the same read
//! to define a structural variant junction between their alignments.

// dependencies
use rust_htslib::bam::{Record as BamRecord, HeaderView};
use crate::genome::Chroms;
use crate::sam::SamRecord;
use super::chrom::get_chrom;

/* -------------------------------------------------------------------------
codifying the breakpoint nodes connected by junctions; also applies to read outer nodes
------------------------------------------------------------------------- */
/// Get a signed 64-bit (isize)-encoded node at one end of a BamRecord alignment.
/// Value is +|- chrom_index << 29 + position (1-based).
pub fn pack_signed_node_aln(
    aln:    &BamRecord, 
    header: &HeaderView, 
    end:    u8, 
    chroms: &Chroms
) -> isize {
    let is_reverse = aln.is_reverse();
    let pos1 = match (end, is_reverse) {
        (5, true)  => aln.cigar().end_pos(),
        (5, false) => aln.pos() + 1,
        (3, true)  => aln.pos() + 1,
        (3, false) => aln.cigar().end_pos(),
        _ => panic!("BamRecord::get_signed_node: `end` must be 5 or 3"),
    };
    let chrom = get_chrom(aln, header, chroms, false);
    SamRecord::pack_signed_node(&chrom.0, pos1 as u32, is_reverse, chroms)
}
/// Get a signed 64-bit (isize)-encoded node based on a BamRecord alignment's
/// chrom and strand and a caller-provided position (1-based) and strand inversion.
/// 
/// This function is useful for calculating a node at a position other than the
/// alignment end itself, e.g., at a projected site position.
pub fn pack_signed_node_at_pos(
    aln:           &BamRecord, 
    header:        &HeaderView, 
    pos1:          u32, 
    invert_strand: bool, 
    chroms:        &Chroms
) -> isize {
    let rev_bool = if invert_strand { false } else { true };
    let is_reverse = if aln.is_reverse() { rev_bool } else { !rev_bool };
    let chrom = get_chrom(aln, header, chroms, false);
    SamRecord::pack_signed_node(&chrom.0, pos1, is_reverse, chroms)
}
/* -------------------------------------------------------------------------
extract the breakpoint nodes connected by junctions
------------------------------------------------------------------------- */
/// Extract one end's chrom, chrom_index, position and strand (is_reverse) 
/// from a BamRecord alignment.
pub fn get_unpacked_node_aln(
    aln:    &BamRecord, 
    header: &HeaderView, 
    end:    u8, 
    chroms: &Chroms
) -> (String, u8, u32, bool) {
    let is_reverse = aln.is_reverse();
    let pos1 = match (end, is_reverse) {
        (5, true)  => aln.cigar().end_pos(),
        (5, false) => aln.pos() + 1,
        (3, true)  => aln.pos() + 1,
        (3, false) => aln.cigar().end_pos(),
        _ => panic!("BamRecord::get_signed_node: `end` must be 5 or 3"),
    };
    let chrom = get_chrom(aln, header, chroms, false);
    (chrom.0, chrom.1, pos1 as u32, is_reverse)
}
