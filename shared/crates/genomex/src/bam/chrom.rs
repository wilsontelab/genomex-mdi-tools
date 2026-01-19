//! Convert a BAM Record chromosome ID (tid) to a chromosome name and index.

// dependencies
use std::str::from_utf8_unchecked;
use rust_htslib::bam::{Record as BamRecord, HeaderView};
use crate::genome::Chroms;

pub fn get_chrom(
    aln:    &BamRecord, 
    header: &HeaderView, 
    chroms: &Chroms,
    allow_unmapped: bool 
) -> (String, u8){
    let tid = aln.tid();
    if tid < 0 {
        if allow_unmapped {
            return ("*".to_string(), 0);
        } else {
            panic!("BamRecord::get_chrom: alignment is unmapped (tid < 0)");
        }
    }
    let chrom = unsafe { from_utf8_unchecked(header.tid2name(tid as u32)) };
    let chrom_index = chroms.index.get(chrom).unwrap_or_else(|| {
        panic!("BamRecord::get_chrom: chrom '{chrom}' not found in Chroms")
    });
    (chrom.to_string(), *chrom_index)
}
