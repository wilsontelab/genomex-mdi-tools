//! Support functions for manipulating unaligned BAM Record Aux tags.

// dependencies
use rust_htslib::bam::{Record as BamRecord, record::Aux};

/* ---------------------------------------------------------------------
type-generic support methods
--------------------------------------------------------------------- */
/// Get the Aux value of a specified tag from a BAM record.
/// Panic if tag is absent.
pub fn get_tag_value<'a>(aln: &'a BamRecord, tag: &'static str) -> Aux<'a> {
    aln.aux(tag.as_bytes()).unwrap_or_else(|e|{
        panic!("Failed to extract {tag} tag: {}", e)
    })
}
/// Get the Aux value of a specified tag from a BAM record.
/// Return None if the tag is not present.
pub fn get_tag_value_opt<'a>(aln: &'a BamRecord, tag: &'static str) -> Option<Aux<'a>> {
    aln.aux(tag.as_bytes()).ok()
}
/* ---------------------------------------------------------------------
type-specific methods that panic on missing tags
--------------------------------------------------------------------- */
/// Parse the value of a specified 'XX:i:' tag from a BAM record as u8.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_u8(aln: &BamRecord, tag: &'static str) -> u8 {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::I32(v) => v as u8,
        Aux::U8(v)   => v,
        Aux::I16(v) => v as u8,
        Aux::U32(v) => v as u8,
        Aux::I8(v)   => v as u8,
        Aux::U16(v) => v as u8,
        _ => panic!("{tag} tag is not u8: {:?}", aux),
    }
}
/// Parse the value of a specified 'XX:i:' tag from a BAM record as u32.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_u32(aln: &BamRecord, tag: &'static str) -> u32 {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::I32(v) => v as u32,
        Aux::U32(v) => v,
        Aux::I16(v) => v as u32,
        Aux::U16(v) => v as u32,
        Aux::U8(v)   => v as u32,
        Aux::I8(v)   => v as u32,
        _ => panic!("{tag} tag is not u32: {:?}", aux),
    }
}
/// Parse the value of a specified 'XX:i:' tag from a BAM record as i32.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_i32(aln: &BamRecord, tag: &'static str) -> i32 {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::I32(v) => v,
        Aux::U32(v) => v as i32,
        Aux::I16(v) => v as i32,
        Aux::U16(v) => v as i32,
        Aux::U8(v)   => v as i32,
        Aux::I8(v)   => v as i32,
        _ => panic!("{tag} tag is not i32: {:?}", aux),
    }
}
/// Parse the value of a specified 'XX:f:' tag from a BAM record as f32.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_f32(aln: &BamRecord, tag: &'static str) -> f32 {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::Float(v) => v,
        Aux::Double(v) => v as f32,
        _ => panic!("{tag} tag is not f32: {:?}", aux),
    }
}
/// Parse the value of a specified 'XX:A:' tag from a BAM record as bool.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_bool(aln: &BamRecord, tag: &'static str) -> bool {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::Char(v) => v != 0,
        Aux::U8(v)   => v != 0,
        _ => panic!("{tag} tag is not a bool-compatible type: {:?}", aux),
    }
}
/// Parse the value of a specified 'XX:Z:' tag from a BAM record as String.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_str(aln: &BamRecord, tag: &'static str) -> String {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::String(v) => v.to_string(),
        _ => panic!("{tag} tag is not a string: {:?}", aux),
    }
}
/* ---------------------------------------------------------------------
type-specific methods that return defaults on missing tags
--------------------------------------------------------------------- */
/// Parse the value of a specified'XX:i:' tag from a BAM record as u8.
/// Return `default` if the tag is not present.
pub fn get_tag_u8_default(aln: &BamRecord, tag: &'static str, default: u8) -> u8 {
    let aux = get_tag_value_opt(aln, tag);
    match aux {
        None => default,
        Some(Aux::I32(v)) => v as u8,
        Some(Aux::U8(v))   => v,
        Some(Aux::U32(v)) => v as u8,
        Some(Aux::I8(v))   => v as u8,
        Some(Aux::I16(v)) => v as u8,
        Some(Aux::U16(v)) => v as u8,
        Some(_) => panic!("{tag} tag is not u8: {:?}", aux),
    }
}
/// Parse the value of a specified'XX:i:' tag from a BAM record as i32.
/// Return `default` if the tag is not present.
pub fn get_tag_i32_default(aln: &BamRecord, tag: &'static str, default: i32) -> i32 {
    let aux = get_tag_value_opt(aln, tag);
    match aux {
        None => default,
        Some(Aux::I32(v)) => v,
        Some(Aux::U32(v)) => v as i32,
        Some(Aux::I16(v)) => v as i32,
        Some(Aux::U16(v)) => v as i32,
        Some(Aux::U8(v))   => v as i32,
        Some(Aux::I8(v))   => v as i32,
        Some(_) => panic!("{tag} tag is not i32: {:?}", aux),
    }
}
/// Parse the value of a specified'XX:i:' tag from a BAM record as u32.
/// Return `default` if the tag is not present.
pub fn get_tag_u32_default(aln: &BamRecord, tag: &'static str, default: u32) -> u32 {
    let aux = get_tag_value_opt(aln, tag);
    match aux {
        None => default,
        Some(Aux::I32(v)) => v as u32,
        Some(Aux::U32(v)) => v,
        Some(Aux::I16(v)) => v as u32,
        Some(Aux::U16(v)) => v as u32,
        Some(Aux::U8(v))   => v as u32,
        Some(Aux::I8(v))   => v as u32,
        Some(_) => panic!("{tag} tag is not u32: {:?}", aux),
    }
}
/// Parse the value of a specified'XX:f:' tag from a BAM record as f32.
/// Return `default` if the tag is not present.
pub fn get_tag_f32_default(aln: &BamRecord, tag: &'static str, default: f32) -> f32 {
    let aux = get_tag_value_opt(aln, tag);
    match aux {
        None => default,
        Some(Aux::Float(v)) => v,
        Some(Aux::Double(v)) => v as f32,
        Some(_) => panic!("{tag} tag is not a float type: {:?}", aux),
    }
}
/// Parse the value of a specified'XX:A:' tag from a BAM record as bool.
/// Return `default` if the tag is not present.
pub fn get_tag_bool_default(aln: &BamRecord, tag: &'static str, default: bool) -> bool {
    let aux = get_tag_value_opt(aln, tag);
    match aux {
        None => default,
        Some(Aux::Char(v)) => v != 0,
        Some(Aux::U8(v))   => v != 0,
        Some(_) => panic!("{tag} tag is not a bool-compatible type: {:?}", aux),
    }
}
/* ---------------------------------------------------------------------
byte array tag methods
--------------------------------------------------------------------- */
/// Parse the value of a specified 'XX:B:C,' tag from a BAM record as Option<Vec<u8>>.
/// Return None if the tag is not present.
pub fn get_tag_u8_vec_opt(aln: &BamRecord, tag: &'static str) -> Option<Vec<u8>> {
    if let Some(aux) = aln.aux(tag.as_bytes()).ok() {
        return match aux {
            Aux::ArrayU8(v) => Some(v.iter().collect()),
            Aux::ArrayI32(v) => Some(v.iter().map(|x| x as u8).collect()),
            _ => panic!("{tag} tag is not ArrayU8: {:?}", aux),
        }
    }
    None
}
/// Parse the value of a specified 'XX:B:i,' tag from a BAM record as Option<Vec<i32>>.
/// Return None if the tag is not present.
pub fn get_tag_i32_vec_opt(aln: &BamRecord, tag: &'static str) -> Option<Vec<i32>> {
    if let Some(aux) = aln.aux(tag.as_bytes()).ok() {
        return match aux {
            Aux::ArrayI32(v) => Some(v.iter().map(|x| x).collect()),
            _ => panic!("{tag} tag is not ArrayI32: {:?}", aux),
        }
    }
    None
}
/// Parse the value of a specified 'XX:B:i,' tag from a BAM record as Vec<i32>.
/// Panic if tag is absent or cannot be parsed as requested.
pub fn get_tag_i32_vec(aln: &BamRecord, tag: &'static str) -> Vec<i32> {
    let aux = get_tag_value(aln, tag);
    match aux {
        Aux::ArrayI32(v) => v.iter().collect(),
        _ => panic!("{tag} tag is not ArrayI32: {:?}", aux),
    }
}
