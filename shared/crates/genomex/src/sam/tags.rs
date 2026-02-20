//! Module with SAM tag constants and utilities.

// dependencies
use std::str::FromStr;
use std::fmt;
use serde::{Deserialize, Serialize};

/// SamFlag struct for working with SAM flag values.
#[derive(Serialize, Deserialize, Debug)]
pub struct SamTags{
    #[serde(deserialize_with = "mdi::record::trailing_to_vec_string")]
    pub tags: Vec<String>,
}
impl SamTags {
    /// Create a new SamTags instance from a vector of tag strings.
    pub fn new(tags: Vec<String>) -> Self {
        SamTags { tags }
    }
    /* -------------------------------------------------------------------------
    tag addition or update methods
    ------------------------------------------------------------------------- */
    /// Add or update a tag in the SamTags Vec<String> using a prefix with
    /// its data type, e.g., either "AS:i" or "AS:i:".
    pub fn set_tag_value(&mut self, prefix: &str, value: &dyn fmt::Display) {
        let mut prefix = prefix.to_string();
        if !prefix.ends_with(':') && !prefix.ends_with(',') {
            prefix = format!("{}:", prefix);
        }
        let mut found = false;
        for i in 0..self.tags.len() {
            if self.tags[i].starts_with(&prefix) {
                self.tags[i] = format!("{}{}", prefix, value);
                found = true;
                break;
            }
        }
        if !found {
            self.tags.push(format!("{}{}", prefix, value)); 
        }
    }
    /* -------------------------------------------------------------------------
    tag retrieval methods
    ------------------------------------------------------------------------- */
    /// Extract a parsed tag value from a SamTags Vec<String> as Some(T),
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., `"AS:"` or `"AS:i:"`. Return None if the tag is not present.
    pub fn get_tag_value_parsed<T>(&self, prefix: &str) -> Option<T>
    where
        T: FromStr,
        <T as FromStr>::Err: fmt::Debug,
    {
        for tag_str in &self.tags {
            if tag_str.starts_with(prefix) {
                // split into 3 parts: TAG, TYPE, VALUE
                //   matched on TAG or TAG:TYPE
                //   ignore TYPE (expect caller to enforce it through T and to handle parse errors)
                //   VALUE may contain colons
                let n_skip = if prefix.starts_with("ML") { 3 } else { 2 };
                let parts = tag_str.splitn(n_skip + 1, ':');
                if let Some(value_str) = parts.skip(n_skip).next() {
                    match value_str.parse::<T>() {
                        Ok(val) => return Some(val),
                        Err(e) => panic!("Failed to parse tag prefix '{}' value '{}': {:?}", prefix, value_str, e),
                    }
                }
            }
        }
        None
    }
    /// Extract a parsed tag value from a SamTags Vec<String> as Some(String),
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., `"AS:"` or `"AS:i:"`. Return None if the tag is not present.
    pub fn get_tag_value(&self, prefix: &str) -> Option<String>{
        self.get_tag_value_parsed::<String>(prefix)
    }
    /* -------------------------------------------------------------------------
    tag removal methods
    ------------------------------------------------------------------------- */
    /// Retain only those tags whose prefix is in the provided &[&str],
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., either `"AS:"` or `"AS:i:"`.
    /// 
    /// Tags on the retention list not found in the SamTags are simply ignored.
    pub fn retain(&mut self, retain_prefixes: &[&str]) {
        self.tags.retain(|tag_str| {
            retain_prefixes.iter().any(|retain_prefix| tag_str.starts_with(retain_prefix))
        });
    }   
    /// Drop tags whose prefix is in the provided &[&str],
    /// using the two-letter tag ID with or without the data type prefix, 
    /// e.g., either `"AS:"` or `"AS:i:"`.
    /// 
    /// Tags on the drop list not found in the SamTags are simply ignored.
    pub fn drop(&mut self, drop_prefixes: &[&str]) {
        self.tags.retain(|tag_str| {
            !drop_prefixes.iter().any(|drop_prefix| tag_str.starts_with(drop_prefix))
        });
    }   
}
