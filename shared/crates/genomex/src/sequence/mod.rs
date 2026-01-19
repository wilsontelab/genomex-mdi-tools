/// Genetic sequence-related utilities and functions, including
/// reverse complement calculation, IUPAC code handling,
/// and sequence alignment.

// modules
mod rev_comp;
mod iupac;
mod smith_waterman;

// re-exports
pub use rev_comp::*;
pub use iupac::*;
pub use smith_waterman::*;
