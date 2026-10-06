//! Reduce3: a fast Rust reimplementation of cctbx Reduce2 (hydrogen addition
//! and optimization of rotatable and flippable groups).
//!
//! The command-line program is a thin layer over this library. To run Reduce3
//! on a model another program has already parsed, implement
//! [`cifsource::CifSource`] for its mmCIF data block and call [`run_cif`]; the
//! result can be streamed back into that program's own document type through
//! [`cifsource::CifSink`] with [`mmcif::write_cif`].

pub mod atominfo;
pub mod atomtypes;
pub mod autolink;
pub mod ccdrestraints;
pub mod cell;
pub mod cif;
pub mod cifsource;
pub mod geom;
pub mod hplace;
pub mod interp;
mod linkdata;
mod rdkit_valence;
pub mod monlib;
pub mod names;
pub mod riding;
pub mod mmcif;
pub mod model;
pub mod movers;
pub mod optimizer;
pub mod pipeline;
pub mod pdbio;
pub mod probe;
pub mod resclass;
mod resclass_data;
mod spacegroup_data;
#[cfg(feature = "refcheck")]
pub mod refcheck;
pub mod world;

pub use cifsource::{CifCell, CifSink, CifSource, CifTable};
pub use monlib::MonLib;
pub use pipeline::{Approach, Output, Params};

/// Run Reduce3 on a parsed mmCIF data block. The block is only read; the
/// result's `structure` can be written with [`mmcif::write_cif`] or
/// [`mmcif::write_mmcif`].
pub fn run_cif<S: CifSource + ?Sized>(block: &S, monlib: &MonLib, params: &Params) -> Result<Output, String> {
    let st = mmcif::structure_from_cif(block)?;
    if st.atoms_size() == 0 {
        return Err("the model has no atoms".to_string());
    }
    pipeline::run(st, monlib, params)
}
