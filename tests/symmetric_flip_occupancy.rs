//! Fixed mode swaps a symmetric side chain's occupancy and anisotropic ADP
//! with its coordinates and B when it flips the atom names (fix 28 in the
//! README). Needs the monomer library (found as the command-line tool finds
//! it); skipped without one.

#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::geom::{v3, Vec3};
use reduce3::{mmcif, pipeline, Params};

/// Five residues of 1QVC chain B. The dihedral CA-CB-CG-OD1 of Asp 290 is
/// 167 degrees, so Reduce3 renames its carboxylate oxygens: OD1 is first at
/// (1.657, 6.904, 22.153) with occupancy 1.00 and B 57.11, OD2 at
/// (1.402, 8.862, 23.129) with occupancy 0.60 and B 54.95.
const MODEL: &str = include_str!("data/1qvc_b289_293.cif");

fn model() -> String {
    let partial = "ATOM 15 O OD2 . ASP B 1 90 ? 1.402 8.862 23.129 0.6 54.95";
    let text = MODEL.replace("ATOM 15 O OD2 . ASP B 1 90 ? 1.402 8.862 23.129 1 54.95", partial);
    assert!(text.contains(partial), "the fixture's OD2 row changed");
    text
}

/// The same model with anisotropic ADPs on the two oxygens (atom ids 14 and 15).
fn anisotropic_model() -> String {
    model()
        + "
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
14 0.80 0.70 0.60 0.01 0.02 0.03
15 0.50 0.40 0.30 -0.01 -0.02 -0.03
"
}

/// Occupancy, B, position and anisotropic ADP of an atom of residue `resseq`.
fn atom(output: &pipeline::Output, resseq: &str, name: &str) -> (f64, f64, Vec3, Option<[f64; 6]>) {
    let st = &output.structure;
    let path = st
        .atom_paths()
        .into_iter()
        .find(|&p| st.residue_group(p).resseq.trim() == resseq && st.atom(p).name.trim() == name)
        .unwrap_or_else(|| panic!("no atom {name} in residue {resseq}"));
    let a = st.atom(path);
    (a.occ, a.b, a.xyz, a.uij)
}

fn run(text: &str, params: &Params) -> Option<pipeline::Output> {
    let library = chem_data::library()?;
    Some(pipeline::run(mmcif::read_mmcif(text).unwrap(), &library, params).unwrap())
}

#[test]
fn flipped_carboxylate_oxygens_keep_their_occupancy() {
    let Some(out) = run(&model(), &Params::default()) else { return };
    let (occ1, b1, p1, _) = atom(&out, "290", "OD1");
    let (occ2, b2, p2, _) = atom(&out, "290", "OD2");
    // the flip put OD1 where OD2 was and the other way round ...
    assert!(p1.dist(v3(1.402, 8.862, 23.129)) < 1e-9, "OD1 at {p1:?}");
    assert!(p2.dist(v3(1.657, 6.904, 22.153)) < 1e-9, "OD2 at {p2:?}");
    // ... and the B factor and occupancy of each site went with it
    assert_eq!((b1, b2), (54.95, 57.11));
    assert_eq!((occ1, occ2), (0.6, 1.0), "occupancy stayed with the atom name, not the site");
}

#[test]
fn flipped_carboxylate_oxygens_keep_their_anisotropic_adp() {
    let Some(out) = run(&anisotropic_model(), &Params::default()) else { return };
    let (occ1, b1, p1, u1) = atom(&out, "290", "OD1");
    let (occ2, b2, p2, u2) = atom(&out, "290", "OD2");
    assert!(p1.dist(v3(1.402, 8.862, 23.129)) < 1e-9 && p2.dist(v3(1.657, 6.904, 22.153)) < 1e-9);
    assert_eq!((occ1, occ2), (0.6, 1.0));
    // the site that was atom 15 (U about 0.4, the lower trace) is now OD1
    let trace = |u: Option<[f64; 6]>| u.map(|u| (u[0] + u[1] + u[2]) / 3.0).unwrap();
    assert!((trace(u1) - 0.4).abs() < 1e-3, "OD1 trace {}", trace(u1));
    assert!((trace(u2) - 0.7).abs() < 1e-3, "OD2 trace {}", trace(u2));
    // an anisotropic atom's B is its equivalent isotropic B, so it follows the ADP
    let eight_pi_sq = 8.0 * std::f64::consts::PI * std::f64::consts::PI;
    assert!((b1 - trace(u1) * eight_pi_sq).abs() < 1e-6 && (b2 - trace(u2) * eight_pi_sq).abs() < 1e-6, "B {b1} {b2}");
}
