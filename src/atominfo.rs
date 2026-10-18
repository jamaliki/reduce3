//! Per-atom contact information for Probe scoring: port of
//! `mmtbx.probe.Helpers.getExtraAtomInfo` (van der Waals or ionic radius,
//! donor/acceptor flags, charge, altloc), including its warning text.

use crate::atomtypes::is_aromatic_acceptor;
use std::fmt::Write;
use crate::interp::{EType, FlatAtoms};
use crate::model::Structure;
use crate::monlib::MonLib;
use crate::probe::AtomInfo;
use crate::resclass::{self, ONE_LETTER_GIVEN_THREE_LETTER};

/// Python `str(float)` for the values that appear in the warnings.
pub fn py_float(v: f64) -> String {
    if !v.is_finite() {
        return if v.is_nan() { "nan".into() } else if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let s = format!("{}", v);
    if s.contains('.') || s.contains('e') { s } else { format!("{}.0", s) }
}

fn py_opt(v: Option<f64>) -> String {
    v.map(py_float).unwrap_or_else(|| "None".into())
}

/// Bonded-neighbor lists (`getBondedNeighborLists`) from bond pairs: each
/// atom's neighbors in increasing atom order, which is the order the sorted
/// cctbx bond proxies produce.
pub fn bonded_lists(n: usize, bonds: impl Iterator<Item = (u32, u32)>) -> Vec<Vec<u32>> {
    let mut v: Vec<Vec<u32>> = vec![Vec::new(); n];
    for (i, j) in bonds {
        v[i as usize].push(j);
        v[j as usize].push(i);
    }
    for l in &mut v {
        l.sort_unstable();
        l.dedup();
    }
    v
}

pub struct ExtraInfo {
    pub info: Vec<AtomInfo>,
    pub warnings: String,
}

/// `getExtraAtomInfo` over all atoms of the structure (flat order).
pub fn extra_atom_info(
    st: &Structure,
    flat: &FlatAtoms,
    etype: &[EType],
    ml: &MonLib,
    bonded: &[Vec<u32>],
    set_polar_hydrogen_radius: bool,
) -> ExtraInfo {
    extra_atom_info_with_mode(st, flat, etype, ml, bonded, set_polar_hydrogen_radius, false, true)
}

/// `extra_atom_info`, with the warning text built only when `report` is set.
pub fn extra_atom_info_report(
    st: &Structure,
    flat: &FlatAtoms,
    etype: &[EType],
    ml: &MonLib,
    bonded: &[Vec<u32>],
    set_polar_hydrogen_radius: bool,
    report: bool,
) -> ExtraInfo {
    extra_atom_info_with_mode(st, flat, etype, ml, bonded, set_polar_hydrogen_radius, false, report)
}

/// Probe validation typing with CCTBX's electron or neutron energy radii.
/// Polar hydrogens retain Probe's 1.05 A override in both modes.
pub fn extra_atom_info_validation(
    st: &Structure,
    flat: &FlatAtoms,
    etype: &[EType],
    ml: &MonLib,
    bonded: &[Vec<u32>],
    nuclear: bool,
) -> ExtraInfo {
    extra_atom_info_with_mode(st, flat, etype, ml, bonded, true, nuclear, true)
}

fn extra_atom_info_with_mode(
    st: &Structure,
    flat: &FlatAtoms,
    etype: &[EType],
    ml: &MonLib,
    bonded: &[Vec<u32>],
    set_polar_hydrogen_radius: bool,
    nuclear: bool,
    report: bool,
) -> ExtraInfo {
    let n = flat.pos.len();
    let mut info = Vec::with_capacity(n);
    let mut w = String::new();
    // The arguments are formatted only when the warning text is wanted.
    macro_rules! note {
        ($($arg:tt)*) => {
            if report {
                let _ = writeln!(w, $($arg)*);
            }
        };
    }
    for k in 0..n {
        let pth = flat.path[k];
        let at = st.atom(pth);
        let ag = st.atom_group(pth);
        let rg = st.residue_group(pth);
        let ch = st.chain(pth);
        let elem = flat.element[k].as_str();
        let name = at.name.trim();
        let full_name = || format!("{} {} {} {}", ch.id, ag.resname.trim(), rg.resseq_as_int(), name);
        let mut e = AtomInfo::default();
        e.vdw_radius = 0.0;
        let alt = ag.altloc.trim();
        e.alt = alt.bytes().next().unwrap_or(0);
        e.charge = at.charge_value() as i8;
        e.is_ion = resclass::element_is_ion(elem);
        // energy-type lookup; a failure takes the exception path of the original
        let entry = match &etype[k] {
            EType::Typed(t) => ml.ener.get(t.as_str()).ok_or_else(|| format!("'{}'", t)),
            EType::NoneType => Err("'None'".to_string()),
            EType::Unexpected => Err("'False'".to_string()),
            // cctbx's energy library has an all-empty entry for the empty type
            EType::Unknown => ml.ener.get("").ok_or_else(|| "''".to_string()),
        };
        let radius: Result<(), String> = entry.as_ref().map_err(|m| m.clone()).and_then(|en| {
            let vdw = energy_radius(en, nuclear);
            if e.is_ion {
                match en.ion_radius {
                    Some(r) => {
                        e.vdw_radius = r;
                        note!( "Using ionic radius for {}: {} (rather than {})", name, py_float(r), py_opt(vdw));
                        Ok(())
                    }
                    None => Err(none_assignment_message()),
                }
            } else {
                match vdw {
                    Some(r) => {
                        e.vdw_radius = r;
                        Ok(())
                    }
                    None => Err(none_assignment_message()),
                }
            }
        });
        if let Err(msg) = radius {
            note!(
                "Warning: Could not find atom info for {} (perhaps interpretation was not run on the model?): keeping some default values: {}",
                full_name(),
                msg
            );
            if e.vdw_radius <= 0.0 {
                if let Some(r) = ml.ener.get(elem).and_then(|x| energy_radius(x, nuclear)).filter(|&r| r != 0.0) {
                    e.vdw_radius = r;
                    note!( "Using element-based VDW radius for {}: {}", full_name(), py_float(r));
                }
            }
            info.push(e);
            continue;
        }
        let en = entry.unwrap();
        if (elem == "C" || elem == "N") && is_aromatic_acceptor(&ag.resname, &at.name) {
            e.is_acceptor = true;
            note!( "Marking {} as an aromatic-ring acceptor", name);
        }
        if elem == "N" && !ONE_LETTER_GIVEN_THREE_LETTER.iter().any(|(r, _)| *r == ag.resname) {
            let has_h = bonded[k].iter().any(|&j| matches!(flat.element[j as usize].as_str(), "H" | "D"));
            if !has_h && !e.is_acceptor {
                e.is_acceptor = true;
                note!( "Marking {} {} as a non-Hydrogen HET acceptor", ag.resname.trim(), name);
            }
        }
        if name.eq_ignore_ascii_case("C") || special_carbonyl(&ag.resname, &at.name) {
            let expected = 1.65;
            if e.vdw_radius != expected {
                note!( "Overriding radius for {}: {} (was {})", name, py_float(expected), py_float(e.vdw_radius));
                e.vdw_radius = expected;
            }
        }
        let polar_h = set_polar_hydrogen_radius
            && matches!(elem, "H" | "D")
            && bonded[k].len() == 1
            && matches!(flat.element[bonded[k][0] as usize].as_str(), "N" | "O" | "S");
        if polar_h && e.vdw_radius != 1.05 {
            note!( "Overriding radius for {}: 1.05 (was {})", name, py_float(e.vdw_radius));
            e.vdw_radius = 1.05;
        }
        let hb = en.hb_type.as_str();
        if matches!(hb, "A" | "B" | "D" | "N" | "H") {
            if hb == "A" || hb == "B" {
                e.is_acceptor = true;
            }
            if hb == "D" || hb == "B" {
                e.is_donor = true;
            }
        } else {
            note!( "Warning: Unrecognized specific H bond type for {}, got {}: keeping default value", full_name(), hb);
        }
        info.push(e);
    }
    ExtraInfo { info, warnings: w }
}

/// `is_special_amino_acid_carbonyl(resname.trim() upper-cased, atom name upper-cased, true)`,
/// compared without copying either name.
fn special_carbonyl(resname: &str, atom_name: &str) -> bool {
    let residues: &[&str] = match atom_name.trim_end() {
        name if name.eq_ignore_ascii_case(" CG") => &["ASP", "ASN", "ASX"],
        name if name.eq_ignore_ascii_case(" CD") => &["GLU", "GLN", "GLX"],
        _ => return false,
    };
    let resname = resname.trim();
    residues.iter().any(|value| resname.eq_ignore_ascii_case(value))
}

fn energy_radius(entry: &crate::monlib::EnerAtom, nuclear: bool) -> Option<f64> {
    if nuclear {
        entry.vdw_radius_neutron.or(entry.vdw_radius)
    } else {
        entry.vdw_radius
    }
}

/// The Boost.Python message for assigning None to a double property.
fn none_assignment_message() -> String {
    "Python argument types in\n    None.None(ExtraAtomInfo, NoneType)\ndid not match C++ signature:\n    None(molprobity::probe::ExtraAtomInfo {lvalue}, double)".to_string()
}
