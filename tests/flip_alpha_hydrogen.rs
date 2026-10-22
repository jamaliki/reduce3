//! Fixed mode places HA again after an Asn, Gln or His flip has docked the
//! side chain about CA (fix 27 in the README). The model is every residue and
//! water of 1A28 chain A within 4 A of the side chains of Asn 689, Gln 725 and
//! His 881, all three of which Reduce3 flips there as in the full entry.

use reduce3::geom::Vec3;
use reduce3::hplace::{self, HPlaceParams, NTermCharge};
use reduce3::model::Structure;
use reduce3::riding::{compute_h_position, RidingCoef};
use reduce3::{mmcif, pipeline, MonLib, Params};
use std::collections::HashMap;

const MODEL: &str = include_str!("data/flips_1a28.cif");

/// Chain, residue number, altloc and atom name.
type Key = (String, String, String, String);

fn key(st: &Structure, p: reduce3::model::AtomPath) -> Key {
    let ag = st.atom_group(p);
    (st.chain(p).id.clone(), st.residue_group(p).resseq.trim().to_string(), ag.altloc.clone(), st.atom(p).name.trim().to_string())
}

fn atoms(st: &Structure) -> Vec<(Key, Vec3)> {
    st.atom_paths().into_iter().map(|p| (key(st, p), st.atom(p).xyz)).collect()
}

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("set REDUCE3_CHEM_DATA for the flip tests")).unwrap()
}

fn run(library: &MonLib, flip_states: &str, reride: bool) -> pipeline::Output {
    let mut params = Params::default();
    params.opt.flip_states = flip_states.into();
    params.opt.reride_flipped_alpha_h = reride;
    pipeline::run(mmcif::read_mmcif(MODEL).unwrap(), library, &params).unwrap()
}

/// The riding coefficients placement gives the model (as the pipeline places
/// it), with each atom's key in their order.
fn riding(library: &MonLib) -> (Vec<Key>, Vec<Option<RidingCoef>>) {
    let mut st = mmcif::read_mmcif(MODEL).unwrap();
    let placed = hplace::place_hydrogens(&mut st, library, &HPlaceParams {
        neutron: false, n_terminal_charge: NTermCharge::ResidueOne, exclude_water: true,
        keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
    });
    (st.atom_paths().into_iter().map(|p| key(&st, p)).collect(), placed.riding)
}

/// Distance of each hydrogen of residue `resseq` in `output` from its riding
/// position, by `compute_h_position` from the output's coordinates.
fn riding_errors(output: &[(Key, Vec3)], keys: &[Key], coef: &[Option<RidingCoef>], resseq: &str) -> HashMap<String, f64> {
    let at: HashMap<&Key, Vec3> = output.iter().map(|(k, x)| (k, *x)).collect();
    let sites: Vec<Vec3> = keys.iter().map(|k| at.get(k).copied().unwrap_or_default()).collect();
    keys.iter().zip(coef).filter(|(k, c)| k.1 == resseq && c.is_some() && at.contains_key(k))
        .map(|(k, c)| (k.3.clone(), compute_h_position(c.as_ref().unwrap(), &sites).unwrap().dist(at[k])))
        .collect()
}

/// Residue numbers and names of the atoms of `new` whose coordinates are not
/// bitwise those of `old`; both must have the same atoms in the same order.
fn moved(old: &[(Key, Vec3)], new: &[(Key, Vec3)]) -> Vec<(String, String)> {
    assert_eq!(old.iter().map(|a| &a.0).collect::<Vec<_>>(), new.iter().map(|a| &a.0).collect::<Vec<_>>());
    let bits = |v: Vec3| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()];
    let mut moved: Vec<_> =
        old.iter().zip(new).filter(|(a, b)| bits(a.1) != bits(b.1)).map(|(a, _)| (a.0.1.clone(), a.0.3.clone())).collect();
    moved.sort();
    moved
}

fn names(atoms: &[(&str, &str)]) -> Vec<(String, String)> {
    atoms.iter().map(|(r, n)| (r.to_string(), n.to_string())).collect()
}

/// The description without its timing lines.
fn report(output: &pipeline::Output) -> Vec<&str> {
    output.description.lines().filter(|l| !l.contains("Time to")).collect()
}

fn flipped(output: &pipeline::Output, mover: &str) -> bool {
    output.description.lines().any(|l| l.contains(mover) && l.contains("pose Flipped") && l.contains("AnglesAdjusted"))
}

#[test]
fn optimizer_flips_place_ha_from_the_docked_side_chain() {
    let library = library();
    let (keys, coef) = riding(&library);
    let new = run(&library, "", true);
    let old = run(&library, "", false);
    for mover in ["AmideFlip at chain A ASN 689 ND2", "AmideFlip at chain A GLN 725 NE2", "HisFlip at chain A HIS 881 NE2"] {
        assert!(flipped(&new, mover), "{mover} is not flipped");
    }
    assert!(!flipped(&new, "ASN 828"));
    // the scores and choices do not change, and nothing but those three HA moves
    assert_eq!(report(&new), report(&old));
    let (new_atoms, old_atoms) = (atoms(&new.structure), atoms(&old.structure));
    assert_eq!(moved(&old_atoms, &new_atoms), names(&[("689", "HA"), ("725", "HA"), ("881", "HA")]));
    for resseq in ["689", "725", "881"] {
        let now = riding_errors(&new_atoms, &keys, &coef, resseq);
        let before = riding_errors(&old_atoms, &keys, &coef, resseq);
        assert!(now["HA"] < 1e-9, "HA of {resseq} is {} A from its riding position", now["HA"]);
        // Reduce2 left HA where the unflipped CB put it
        assert!(before["HA"] > 0.01 && before["HA"] < 0.05, "HA of {resseq} was {} A off", before["HA"]);
        // the CB and CG hydrogens turn with the rest of the side chain
        for (name, error) in &now {
            if name.starts_with("HB") || name.starts_with("HG") {
                assert!(*error < 1e-9, "{name} of {resseq} is {error} A from its riding position");
            }
        }
    }
}

#[test]
fn forced_flips_place_ha_only_when_the_side_chain_is_docked() {
    let library = library();
    let (keys, coef) = riding(&library);
    // Asn 828, which the optimizer leaves, flipped and docked when Movers are
    // placed; His 881 docked by the fix-up of its fixed Mover; Gln 725 flipped
    // without the fix-up and Asn 689 kept unflipped, so neither is docked
    let states = "1 . A ASN 828 Flipped AnglesAdjusted,1 . A HIS 881 Flipped AnglesAdjusted,\
        1 . A GLN 725 Flipped,1 . A ASN 689 Unflipped";
    let new = run(&library, states, true);
    let old = run(&library, states, false);
    assert_eq!(report(&new), report(&old));
    let (new_atoms, old_atoms) = (atoms(&new.structure), atoms(&old.structure));
    assert_eq!(moved(&old_atoms, &new_atoms), names(&[("828", "HA"), ("881", "HA")]));
    for resseq in ["828", "881"] {
        assert!(riding_errors(&new_atoms, &keys, &coef, resseq)["HA"] < 1e-9);
    }
}
