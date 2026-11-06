//! Fixed mode keeps a hydroxyl, thiol or ammonium hydrogen at least 1.5 A from
//! every heavy atom that is not its parent when some turn about its bond does
//! (fix 29 in the README). Reduce2 turned such a hydrogen straight at a close
//! acceptor and left it as near as 1.2 A. Needs the monomer library (found as
//! the command-line tool finds it); skipped without one.
//!
//! The first test is a real case, five residues of 1QVC chain B. The others
//! put a water, or several, next to a residue of that slice.

#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::geom::{rotate_deg_axis_dir, v3, Vec3};
use reduce3::hplace::NTermCharge;
use reduce3::{mmcif, pipeline, MonLib, Params};

const SLICE: &str = include_str!("data/1qvc_b289_293.cif");

/// Where a hydrogen may not be: the clearance the guard promises, less rounding.
const CLEAR: f64 = 1.5 - 1e-6;
/// Closer than this is the defect: a hydrogen aimed at an acceptor 2.1 A or
/// less from its parent.
const CLOSE: f64 = 1.4;

// Coordinates of the slice's Thr 289 OG1 and CB, and Ser 292 OG and CB.
const THR_OG1: Vec3 = v3(4.967, 14.184, 25.034);
const THR_CB: Vec3 = v3(4.203, 13.753, 23.903);
const THR_N: Vec3 = v3(6.192, 12.635, 22.959);
const SER_CB: Vec3 = v3(-2.224, 11.447, 25.877);
const SER_OG: Vec3 = v3(-2.349, 10.21, 25.206);

fn run(text: &str, params: &Params) -> Option<pipeline::Output> {
    let root = chem_data::locate()?;
    let library = MonLib::load(&root).unwrap();
    Some(pipeline::run(mmcif::read_mmcif(text).unwrap(), &library, params).unwrap())
}

fn params(clearance: f64) -> Params {
    let mut p = Params::default();
    p.opt.rotor_clearance = clearance;
    p
}

/// The slice with waters added: (position, occupancy, B).
fn with_waters(text: &str, waters: &[(Vec3, f64, f64)]) -> String {
    let mut rows = String::new();
    for (k, (xyz, occ, b)) in waters.iter().enumerate() {
        let (id, seq) = (100 + k, 101 + k);
        rows += &format!(
            "HETATM {id} O O . HOH B 2 {seq} ? {:.3} {:.3} {:.3} {occ} {b} ? {seq} B 1\n",
            xyz.x, xyz.y, xyz.z
        );
    }
    let at = text.find("ATOM 34 ").unwrap();
    let end = at + text[at..].find('\n').unwrap() + 1;
    let mut out = text.to_string();
    out.insert_str(end, &rows);
    let out = out.replace("THR .\n", "THR .\nHOH .\n");
    assert!(out.contains("HOH .\n") || waters.is_empty());
    out
}

/// The slice's Ser 292 turned into a cysteine, its SG one bond length out
/// along CB-OG.
fn cysteine_slice() -> String {
    let sg = SER_CB + (SER_OG - SER_CB).normalize() * 1.81;
    let ser_og = "ATOM 30 O OG . SER B 1 92 ? -2.349 10.21 25.206 1 50.62 ? 292 B 1";
    let cys_sg = format!("ATOM 30 S SG . CYS B 1 92 ? {:.3} {:.3} {:.3} 1 50.62 ? 292 B 1", sg.x, sg.y, sg.z);
    assert!(SLICE.contains(ser_og));
    let text = SLICE.replace(ser_og, &cys_sg);
    let text = text.replace(" SER B 1 92 ", " CYS B 1 92 ").replace("SER .\n", "CYS .\n");
    text.replace("_atom_type.symbol\nC\nN\nO\n", "_atom_type.symbol\nC\nN\nO\nS\n")
}

fn position(output: &pipeline::Output, resseq: &str, name: &str) -> Vec3 {
    let st = &output.structure;
    let path = st
        .atom_paths()
        .into_iter()
        .find(|&p| st.residue_group(p).resseq.trim() == resseq && st.atom(p).name.trim() == name)
        .unwrap_or_else(|| panic!("no atom {name} in residue {resseq}"));
    st.atom(path).xyz
}

/// The nearest heavy atom to atom `name` of residue `resseq` other than the
/// atoms of that residue named in `own` (its parent and the atom behind it),
/// and its distance.
fn nearest_heavy(output: &pipeline::Output, resseq: &str, name: &str, own: &[&str]) -> (String, f64) {
    let st = &output.structure;
    let h = position(output, resseq, name);
    st.atom_paths()
        .into_iter()
        .filter(|&p| st.atom(p).element.trim() != "H" && st.atom(p).element.trim() != "D")
        .filter(|&p| !(st.residue_group(p).resseq.trim() == resseq && own.contains(&st.atom(p).name.trim())))
        .map(|p| (format!("{} {}", st.residue_group(p).resseq.trim(), st.atom(p).name.trim()), st.atom(p).xyz.dist(h)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap()
}

/// Every atom's name and position, in output order.
fn coordinates(output: &pipeline::Output) -> Vec<(String, [f64; 3])> {
    let st = &output.structure;
    st.atom_paths()
        .into_iter()
        .map(|p| (format!("{} {}", st.residue_group(p).resseq.trim(), st.atom(p).name.trim()), st.atom(p).xyz.as_array()))
        .collect()
}

/// The atoms whose positions differ, by more than the rounding a Mover's
/// rotation of its own parent atom leaves.
fn moved(a: &[(String, [f64; 3])], b: &[(String, [f64; 3])]) -> Vec<String> {
    assert_eq!(a.len(), b.len());
    let far = |p: &[f64; 3], q: &[f64; 3]| (0..3).any(|k| (p[k] - q[k]).abs() > 1e-9);
    a.iter().zip(b).filter(|(x, y)| x.0 != y.0 || far(&x.1, &y.1)).map(|(x, _)| x.0.clone()).collect()
}

/// A point `distance` from `at`, 109.5 degrees from the direction to `from`
/// and `torsion` degrees round the bond: where the acceptor of a hydroxyl
/// lies when its hydrogen is turned the opposite way.
fn around(at: Vec3, from: Vec3, torsion: f64, distance: f64) -> Vec3 {
    let axis = (at - from).normalize();
    let side = axis.cross(v3(0.3, 0.5, 0.8)).normalize();
    let angle = 70.5_f64.to_radians();
    let point = at + (axis * angle.cos() + side * angle.sin()) * distance;
    rotate_deg_axis_dir(at, axis, point, torsion)
}

#[test]
fn a_real_serine_hydroxyl_hydrogen_is_not_aimed_at_the_carbonyl_oxygen_it_can_avoid() {
    // Ser 292 OG is 2.06 A from the carbonyl O of Gln 291. Reduce2 puts HG
    // 1.27 A from that oxygen, although a turn that leaves it 1.87 A away exists.
    let Some(clear) = run(SLICE, &params(1.5)) else { return };
    let (who, d) = nearest_heavy(&clear, "292", "HG", &["OG", "CB"]);
    assert!(d >= CLEAR, "HG is {d:.3} A from {who}");
    // 1.5 A is the default
    let default = run(SLICE, &Params::default()).unwrap();
    assert_eq!(moved(&coordinates(&default), &coordinates(&clear)), Vec::<String>::new());
    let off = run(SLICE, &params(0.0)).unwrap();
    let (who, d) = nearest_heavy(&off, "292", "HG", &["OG", "CB"]);
    assert!(d < CLOSE && who == "291 O", "without the guard HG is {d:.3} A from {who}");
    // compat mode reproduces Reduce2 and keeps its pose
    let compat = run(SLICE, &Params { compat: true, ..Params::default() }).unwrap();
    let (who, d) = nearest_heavy(&compat, "292", "HG", &["OG", "CB"]);
    assert!(d < CLOSE && who == "291 O", "in compat mode HG is {d:.3} A from {who}");
    // nothing else moves: this slice has no other rotor with a blocked turn
    assert_eq!(moved(&coordinates(&clear), &coordinates(&off)), ["292 HG"]);
}

#[test]
fn a_threonine_hydroxyl_hydrogen_keeps_clear_of_a_water_at_the_oxygen() {
    // a water 2.0 A from OG1, where the hydrogen could point straight at it
    let text = with_waters(SLICE, &[(around(THR_OG1, THR_CB, 0.0, 2.0), 1.0, 20.0)]);
    let Some(clear) = run(&text, &params(1.5)) else { return };
    let (who, d) = nearest_heavy(&clear, "289", "HG1", &["OG1", "CB"]);
    assert!(d >= CLEAR, "HG1 is {d:.3} A from {who}");
    let off = run(&text, &params(0.0)).unwrap();
    let (who, d) = nearest_heavy(&off, "289", "HG1", &["OG1", "CB"]);
    assert!(d < CLOSE && who == "101 O", "without the guard HG1 is {d:.3} A from {who}");
}

#[test]
fn a_water_the_scoring_ignores_still_keeps_the_hydrogen_out() {
    // A water with B 60 is left out of the scoring, so Reduce2 does not see it.
    // Put one right where the hydrogen goes without it.
    let Some(bare) = run(SLICE, &params(0.0)) else { return };
    let (og1, hg1) = (position(&bare, "289", "OG1"), position(&bare, "289", "HG1"));
    let water = og1 + (hg1 - og1).normalize() * 2.0;
    let text = with_waters(SLICE, &[(water, 1.0, 60.0)]);
    let off = run(&text, &params(0.0)).unwrap();
    assert_eq!(position(&off, "289", "HG1"), hg1, "the scoring is meant not to see the water");
    assert!(position(&off, "289", "HG1").dist(water) < CLOSE);
    let clear = run(&text, &params(1.5)).unwrap();
    let (who, d) = nearest_heavy(&clear, "289", "HG1", &["OG1", "CB"]);
    assert!(d >= CLEAR, "HG1 is {d:.3} A from {who}");
}

#[test]
fn ammonium_hydrogens_keep_clear_of_a_water() {
    // the first residue of the chain as an N-terminus, a water 2.1 A from N
    // straight out along one of the hydrogens Reduce2 places
    let mut p = params(0.0);
    p.n_terminal_charge = NTermCharge::FirstInChain;
    let Some(bare) = run(SLICE, &p) else { return };
    let h1 = position(&bare, "289", "H1");
    let water = THR_N + (h1 - THR_N).normalize() * 2.1;
    let text = with_waters(SLICE, &[(water, 1.0, 20.0)]);
    let off = run(&text, &p).unwrap();
    let closest = |out: &pipeline::Output| {
        ["H1", "H2", "H3"].map(|h| position(out, "289", h).dist(water)).into_iter().fold(f64::INFINITY, f64::min)
    };
    assert!(closest(&off) < CLOSE, "without the guard the nearest NH3 hydrogen is {:.3} A from the water", closest(&off));
    p.opt.rotor_clearance = 1.5;
    let clear = run(&text, &p).unwrap();
    assert!(closest(&clear) >= CLEAR, "the nearest NH3 hydrogen is {:.3} A from the water", closest(&clear));
    // the group is still staggered on the bond: the N-H bonds are unchanged
    for h in ["H1", "H2", "H3"] {
        let bond = position(&clear, "289", h).dist(THR_N);
        assert!((bond - position(&off, "289", h).dist(THR_N)).abs() < 1e-6);
    }
}

#[test]
fn a_thiol_hydrogen_keeps_clear_of_a_water() {
    let slice = cysteine_slice();
    let Some(bare) = run(&slice, &params(0.0)) else { return };
    let (sg, hg) = (position(&bare, "292", "SG"), position(&bare, "292", "HG"));
    let water = sg + (hg - sg).normalize() * 2.4;
    let text = with_waters(&slice, &[(water, 1.0, 20.0)]);
    let off = run(&text, &params(0.0)).unwrap();
    assert!(position(&off, "292", "HG").dist(water) < CLOSE);
    let clear = run(&text, &params(1.5)).unwrap();
    let (who, d) = nearest_heavy(&clear, "292", "HG", &["SG", "CB"]);
    assert!(d >= CLEAR, "HG is {d:.3} A from {who}");
}

#[test]
fn a_hydrogen_with_no_clear_turn_is_placed_as_without_the_guard() {
    // six waters 2.0 A from Ser 292 OG, round the cone its hydrogen sweeps:
    // every turn leaves HG within 1.5 A of one of them. Nothing in the slice
    // is near enough to give another rotor a blocked turn.
    let waters: Vec<_> = (0..6).map(|k| (around(SER_OG, SER_CB, 60.0 * k as f64, 2.0), 1.0, 20.0)).collect();
    let text = with_waters(SLICE, &waters);
    let Some(clear) = run(&text, &params(1.5)) else { return };
    let off = run(&text, &params(0.0)).unwrap();
    let (who, d) = nearest_heavy(&off, "292", "HG", &["OG", "CB"]);
    assert!(d < CLOSE, "the case is meant to leave no clear turn, but HG is {d:.3} A from {who}");
    assert_eq!(moved(&coordinates(&clear), &coordinates(&off)), Vec::<String>::new());
}

#[test]
fn an_ordinary_hydrogen_bond_is_left_alone() {
    // a water 2.8 A from OG1, a normal hydrogen bond: no turn is blocked
    let text = with_waters(SLICE, &[(around(THR_OG1, THR_CB, 0.0, 2.8), 1.0, 20.0)]);
    let Some(clear) = run(&text, &params(1.5)) else { return };
    let off = run(&text, &params(0.0)).unwrap();
    assert_eq!(moved(&coordinates(&clear), &coordinates(&off)), ["292 HG"]);
    assert!(position(&clear, "289", "HG1").dist(position(&clear, "101", "O")) > 1.6);
}
