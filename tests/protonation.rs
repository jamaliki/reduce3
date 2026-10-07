//! Fixed mode's protonation fixes: no histidine without a ring hydrogen, and
//! no thiol hydrogen on a cysteine bound to a metal other than zinc. Needs the
//! monomer library (found as the command-line tool finds it); skipped without one.

use reduce3::{mmcif, pipeline, MonLib, Params};

const HEADER: &str = "\
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
";

/// Cys 507 of ascorbate oxidase (1ASQ) and its type-1 copper, 2.08 A from SG.
const CYS_CU: &str = "\
ATOM 11 N N . CYS A 1 . ? 38.735 18.444 -0.540 1.00 6.00 ? 507 CYS A N 1
ATOM 12 C CA . CYS A 1 . ? 39.040 18.863 0.814 1.00 7.00 ? 507 CYS A CA 1
ATOM 13 C C . CYS A 1 . ? 40.525 18.522 1.116 1.00 8.78 ? 507 CYS A C 1
ATOM 14 O O . CYS A 1 . ? 40.966 17.404 0.857 1.00 6.00 ? 507 CYS A O 1
ATOM 15 C CB . CYS A 1 . ? 38.131 18.122 1.776 1.00 6.01 ? 507 CYS A CB 1
ATOM 16 S SG . CYS A 1 . ? 38.425 18.518 3.507 1.00 6.18 ? 507 CYS A SG 1
HETATM 27 CU CU . CU A 1 . ? 36.896 17.792 4.719 0.94 10.53 ? 554 CU A CU 1
";

/// Cys 57 of 1B4L in both copies of its biological assembly (two chains related by
/// a crystallographic two-fold), written with the crystal's cell and space group.
const CYS_TWO_COPIES: &str = "\
data_1b4l_cys57
_cell.length_a 118.900
_cell.length_b 118.900
_cell.length_c 75.200
_cell.angle_alpha 90.00
_cell.angle_beta 90.00
_cell.angle_gamma 120.00
_symmetry.space_group_name_H-M 'H 3 2'
";

const CYS_TWO_COPIES_ATOMS: &str = "\
ATOM 1 N N . CYS A1 1 . ? 0.589 13.997 29.886 1.00 17.42 ? 57 CYS A1 N 1
ATOM 2 C CA . CYS A1 1 . ? 1.614 14.728 29.150 1.00 16.50 ? 57 CYS A1 CA 1
ATOM 3 C C . CYS A1 1 . ? 1.022 15.962 28.481 1.00 16.52 ? 57 CYS A1 C 1
ATOM 4 O O . CYS A1 1 . ? 1.743 16.771 27.905 1.00 17.36 ? 57 CYS A1 O 1
ATOM 5 C CB . CYS A1 1 . ? 2.269 13.833 28.098 1.00 15.33 ? 57 CYS A1 CB 1
ATOM 6 S SG . CYS A1 1 . ? 3.168 12.443 28.833 1.00 14.05 ? 57 CYS A1 SG 1
ATOM 7 N N . CYS A2 1 . ? 11.827 7.509 45.314 1.00 17.42 ? 57 CYS A2 N 1
ATOM 8 C CA . CYS A2 1 . ? 11.948 8.762 46.050 1.00 16.50 ? 57 CYS A2 CA 1
ATOM 9 C C . CYS A2 1 . ? 13.312 8.866 46.719 1.00 16.52 ? 57 CYS A2 C 1
ATOM 10 O O . CYS A2 1 . ? 13.653 9.895 47.295 1.00 17.36 ? 57 CYS A2 O 1
ATOM 11 C CB . CYS A2 1 . ? 10.845 8.882 47.102 1.00 15.33 ? 57 CYS A2 CB 1
ATOM 12 S SG . CYS A2 1 . ? 9.192 8.965 46.367 1.00 14.05 ? 57 CYS A2 SG 1
";

/// Run fixed mode; the names of residue `resseq`'s atoms (None without chem_data).
fn atom_names(model: &str, resseq: &str) -> Option<Vec<String>> {
    let Some(root) = MonLib::locate(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return None;
    };
    let ml = MonLib::load(&root).unwrap();
    let out = pipeline::run(mmcif::read_mmcif(model).unwrap(), &ml, &Params::default()).unwrap();
    let mut names = Vec::new();
    for chain in &out.structure.models[0].chains {
        for rg in chain.residue_groups.iter().filter(|rg| rg.resseq.trim() == resseq) {
            for ag in &rg.atom_groups {
                names.extend(ag.atoms.iter().map(|a| a.name.trim().to_string()));
            }
        }
    }
    Some(names)
}

#[test]
fn a_cysteine_on_copper_is_a_thiolate() {
    let Some(names) = atom_names(&format!("data_cys_cu\n{HEADER}{CYS_CU}"), "507") else { return };
    assert!(names.iter().any(|n| n == "SG"), "{names:?}");
    assert!(!names.iter().any(|n| n == "HG"), "Cys 507 keeps HG: {names:?}");
}

#[test]
fn a_histidine_keeps_a_ring_hydrogen() {
    // 6TAE His 239 and everything within 6 A. Water-rich surroundings made the
    // state with neither ring hydrogen score best (13.00 against 10.76 for the
    // HIE that the neutron data show).
    let Some(names) = atom_names(include_str!("data/his_6tae_239.cif"), "239") else { return };
    assert!(names.iter().any(|n| n == "HD1" || n == "HE2"), "His 239 has no ring hydrogen: {names:?}");
}

#[test]
fn copies_in_the_model_are_not_symmetry_partners() {
    // The two-fold maps each copy's SG onto the other's, so the symmetry search
    // found a "disulfide" at 0 A for both and removed their HG; the SGs are 19 A apart.
    let model = format!("{CYS_TWO_COPIES}{HEADER}{CYS_TWO_COPIES_ATOMS}");
    let Some(names) = atom_names(&model, "57") else { return };
    assert_eq!(names.iter().filter(|n| *n == "HG").count(), 2, "{names:?}");
}
