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
