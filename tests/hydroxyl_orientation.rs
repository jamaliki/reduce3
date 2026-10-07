//! Fixed mode's hydroxyl preferences: with nothing to hydrogen-bond to, a
//! tyrosine's hydroxyl hydrogen lies in the ring plane and a carboxylic acid's
//! hydrogen is syn to its carbonyl oxygen. Needs the monomer library (found as
//! the command-line tool finds it); skipped without one.

use reduce3::geom::{dihedral_deg, Vec3};
use reduce3::{mmcif, pipeline, MonLib, Params};

/// Tyr 7 and the 4-mercuribenzoic acid (MBO) of 2QO8, about 25 A apart.
const MODEL: &str = "\
data_test
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
ATOM   45 N  N   . TYR A 1 6 ? 4.981  -6.717  10.785 1.00 12.74 ? 7   TYR A N   1
ATOM   46 C  CA  . TYR A 1 6 ? 4.621  -7.832  11.650 1.00 13.75 ? 7   TYR A CA  1
ATOM   47 C  C   . TYR A 1 6 ? 5.516  -9.043  11.418 1.00 16.26 ? 7   TYR A C   1
ATOM   48 O  O   . TYR A 1 6 ? 5.211  -10.147 11.865 1.00 17.99 ? 7   TYR A O   1
ATOM   49 C  CB  . TYR A 1 6 ? 3.154  -8.214  11.434 1.00 12.91 ? 7   TYR A CB  1
ATOM   50 C  CG  . TYR A 1 6 ? 2.185  -7.141  11.878 1.00 10.04 ? 7   TYR A CG  1
ATOM   51 C  CD1 . TYR A 1 6 ? 1.745  -6.156  10.994 1.00 9.16  ? 7   TYR A CD1 1
ATOM   52 C  CD2 . TYR A 1 6 ? 1.731  -7.094  13.195 1.00 9.58  ? 7   TYR A CD2 1
ATOM   53 C  CE1 . TYR A 1 6 ? 0.875  -5.148  11.410 1.00 7.70  ? 7   TYR A CE1 1
ATOM   54 C  CE2 . TYR A 1 6 ? 0.866  -6.091  13.623 1.00 8.82  ? 7   TYR A CE2 1
ATOM   55 C  CZ  . TYR A 1 6 ? 0.441  -5.122  12.728 1.00 7.68  ? 7   TYR A CZ  1
ATOM   56 O  OH  . TYR A 1 6 ? -0.415 -4.135  13.154 1.00 8.87  ? 7   TYR A OH  1
HETATM 2118 HG HG  . MBO C 3 . ? -1.917 10.664 24.383 1.00 18.04 ? 263 MBO B HG  1
HETATM 2119 C  CE1 . MBO C 3 . ? -1.257 11.366 26.177 1.00 16.65 ? 263 MBO B CE1 1
HETATM 2120 C  CE2 . MBO C 3 . ? -0.828 12.707 26.298 1.00 18.38 ? 263 MBO B CE2 1
HETATM 2121 C  CE3 . MBO C 3 . ? -0.342 13.190 27.535 1.00 18.98 ? 263 MBO B CE3 1
HETATM 2122 C  CE4 . MBO C 3 . ? -0.276 12.333 28.665 1.00 18.43 ? 263 MBO B CE4 1
HETATM 2123 C  CE5 . MBO C 3 . ? -0.710 10.991 28.535 1.00 18.98 ? 263 MBO B CE5 1
HETATM 2124 C  CE6 . MBO C 3 . ? -1.196 10.510 27.297 1.00 18.79 ? 263 MBO B CE6 1
HETATM 2125 C  CZ  . MBO C 3 . ? 0.221  12.811 29.914 1.00 18.80 ? 263 MBO B CZ  1
HETATM 2126 O  OZ1 . MBO C 3 . ? 0.601  13.972 30.027 1.00 17.85 ? 263 MBO B OZ1 1
HETATM 2127 O  OZ2 . MBO C 3 . ? 0.280  12.067 30.897 1.00 19.62 ? 263 MBO B OZ2 1
";

fn run(params: &Params) -> Option<Vec<(String, Vec3)>> {
    let Some(root) = MonLib::locate(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return None;
    };
    let ml = MonLib::load(&root).unwrap();
    let out = pipeline::run(mmcif::read_mmcif(MODEL).unwrap(), &ml, params).unwrap();
    let mut atoms = Vec::new();
    for chain in &out.structure.models[0].chains {
        for rg in &chain.residue_groups {
            for ag in &rg.atom_groups {
                atoms.extend(ag.atoms.iter().map(|a| (format!("{} {}", ag.resname, a.name.trim()), a.xyz)));
            }
        }
    }
    Some(atoms)
}

fn dihedral(atoms: &[(String, Vec3)], names: [&str; 4]) -> f64 {
    let at = |n: &str| atoms.iter().find(|(k, _)| k == n).unwrap_or_else(|| panic!("no atom {n}")).1;
    dihedral_deg(at(names[0]), at(names[1]), at(names[2]), at(names[3])).unwrap()
}

#[test]
fn hydroxyls_prefer_the_plane_and_acids_syn() {
    let Some(atoms) = run(&Params::default()) else { return };
    // the phenol hydrogen in the ring plane, on either side
    let tyr = dihedral(&atoms, ["TYR CE1", "TYR CZ", "TYR OH", "TYR HH"]).abs();
    assert!(tyr < 10.0 || tyr > 170.0, "TYR HH {tyr} degrees from CE1");
    // the acid hydrogen syn to the carbonyl oxygen
    let acid = dihedral(&atoms, ["MBO OZ1", "MBO CZ", "MBO OZ2", "MBO HZ2"]).abs();
    assert!(acid < 10.0, "MBO HZ2 {acid} degrees from syn");
}

#[test]
fn without_the_preferences_the_acid_is_not_syn() {
    // Reduce2's rotator has no preference: alone, this acid's hydrogen ends
    // up on the anti side (about 150 degrees), so the syn result above comes
    // from the preference terms.
    let mut params = Params::default();
    params.opt.acid_syn_preference = 0.0;
    params.opt.planar_hydroxyl_preference = 0.0;
    let Some(atoms) = run(&params) else { return };
    let acid = dihedral(&atoms, ["MBO OZ1", "MBO CZ", "MBO OZ2", "MBO HZ2"]).abs();
    assert!(acid > 90.0, "MBO HZ2 {acid} degrees from syn");
}
