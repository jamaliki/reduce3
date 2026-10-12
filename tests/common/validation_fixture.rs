use reduce3::geom::v3;
use reduce3::model::{Atom, AtomGroup, Chain, Model, ResidueGroup, Structure};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Fixture(pub PathBuf);
impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

pub fn library() -> Fixture {
    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);
    let root = std::env::temp_dir().join(format!(
        "reduce3-validation-{}-{}",
        std::process::id(),
        NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(root.join("mon_lib/list")).unwrap();
    std::fs::create_dir_all(root.join("geostd/list")).unwrap();
    std::fs::create_dir_all(root.join("geostd/a")).unwrap();
    std::fs::write(root.join("geostd/list/geostd_list.cif"), "data_empty\n").unwrap();
    std::fs::write(root.join("geostd/ener_lib.cif"), "data_energy\n").unwrap();
    let list = r#"data_link_list
loop_
_chem_link.id
_chem_link.comp_id_1
_chem_link.mod_id_1
_chem_link.group_comp_1
_chem_link.comp_id_2
_chem_link.mod_id_2
_chem_link.group_comp_2
_chem_link.name
TRANS . CHANGE-CO peptide . CHANGE-NH peptide TRANS
data_link_TRANS
loop_
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
1 C 2 N 1.33 .02
data_mod_CHANGE-NH
loop_
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_value_dist
_chem_mod_bond.new_value_dist_esd
change N H .914 .010
data_mod_CHANGE-CO
loop_
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_value_dist
_chem_mod_bond.new_value_dist_esd
change C O 1.20 .010
"#;
    std::fs::write(root.join("mon_lib/list/mon_lib_list.cif"), list).unwrap();
    std::fs::write(
        root.join("geostd/list/mon_lib_list.cif"),
        list.replace(".914 .010", ".860 .020"),
    )
    .unwrap();
    std::fs::write(
        root.join("geostd/a/data_ALA.cif"),
        r#"data_comp_ALA
loop_
_chem_comp.id
_chem_comp.group
ALA peptide
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
N N NH1
CA C CH1
C C C
O O O
H H HNH1
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
N CA single 1.458 .019
CA C single 1.525 .021
C O double 1.230 .020
N H single .860 .020
"#,
    )
    .unwrap();
    Fixture(root)
}

pub fn one_residue(name: &str, atoms: &[(&str, &str)]) -> Structure {
    Structure {
        models: vec![Model {
            id: "1".into(),
            chains: vec![Chain {
                id: "A".into(),
                residue_groups: vec![ResidueGroup {
                    resseq: "1".into(),
                    icode: "".into(),
                    link_to_previous: false,
                    atom_groups: vec![AtomGroup {
                        altloc: "".into(),
                        resname: name.into(),
                        atoms: atoms
                            .iter()
                            .enumerate()
                            .map(|(i, (name, element))| {
                                Atom::new(name, element, v3(i as f64, 0., 0.))
                            })
                            .collect(),
                    }],
                }],
            }],
        }],
        ..Default::default()
    }
}

pub fn optional_library(
    fixture: &Fixture,
    code: &str,
    atoms: &[(&str, &str)],
    bonds: &str,
    angles: &str,
) {
    let directory = fixture
        .0
        .join(format!("geostd/{}", code[..1].to_ascii_lowercase()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut data = format!(
        "data_comp_{code}\nloop_\n_chem_comp.id\n_chem_comp.group\n{code} non-polymer\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.type_energy\n"
    );
    for (name, element) in atoms {
        data.push_str(&format!("{name} {element} {element}\n"));
    }
    data.push_str("loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.type\n_chem_comp_bond.value_dist\n_chem_comp_bond.value_dist_esd\n");
    data.push_str(bonds);
    data.push_str("loop_\n_chem_comp_angle.atom_id_1\n_chem_comp_angle.atom_id_2\n_chem_comp_angle.atom_id_3\n_chem_comp_angle.value_angle\n_chem_comp_angle.value_angle_esd\n");
    data.push_str(angles);
    std::fs::write(directory.join(format!("data_{code}.cif")), data).unwrap();
}

pub fn structure() -> Structure {
    let residues = (0..2)
        .map(|i| {
            let offset = i as f64 * 3.33;
            let atoms = [
                ("N", "N", 0.),
                ("CA", "C", 1.),
                ("C", "C", 2.),
                ("O", "O", 2.5),
                ("H", "H", -0.86),
            ]
            .iter()
            .map(|(name, element, x)| Atom::new(name, element, v3(x + offset, 0., 0.)))
            .collect();
            ResidueGroup {
                resseq: (i + 1).to_string(),
                icode: "".into(),
                link_to_previous: i > 0,
                atom_groups: vec![AtomGroup {
                    altloc: "".into(),
                    resname: "ALA".into(),
                    atoms,
                }],
            }
        })
        .collect();
    Structure {
        models: vec![Model {
            id: "1".into(),
            chains: vec![Chain {
                id: "A".into(),
                residue_groups: residues,
            }],
        }],
        ..Default::default()
    }
}
