use reduce3::geom::v3;
use reduce3::interp::{target_uncertainty::{self, BondUncertainty, BondUncertaintyReason as Reason},
    FlatAtoms, ORIGIN_COVALENT, ORIGIN_MISC, ORIGIN_SS};
use reduce3::model::{Atom, AtomGroup, AtomLabel, AtomPath, LinkRecord, Structure};

fn fixture() -> Structure {
    reduce3::mmcif::read_mmcif(concat!(
        "data_targets\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "C CB . CYS A 1 0 0 0 1 10\nS SG . CYS A 1 1.8 0 0 1 10\n",
        "S SG . CYS A 2 3.83 0 0 1 10\nC CB . CYS A 2 5.63 0 0 1 10\n",
        "C C1 . LIG B 602 13.8 0 0 1 10\nC C2 . LIG B 602 20 0 0 1 10\n",
        "N N . GLY A 3 30 0 0 1 10\nC C . ALA A 4 31.33 0 0 1 10\n",
        "O OXT . ALA A 4 32.53 0 0 1 10\nO O . ALA A 4 31.33 1.2 0 1 10\n",
        "C CA . ALA A 4 31.33 0 1.5 1 10\nO OXT . ALA A 5 40 0 0 1 10\n",
    )).unwrap()
}

fn site(st: &Structure, seq: &str, name: &str) -> u32 {
    st.atom_paths().into_iter().position(|p| st.residue_group(p).resseq.trim() == seq
        && st.atom(p).name.trim() == name).unwrap() as u32
}

fn source(st: &Structure, index: u32) -> u32 { st.atom(st.atom_paths()[index as usize]).src }

fn group_mut(st: &mut Structure, p: AtomPath) -> &mut AtomGroup {
    &mut st.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize].atom_groups[p.ag as usize]
}

fn declare(st: &mut Structure, a: u32, b: u32, kind: &str) {
    let flat = FlatAtoms::from_structure(st);
    let label = |index: u32| {
        let p = flat.path[index as usize];
        AtomLabel { name: st.atom(p).name.clone(), resname: st.atom_group(p).resname.clone(),
            resseq: st.residue_group(p).resseq.clone(), chain: st.chain(p).id.clone(),
            altloc: st.atom_group(p).altloc.clone(), icode: st.residue_group(p).icode.clone() }
    };
    st.links.push(LinkRecord { atom1: label(a), atom2: label(b), kind: kind.into(), distance: None, order: None });
}

fn report(st: &Structure, index: u32, hydrogens_unknown: bool, reason: Reason) -> BondUncertainty {
    BondUncertainty { parent_source: source(st, index), hydrogens_unknown, reason }
}

#[test]
fn accepted_intrinsic_records_and_duplicate_disulfides_keep_known_targets() {
    let mut st = fixture();
    let (sg, cb, other) = (site(&st, "1", "SG"), site(&st, "1", "CB"), site(&st, "2", "SG"));
    declare(&mut st, sg, cb, "covale");
    declare(&mut st, sg, other, "disulf");
    let mut duplicate = st.links[1].clone();
    std::mem::swap(&mut duplicate.atom1, &mut duplicate.atom2);
    st.links.extend([st.links[1].clone(), duplicate]);
    assert!(target_uncertainty::inspect(&st, &[(sg, cb, ORIGIN_COVALENT), (sg, other, ORIGIN_SS)]).is_empty());
}

#[test]
fn ordinary_accepted_peptides_and_unrecorded_contacts_keep_known_targets() {
    let mut st = fixture();
    let (n, c) = (site(&st, "3", "N"), site(&st, "4", "C"));
    declare(&mut st, c, n, "covale");
    assert!(target_uncertainty::inspect(&st, &[(n, c, ORIGIN_COVALENT)]).is_empty());
    st.links.clear();
    assert!(target_uncertainty::inspect(&st, &[]).is_empty());
    declare(&mut st, c, n, "metalc");
    assert!(target_uncertainty::inspect(&st, &[]).is_empty());
}

#[test]
fn d_amino_acid_peptides_need_dictionary_or_named_peptide_provenance() {
    for origin in [ORIGIN_COVALENT, 35, 63, 65, 66, 70, 79, ORIGIN_MISC, 83] {
        let mut st = fixture();
        let (n, c) = (site(&st, "3", "N"), site(&st, "4", "C"));
        let p = st.atom_paths()[n as usize];
        group_mut(&mut st, p).resname = "DTY".into();
        declare(&mut st, c, n, "covale");
        let result = target_uncertainty::inspect(&st, &[(n, c, origin)]);
        if origin == ORIGIN_MISC || origin == 83 {
            assert_eq!(result, [report(&st, c, false, Reason::UnsupportedResidueEnvironment)]);
        } else { assert!(result.is_empty(), "origin {origin}"); }
    }
}

#[test]
fn multiple_external_sulfur_partners_mark_both_ends_of_every_edge() {
    for first_accepted in [false, true] {
        let mut st = fixture();
        let (sg, other, ligand, cb) = (site(&st, "1", "SG"), site(&st, "2", "SG"),
            site(&st, "602", "C1"), site(&st, "1", "CB"));
        declare(&mut st, sg, cb, "covale");
        declare(&mut st, sg, other, "disulf");
        declare(&mut st, sg, ligand, "covale");
        st.links.push(st.links[2].clone());
        let mut graph = vec![(sg, cb, ORIGIN_COVALENT), (sg, ligand, ORIGIN_MISC)];
        if first_accepted { graph.push((sg, other, ORIGIN_SS)); }
        let result = target_uncertainty::inspect(&st, &graph);
        assert_eq!(result, [report(&st, sg, true, Reason::MultipleSulfurPartners),
            report(&st, other, true, Reason::MultipleSulfurPartners),
            report(&st, ligand, true, Reason::MultipleSulfurPartners)]);
        st.links.reverse();
        assert_eq!(target_uncertainty::inspect(&st, &graph), result);
    }
}

#[test]
fn rejected_distant_links_report_source_rows_without_changing_atoms_or_bonds() {
    let mut st = fixture();
    let (sg, ligand) = (site(&st, "1", "SG"), site(&st, "602", "C1"));
    declare(&mut st, sg, ligand, "covale");
    let before = format!("{st:?}");
    let bonds = vec![(site(&st, "1", "CB"), sg, ORIGIN_COVALENT)];
    let saved = bonds.clone();
    assert_eq!(target_uncertainty::inspect(&st, &bonds), [
        report(&st, sg, true, Reason::UninterpretedBond), report(&st, ligand, true, Reason::UninterpretedBond)]);
    assert_eq!(format!("{st:?}"), before);
    assert_eq!(bonds, saved);
    assert!(st.atom(st.atom_paths()[sg as usize]).xyz.dist(st.atom(st.atom_paths()[ligand as usize]).xyz) > 10.0);
}

#[test]
fn accepted_custom_protein_edges_mask_chemical_state_without_masking_supported_h() {
    let mut st = fixture();
    let (c, ligand, oxt) = (site(&st, "4", "C"), site(&st, "602", "C1"), site(&st, "4", "OXT"));
    declare(&mut st, c, ligand, "covale");
    let result = target_uncertainty::inspect(&st, &[(c, ligand, ORIGIN_MISC), (c, oxt, ORIGIN_COVALENT)]);
    assert_eq!(result, [report(&st, c, false, Reason::UnsupportedResidueEnvironment),
        report(&st, oxt, false, Reason::UnsupportedResidueEnvironment)]);
    let sg = site(&st, "1", "SG");
    st.links.clear();
    declare(&mut st, sg, ligand, "covale");
    assert!(target_uncertainty::inspect(&st, &[(sg, ligand, ORIGIN_MISC)]).is_empty());
}

#[test]
fn unresolved_backbone_c_marks_only_its_graph_bonded_terminal_oxygen() {
    for origin in [ORIGIN_COVALENT, ORIGIN_MISC] {
        let mut st = fixture();
        let (c, ligand, oxt, oxygen) = (site(&st, "4", "C"), site(&st, "602", "C1"),
            site(&st, "4", "OXT"), site(&st, "4", "O"));
        declare(&mut st, c, ligand, "covale");
        let mut expected = vec![report(&st, ligand, true, Reason::UninterpretedBond),
            report(&st, c, true, Reason::UninterpretedBond)];
        if origin == ORIGIN_COVALENT { expected.push(report(&st, oxt, true, Reason::UninterpretedBond)); }
        expected.sort_by_key(|r| r.parent_source);
        assert_eq!(target_uncertainty::inspect(&st, &[(c, oxt, origin), (c, oxygen, ORIGIN_COVALENT)]), expected);
    }
}

#[test]
fn linked_sulfur_distinguishes_observed_and_generated_graph_attached_h() {
    for generated in [false, true] {
        for attached in [false, true] {
            let mut st = fixture();
            let sg_path = st.atom_paths()[site(&st, "1", "SG") as usize];
            let mut h = Atom::new("HG", "H", v3(1.8, 1.3, 0.));
            if !generated { h.src = 100; }
            group_mut(&mut st, sg_path).atoms.push(h);
            let (sg, other, h) = (site(&st, "1", "SG"), site(&st, "2", "SG"), site(&st, "1", "HG"));
            declare(&mut st, sg, other, "disulf");
            let mut graph = vec![(sg, other, ORIGIN_SS)];
            if attached { graph.push((sg, h, ORIGIN_COVALENT)); }
            let expected = if attached {
                vec![report(&st, sg, generated, Reason::UnsupportedResidueEnvironment)]
            } else { Vec::new() };
            assert_eq!(target_uncertainty::inspect(&st, &graph), expected);
        }
    }
}

#[test]
fn reports_use_original_rows_and_never_report_generated_heavy_atoms() {
    let mut st = fixture();
    let (sg, ligand) = (site(&st, "1", "SG"), site(&st, "602", "C1"));
    let paths = st.atom_paths();
    st.atom_mut(paths[sg as usize]).src = 7001;
    st.atom_mut(paths[ligand as usize]).src = Atom::NEW;
    declare(&mut st, sg, ligand, "covale");
    assert_eq!(target_uncertainty::inspect(&st, &[]), [BondUncertainty {
        parent_source: 7001, hydrogens_unknown: true, reason: Reason::UninterpretedBond,
    }]);
}

#[test]
fn self_links_and_cross_model_endpoints_remain_ambiguous() {
    for cross_model in [false, true] {
        let mut st = fixture();
        let (sg, other) = (site(&st, "1", "SG"), site(&st, "2", "SG"));
        declare(&mut st, sg, if cross_model { other } else { sg }, "disulf");
        let mut expected = vec![report(&st, sg, true, Reason::AmbiguousEndpoint)];
        if cross_model {
            expected.push(report(&st, other, true, Reason::AmbiguousEndpoint));
            let residue = st.models[0].chains[0].residue_groups.remove(1);
            let mut model = st.models[0].clone();
            model.id = "2".into();
            model.chains.truncate(1);
            model.chains[0].residue_groups = vec![residue];
            st.models.push(model);
        }
        assert_eq!(target_uncertainty::inspect(&st, &[]), expected);
    }
}

#[test]
fn exact_source_products_preserve_supported_h_targets() {
    let ml = reduce3::MonLib::load(&reduce3::MonLib::locate(None).expect("set REDUCE3_CHEM_DATA")).unwrap();
    for (name, cif) in [("9B5A", include_str!("data/covalent_9b5a.cif")), ("6WGN", include_str!("data/covalent_6wgn.cif")),
        ("8QSC", include_str!("data/covalent_8qsc.cif")), ("3PX7", include_str!("data/ptr_dt_3px7.cif"))] {
        let mut st = reduce3::mmcif::read_mmcif(cif).unwrap();
        let resolves = |label: &AtomLabel| st.atom_paths().into_iter().any(|p|
            st.atom(p).name.trim() == label.name.trim() && st.atom_group(p).resname.trim() == label.resname.trim()
                && st.chain(p).id == label.chain && st.residue_group(p).resseq.trim() == label.resseq.trim());
        let complete = st.links.iter().filter(|link| resolves(&link.atom1) && resolves(&link.atom2)).cloned().collect();
        st.links = complete;
        let placed = reduce3::hplace::place_hydrogens(&mut st, &ml, &reduce3::hplace::HPlaceParams {
            neutron: false, n_terminal_charge: reduce3::hplace::NTermCharge::ResidueOne,
            exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
        });
        let result = target_uncertainty::inspect(&st, &placed.bonds);
        if name == "6WGN" {
            assert_eq!(result, [report(&st, site(&st, "0", "CA"), false, Reason::UnsupportedResidueEnvironment)]);
        } else if name == "9B5A" {
            assert_eq!(result, [report(&st, site(&st, "75", "C"), false, Reason::UnsupportedResidueEnvironment)]);
        } else { assert!(result.is_empty(), "{name}: {result:?}"); }
    }
}

#[test]
fn stale_or_ambiguous_endpoints_keep_known_source_rows_and_a_fixed_reason() {
    for invalid in 0..4 {
        let mut st = fixture();
        let (sg, other, ligand) = (site(&st, "1", "SG"), site(&st, "2", "SG"), site(&st, "602", "C1"));
        declare(&mut st, sg, other, "disulf");
        declare(&mut st, sg, ligand, "covale");
        match invalid {
            0 => st.links[1].atom2.resseq = "999".into(),
            1 => st.links[1].atom2.chain = "missing".into(),
            2 => {
                let p = st.atom_paths()[ligand as usize];
                let mut duplicate = st.atom(p).clone();
                duplicate.src = 101;
                group_mut(&mut st, p).atoms.push(duplicate);
            }
            _ => {
                let p = st.atom_paths()[ligand as usize];
                st.atom_mut(p).src = Atom::NEW;
                st.links[1].atom2.resseq = "999".into();
            }
        }
        let result = target_uncertainty::inspect(&st, &[]);
        let sulfur = result.iter().find(|r| r.parent_source == source(&st, sg)).unwrap();
        assert_eq!(sulfur.reason, Reason::AmbiguousEndpoint);
        assert!(sulfur.hydrogens_unknown);
        assert!(result.iter().all(|r| r.parent_source != Atom::NEW));
        st.links.reverse();
        assert_eq!(target_uncertainty::inspect(&st, &[]), result);
    }
}
