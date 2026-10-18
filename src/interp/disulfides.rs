use super::*;

fn matches(st: &Structure, flat: &FlatAtoms, i: u32, label: &AtomLabel) -> bool {
    label.matches(st, flat.path[i as usize])
}

/// Keep validated deposited disulfides even when one sulfur is near a metal.
pub(super) fn add_explicit(it: &mut Interp, st: &Structure, flat: &FlatAtoms) -> FxHashSet<u32> {
    let mut pairs = Vec::new();
    for link in st.links.iter().filter(|l| l.kind == "disulf") {
        let candidates = |label: &AtomLabel| (0..flat.pos.len() as u32).filter(|&i| matches(st, flat, i, label)).collect::<Vec<_>>();
        let (first, second) = (candidates(&link.atom1), candidates(&link.atom2));
        let ([a], [b]) = (&first[..], &second[..]) else { continue };
        if a == b || flat.path[*a as usize].model != flat.path[*b as usize].model {
            continue;
        }
        let sulfur = |i: u32| flat.element[i as usize] == "S" && flat.name[i as usize].trim() == "SG"
            && st.atom_group(flat.path[i as usize]).resname.trim() == "CYS";
        let (alt_a, alt_b) = (&flat.altloc[*a as usize], &flat.altloc[*b as usize]);
        let distance = flat.pos[*a as usize].dist(flat.pos[*b as usize]);
        if sulfur(*a) && sulfur(*b) && distance.is_finite() && distance > 1.5 && distance <= 3.0
            && (alt_a.is_empty() || alt_b.is_empty() || alt_a == alt_b) {
            pairs.push((*a, *b));
        }
    }
    let mut explicit = FxHashSet::default();
    for &(a, b) in &pairs {
        if [a, b].iter().any(|i| pairs.iter().filter(|(x, y)| x == i || y == i).count() != 1) {
            continue;
        }
        it.add_bond(a, b, 2.031, ORIGIN_SS);
        explicit.extend([a, b]);
    }
    explicit
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (Structure, FlatAtoms, Interp) {
        let mut st = crate::mmcif::read_mmcif(concat!(
            "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n",
            "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
            "S SG . CYS A 188 0 0 0 1 10\nS SG . CYS A 192 2.076 0 0 1 10\nFE FE4 . SF4 A 1007 2.076 2.31 0 1 10\n",
        )).unwrap();
        let atom = |seq: &str| AtomLabel { name: "SG".into(), resname: "CYS".into(), chain: "A".into(), resseq: seq.into(), altloc: String::new(), icode: String::new() };
        st.links.push(LinkRecord { atom1: atom("188"), atom2: atom("192"), kind: "disulf".into(), order: None, distance: Some(2.076) });
        let flat = FlatAtoms::from_structure(&st);
        (st, flat, Interp::empty(3))
    }

    #[test]
    fn deposited_disulfide_overrides_metal_exclusion_in_fixed_mode() {
        for compat in [false, true] {
            let (st, flat, mut it) = fixture();
            super::super::add_disulfides(&mut it, &st, &flat, &super::super::SpecialAtoms::new(&flat), &[0, 1], compat);
            assert_eq!(it.has_bond(0, 1), !compat);
            assert_eq!(flat.pos[0], st.atom(flat.path[0]).xyz);
        }
    }

    #[test]
    fn stale_ambiguous_or_implausible_deposited_links_do_not_add_bonds() {
        for invalid in 0..6 {
            let (mut st, mut flat, mut it) = fixture();
            match invalid {
                0 => st.links[0].atom1.chain = "B".into(),
                1 => st.links.push(st.links[0].clone()),
                2 => flat.pos[1] = Vec3::from_array([8.0, 0.0, 0.0]),
                3 => st.links[0].atom1.icode = "A".into(),
                4 => flat.element[1] = "O".into(),
                5 => { flat.altloc[0] = "A".into(); flat.altloc[1] = "B".into(); }
                _ => unreachable!(),
            }
            assert!(add_explicit(&mut it, &st, &flat).is_empty());
            assert!(!it.has_bond(0, 1));
        }
    }
}
