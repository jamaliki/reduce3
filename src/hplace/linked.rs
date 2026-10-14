use super::*;

fn residue(flat: &FlatAtoms, i: u32) -> (u32, u32, u32) {
    let p = flat.path[i as usize];
    (p.model, p.chain, p.rg)
}

fn neighbors(it: &Interp, i: u32) -> Vec<u32> {
    it.bonds.iter().filter_map(|b| {
        if b.i == i { Some(b.j) } else if b.j == i { Some(b.i) } else { None }
    }).collect()
}

pub(super) fn tertiary_amide_h(it: &Interp, st: &Structure, flat: &FlatAtoms, h: u32) -> bool {
    if !flat.is_h[h as usize] {
        return false;
    }
    let parents: Vec<u32> = neighbors(it, h).into_iter().filter(|&i| !flat.is_h[i as usize]).collect();
    let [parent] = parents[..] else { return false };
    if flat.name[parent as usize].trim() != "N"
        || !resclass::get_class(st.atom_group(flat.path[parent as usize]).resname.trim()).is_amino_acid() {
        return false;
    }
    neighbors(it, parent).into_iter().filter(|&i| {
        !flat.is_h[i as usize] && (flat.altloc[i as usize].is_empty() || flat.altloc[i as usize] == flat.altloc[h as usize])
    }).count() >= 3
}

pub(super) fn add_replacement_angles(it: &mut Interp, st: &Structure, flat: &FlatAtoms,
    bond_order: &impl Fn(&Comp, &crate::monlib::CompBond, &str) -> Option<u8>) -> Vec<(u32, u32, String)> {
    if !st.links.iter().any(|link| link.kind == "covale") { return Vec::new(); }
    let resolver = EndpointResolver::new(st, &flat.path);
    let mut pairs = Vec::new();
    let mut annotated = FxHashSet::default();
    for link in &st.links {
        if link.kind != "covale" {
            continue;
        }
        let first: Vec<u32> = resolver.matching_indices(&link.atom1).collect();
        let second: Vec<u32> = resolver.matching_indices(&link.atom2).collect();
        annotated.extend(first.iter().chain(&second).copied());
        if link.order.as_deref().is_some_and(|v| !v.eq_ignore_ascii_case("sing") && !v.eq_ignore_ascii_case("single")) {
            continue;
        }
        let ([a], [b]) = (&first[..], &second[..]) else { continue };
        if flat.path[*a as usize].model == flat.path[*b as usize].model && it.has_bond(*a, *b) {
            pairs.extend([(*a, *b), (*b, *a)]);
        }
    }
    let mut replacements = Vec::new();
    for &(parent, partner) in &pairs {
        if pairs.iter().filter(|&&(p, _)| p == parent).count() != 1 {
            continue;
        }
        if let Some(leaving) = replacement_angles(it, st, flat, parent, partner, None, bond_order) {
            replacements.push((parent, partner, leaving));
        }
    }
    let templates: Vec<_> = replacements.iter().filter(|&&(p, q, _)| product_bond(it, st, flat, p, q, bond_order)).cloned().collect();
    let bonds: Vec<_> = it.bonds.iter().filter(|b| b.origin == interp::ORIGIN_MISC).map(|b| (b.i, b.j)).collect();
    for (a, b) in bonds {
        for (parent, partner) in [(a, b), (b, a)] {
            if annotated.contains(&parent) || !product_bond(it, st, flat, parent, partner, bond_order) { continue; }
            let matches: Vec<_> = templates.iter().filter(|&&(p, q, _)| same_site(it, parent, p) && same_site(it, partner, q)).collect();
            let Some(first) = matches.first() else { continue };
            if matches.iter().any(|t| t.2 != first.2) { continue; }
            if let Some(leaving) = replacement_angles(it, st, flat, parent, partner, Some(&first.2), bond_order) {
                replacements.push((parent, partner, leaving));
            }
        }
    }
    replacements
}

fn same_site(it: &Interp, a: u32, b: u32) -> bool {
    let (a, b) = (&it.dict[a as usize], &it.dict[b as usize]);
    let (Some(ac), Some(bc)) = (&a.comp, &b.comp) else { return false };
    a.dict_id.is_some() && a.dict_id == b.dict_id && ac.id == bc.id && ac.source == bc.source
        && ac.applied_mods == bc.applied_mods && ac.from_ccd == bc.from_ccd
}

fn product_bond(it: &Interp, st: &Structure, flat: &FlatAtoms, parent: u32, partner: u32,
    bond_order: &impl Fn(&Comp, &crate::monlib::CompBond, &str) -> Option<u8>) -> bool {
    let sulfur = &it.dict[partner as usize];
    let Some(comp) = &sulfur.comp else { return false };
    let heavy: Vec<_> = neighbors(it, partner).into_iter().filter(|&i| !flat.is_h[i as usize]).collect();
    let local: Vec<_> = heavy.iter().copied().filter(|&i| residue(flat, i) == residue(flat, partner)).collect();
    let [beta] = local[..] else { return false };
    let single = comp.bonds.iter().filter(|b| (b.a1 == "SG" && b.a2 == "CB") || (b.a1 == "CB" && b.a2 == "SG")).collect::<Vec<_>>();
    let [single] = single[..] else { return false };
    let parent_path = flat.path[parent as usize];
    let partner_path = flat.path[partner as usize];
    flat.element[parent as usize] == "C" && flat.element[partner as usize] == "S"
        && comp.id == "CYS" && sulfur.dict_id.as_deref() == Some("SG")
        && it.dict[beta as usize].dict_id.as_deref() == Some("CB")
        && bond_order(comp, single, &st.atom_group(flat.path[partner as usize]).resname) == Some(2)
        && parent_path.model == partner_path.model && residue(flat, parent) != residue(flat, partner)
        && (flat.altloc[parent as usize].is_empty() || flat.altloc[partner as usize].is_empty()
            || flat.altloc[parent as usize] == flat.altloc[partner as usize])
        && (1.5..=2.3).contains(&flat.pos[parent as usize].dist(flat.pos[partner as usize]))
        && heavy.len() == 2 && heavy.contains(&parent)
}

fn replacement_angles(it: &mut Interp, st: &Structure, flat: &FlatAtoms, parent: u32, partner: u32,
    required_leaving: Option<&str>, bond_order: &impl Fn(&Comp, &crate::monlib::CompBond, &str) -> Option<u8>) -> Option<String> {
    let info = &it.dict[parent as usize];
    let (Some(comp), Some(id)) = (&info.comp, &info.dict_id) else { return None };
    let Some(atom) = comp.atom(id) else { return None };
    let bonded = neighbors(it, parent);
    if atom.type_symbol != "C" || bonded.len() != 4
        || bonded.iter().filter(|&&i| !flat.is_h[i as usize]).count() != 2 {
        return None;
    }
    let expected: Vec<&str> = comp.bonds.iter().filter_map(|b| {
        let other = if &b.a1 == id { &b.a2 } else if &b.a2 == id { &b.a1 } else { return None };
        (!comp.is_h(other)).then_some(other.as_str())
    }).collect();
    if expected.len() != 2 || comp.bonds.iter().any(|b| {
        (&b.a1 == id || &b.a2 == id) && bond_order(comp, b, &st.atom_group(flat.path[parent as usize]).resname) != Some(2)
    }) {
        return None;
    }
    let local: FxHashMap<&str, u32> = bonded.iter().filter_map(|&i| {
        (residue(flat, i) == residue(flat, parent)).then(|| it.dict[i as usize].dict_id.as_deref().map(|id| (id, i))).flatten()
    }).collect();
    let missing: Vec<&str> = expected.into_iter().filter(|name| !local.contains_key(name)).collect();
    let [leaving] = missing[..] else { return None };
    let terminal = comp.bonds.iter().filter(|b| {
        (b.a1 == leaving && !comp.is_h(&b.a2)) || (b.a2 == leaving && !comp.is_h(&b.a1))
    }).count() == 1;
    if !terminal || residue(flat, parent) == residue(flat, partner)
        || required_leaving.is_some_and(|name| name != leaving)
        || !(1.2..=2.3).contains(&flat.pos[parent as usize].dist(flat.pos[partner as usize])) {
        return None;
    }
    let resolve = |name: &str| if name == leaving { Some(partner) } else { local.get(name).copied() };
    let angles: Vec<(u32, u32, f64, f64)> = comp.angles.iter().filter_map(|a| {
        if &a.a2 != id || (a.a1 != leaving && a.a3 != leaving) { return None }
        let (Some(i), Some(k), Some(v), Some(e)) = (resolve(&a.a1), resolve(&a.a3), a.value, a.esd) else { return None };
        (v.is_finite() && v > 0.0 && v < 180.0 && e.is_finite() && e > 0.0).then_some((i, k, v, e))
    }).collect();
    if angles.len() != 3 {
        return None;
    }
    let leaving = leaving.to_owned();
    for (i, k, value, sigma) in angles {
        it.push_angle_unchecked(i, parent, k, value, sigma);
    }
    Some(leaving)
}

pub(super) fn prochiral_atoms<'a>(st: &'a Structure, flat: &FlatAtoms, paths: &[AtomPath], replacements: &'a [(u32, u32, String)]) -> FxHashMap<&'a str, AtomPath> {
    let mut atoms: FxHashMap<&str, AtomPath> = paths.iter().map(|&p| (st.atom(p).name.trim(), p)).collect();
    for &(parent, partner, ref leaving) in replacements {
        if paths.contains(&flat.path[parent as usize]) {
            atoms.insert(leaving, flat.path[partner as usize]);
        }
    }
    atoms
}
