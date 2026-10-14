use super::*;

fn matches(st: &Structure, flat: &FlatAtoms, i: u32, atom: &AtomLabel) -> bool {
    let p = flat.path[i as usize];
    flat.name[i as usize].trim() == atom.name.trim()
        && st.atom_group(p).resname.trim() == atom.resname.trim()
        && st.chain(p).id == atom.chain
        && st.residue_group(p).resseq.trim() == atom.resseq.trim()
        && st.residue_group(p).icode.trim() == atom.icode.trim()
        && (atom.altloc.is_empty() || flat.altloc[i as usize] == atom.altloc)
}

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

pub(super) fn add_replacement_angles(it: &mut Interp, st: &Structure, flat: &FlatAtoms) -> Vec<(u32, u32, String)> {
    let mut pairs = Vec::new();
    for link in &st.links {
        if link.kind != "covale" {
            continue;
        }
        let first: Vec<u32> = (0..flat.pos.len() as u32).filter(|&i| matches(st, flat, i, &link.atom1)).collect();
        let second: Vec<u32> = (0..flat.pos.len() as u32).filter(|&i| matches(st, flat, i, &link.atom2)).collect();
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
        if let Some(leaving) = replacement_angles(it, flat, parent, partner) {
            replacements.push((parent, partner, leaving));
        }
    }
    replacements
}

fn replacement_angles(it: &mut Interp, flat: &FlatAtoms, parent: u32, partner: u32) -> Option<String> {
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
        (&b.a1 == id || &b.a2 == id) && !matches!(b.type_.as_str(), "single" | "SING")
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
    if !terminal || residue(flat, parent) == residue(flat, partner) {
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
