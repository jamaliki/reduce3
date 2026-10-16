use super::*;
use crate::monlib::CompBond;

fn order(value: &str) -> Option<u8> {
    if ["single", "sing"].iter().any(|v| value.eq_ignore_ascii_case(v)) { Some(2) }
    else if ["double", "doub"].iter().any(|v| value.eq_ignore_ascii_case(v)) { Some(4) }
    else if ["triple", "trip"].iter().any(|v| value.eq_ignore_ascii_case(v)) { Some(6) }
    else if ["deloc", "delo", "aromatic", "arom"].iter().any(|v| value.eq_ignore_ascii_case(v)) { Some(3) }
    else { None }
}

pub(crate) fn source_order(ml: &MonLib, comp: &Comp, bond: &CompBond) -> Option<u8> {
    if !bond.type_.eq_ignore_ascii_case("coval") { return order(&bond.type_); }
    let ccd = ml.ccd(&comp.id)?;
    let found: Vec<_> = ccd.bonds.iter().filter(|(a, b, _)|
        (a == &bond.a1 && b == &bond.a2) || (a == &bond.a2 && b == &bond.a1)).collect();
    let [(_, _, value)] = found[..] else { return None };
    order(value)
}

fn heavy_valence(it: &Interp, flat: &FlatAtoms, parent: u32, bond_order: &impl Fn(&Comp, &CompBond) -> Option<u8>) -> Option<u8> {
    let info = &it.dict[parent as usize];
    let (Some(comp), Some(id)) = (&info.comp, &info.dict_id) else { return None };
    let mut value = 0u8;
    for bond in &it.bonds {
        let other = if bond.i == parent { bond.j } else if bond.j == parent { bond.i } else { continue };
        if flat.is_h[other as usize] { continue; }
        let same = flat.path[parent as usize].model == flat.path[other as usize].model
            && flat.path[parent as usize].chain == flat.path[other as usize].chain
            && flat.path[parent as usize].rg == flat.path[other as usize].rg;
        let weight = if same {
            let other_id = it.dict[other as usize].dict_id.as_ref()?;
            let source: Vec<_> = comp.bonds.iter().filter(|b| (&b.a1 == id && &b.a2 == other_id) || (&b.a2 == id && &b.a1 == other_id)).collect();
            let [source] = source[..] else { return None };
            bond_order(comp, source)?
        } else if bond.origin != ORIGIN_METAL { 2 } else { continue };
        value = value.checked_add(weight)?;
    }
    Some(value)
}

/// Honor declared thioethers before automatic linking can consume an endpoint.
pub(super) fn add_thioethers(it: &mut Interp, st: &Structure, flat: &FlatAtoms, bond_order: &impl Fn(&Comp, &CompBond) -> Option<u8>) {
    let mut pairs = Vec::new();
    let mut seen = FxHashSet::default();
    for link in &st.links {
        if link.kind != "covale" || link.order.as_deref().is_some_and(|v| order(v) != Some(2)) { continue; }
        let sulfur = |a: &AtomLabel| a.resname.trim() == "CYS" && a.name.trim() == "SG";
        let (s, c) = if sulfur(&link.atom1) { (&link.atom1, &link.atom2) }
            else if sulfur(&link.atom2) { (&link.atom2, &link.atom1) } else { continue };
        let candidates = |a: &AtomLabel| (0..flat.pos.len() as u32).filter(|&i| a.matches(st, flat.path[i as usize])).collect::<Vec<_>>();
        let (first, second) = (candidates(s), candidates(c));
        let ([s], [c]) = (&first[..], &second[..]) else { continue };
        let (sp, cp) = (flat.path[*s as usize], flat.path[*c as usize]);
        let distance = flat.pos[*s as usize].dist(flat.pos[*c as usize]);
        if sp.model != cp.model || (sp.chain == cp.chain && sp.rg == cp.rg)
            || flat.element[*s as usize] != "S" || flat.element[*c as usize] != "C"
            || !(1.5..=2.3).contains(&distance)
            || (!flat.altloc[*s as usize].is_empty() && !flat.altloc[*c as usize].is_empty()
                && flat.altloc[*s as usize] != flat.altloc[*c as usize]) { continue; }
        if seen.insert((*s, *c)) { pairs.push((*s, *c, distance)); }
    }
    let allowed: Vec<_> = pairs.iter().copied().filter(|&(s, c, _)| {
        let sulfur_links = pairs.iter().filter(|&&(a, b, _)| a == s && !it.has_bond(a, b)).count();
        let carbon_links = pairs.iter().filter(|&&(a, b, _)| b == c && !it.has_bond(a, b)).count();
        sulfur_links <= 1 && heavy_valence(it, flat, s, bond_order).is_some_and(|v| usize::from(v) + 2 * sulfur_links <= 4)
            && heavy_valence(it, flat, c, bond_order).is_some_and(|v| usize::from(v) + 2 * carbon_links <= 8)
    }).collect();
    for (s, c, distance) in allowed { it.add_bond(s, c, distance, ORIGIN_MISC); }
}

#[cfg(test)]
mod tests;
