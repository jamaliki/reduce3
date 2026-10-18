use super::*;
use std::sync::Arc;

fn order(value: &str) -> u8 {
    if value.eq_ignore_ascii_case("single") || value.eq_ignore_ascii_case("SING") { 1 }
    else if value.eq_ignore_ascii_case("double") || value.eq_ignore_ascii_case("DOUB") { 2 }
    else if ["deloc", "DELO", "aromatic", "AROM"].iter().any(|v| value.eq_ignore_ascii_case(v)) { 3 }
    else { 0 }
}

/// Require complete source-defined trigonal carbon connectivity.
pub(super) fn dictionary_carbon(it: &Interp, at: &RidingAtoms, neighbors: &Adjacency, nb: &Nb) -> bool {
    let Some((parent, _)) = nb.a0 else { return false };
    let parent = parent as usize;
    let info = &it.dict[parent];
    let (Some(comp), Some(id)) = (&info.comp, &info.dict_id) else { return false };
    if comp.atom(id).is_none_or(|a| a.type_symbol != "C") || neighbors.of(parent as u32).len() != 3 {
        return false;
    }
    let mut expected = FxHashMap::default();
    for bond in &comp.bonds {
        let other = if &bond.a1 == id { &bond.a2 } else if &bond.a2 == id { &bond.a1 } else { continue };
        if expected.insert(other.as_str(), order(&bond.type_)).is_some() { return false; }
    }
    if expected.len() != 3 { return false; }
    let mut hydrogen = 0;
    let mut heavy = Vec::with_capacity(2);
    for &neighbor in neighbors.of(parent as u32) {
        let i = neighbor as usize;
        let info = &it.dict[i];
        let (Some(source), Some(name)) = (&info.comp, &info.dict_id) else { return false };
        if !Arc::ptr_eq(source, comp) || it.dict[i].instance != it.dict[parent].instance { return false; }
        let Some(bond_order) = expected.remove(name.as_str()) else { return false };
        if at.is_h[i] {
            if neighbor != nb.ih || bond_order != 1 || !comp.is_h(name) { return false; }
            hydrogen += 1;
        } else {
            if comp.is_h(name) { return false; }
            heavy.push(bond_order);
        }
    }
    heavy.sort_unstable();
    hydrogen == 1 && expected.is_empty() && matches!(heavy.as_slice(), [1, 2] | [3, 3])
}

#[cfg(test)]
mod tests;
