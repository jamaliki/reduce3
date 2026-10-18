//! Validate explicit FE/NI hydrides against their selected source dictionary.

use super::{Params, mapping};
use crate::cif;
use crate::hplace;
use crate::model::{Atom, AtomGroup, AtomPath, Structure};
use crate::monlib::{Comp, MonLib};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock, PoisonError};

#[derive(Default)]
pub(super) struct Bridges {
    sites: FxHashSet<(String, [u32; 2])>,
}

fn ordered(mut rows: [u32; 2]) -> [u32; 2] { rows.sort_unstable(); rows }

fn same_group(first: AtomPath, second: AtomPath) -> bool {
    (first.model, first.chain, first.rg, first.ag) == (second.model, second.chain, second.rg, second.ag)
}

fn metal(element: &str) -> bool { matches!(element, "FE" | "NI") }

fn dictionary_parents(comp: &Comp, name: &str) -> Option<[String; 2]> {
    if !comp.is_h(name) { return None; }
    let mut parents = Vec::new();
    for bond in &comp.bonds {
        let other = if bond.a1 == name { &bond.a2 }
            else if bond.a2 == name { &bond.a1 } else { continue };
        if !metal(&comp.atom(other)?.type_symbol)
            || !["single", "sing", "coval"].iter().any(|value| bond.type_.eq_ignore_ascii_case(value)) { return None; }
        parents.push(other.clone());
    }
    parents.sort();
    parents.dedup();
    parents.try_into().ok()
}

fn source_definition(comp: &Comp, name: &str, parents: &[String; 2], source: &str) -> bool {
    let document = cif::parse(source);
    let Some(block) = document.blocks.iter().find(|block| block.category("_chem_comp_atom").is_some_and(|atoms| {
        (0..atoms.nrows()).any(|row| atoms.get_tag(row, "comp_id") == Some(comp.id.as_str()))
    })) else { return false };
    let Some(atoms) = block.category("_chem_comp_atom") else { return false };
    let rows = (0..atoms.nrows()).filter(|&row| atoms.get_tag(row, "comp_id") == Some(comp.id.as_str())
        && atoms.get_tag(row, "atom_id") == Some(name)).collect::<Vec<_>>();
    let [row] = rows[..] else { return false };
    if atoms.get_tag(row, "charge").and_then(|value| value.parse::<i32>().ok()) != Some(-1)
        || !matches!(atoms.get_tag(row, "type_symbol"), Some("H" | "D")) { return false; }
    let Some(bonds) = block.category("_chem_comp_bond") else { return false };
    let mut neighbors = Vec::new();
    for row in 0..bonds.nrows() {
        if bonds.get_tag(row, "comp_id") != Some(comp.id.as_str()) { continue; }
        let (Some(first), Some(second)) = (bonds.get_tag(row, "atom_id_1"), bonds.get_tag(row, "atom_id_2")) else { return false };
        let other = if first == name { second } else if second == name { first } else { continue };
        if !bonds.get_tag(row, "value_order").is_some_and(|value| value.eq_ignore_ascii_case("sing"))
            || !bonds.get_tag(row, "pdbx_metal_coordination_flag").is_some_and(|value| value.eq_ignore_ascii_case("Y")) { return false; }
        neighbors.push(other.to_owned());
    }
    neighbors.sort();
    neighbors.dedup();
    neighbors == parents
}

type CacheKey = (PathBuf, String, Vec<String>, bool, String, [String; 2]);
static DEFINITIONS: OnceLock<Mutex<FxHashMap<CacheKey, bool>>> = OnceLock::new();

fn source_valid(comp: &Comp, name: &str, parents: &[String; 2]) -> Result<bool, String> {
    let key = (comp.source.clone(), comp.id.clone(), comp.applied_mods.clone(), comp.from_ccd, name.into(), parents.clone());
    let cache = DEFINITIONS.get_or_init(Mutex::default);
    if let Some(valid) = cache.lock().unwrap_or_else(PoisonError::into_inner).get(&key) { return Ok(*valid); }
    let source = std::fs::read_to_string(&comp.source).map_err(|error| format!("Cannot read bridge hydrogen dictionary {}: {error}", comp.id))?;
    let valid = source_definition(comp, name, parents, &source);
    cache.lock().unwrap_or_else(PoisonError::into_inner).insert(key, valid);
    Ok(valid)
}

impl Bridges {
    pub fn read(st: &Structure, ml: &MonLib, params: &Params) -> Result<Self, String> {
        let mut result = Self::default();
        for (name, sources) in &params.allowed_unplaced_bridge_hydrogens {
            let rows = ordered(*sources);
            if rows[0] == rows[1] || rows.contains(&Atom::NEW) || !result.sites.insert((name.clone(), rows)) {
                return Err("Allowed bridge hydrogens need a unique H name and two distinct original metal rows".into());
            }
            validate_site(st, ml, name, rows)?;
        }
        Ok(result)
    }

    pub fn is_empty(&self) -> bool { self.sites.is_empty() }

    pub fn sites(&self) -> impl Iterator<Item = (&str, &[u32; 2])> {
        self.sites.iter().map(|(name, rows)| (name.as_str(), rows))
    }

    pub fn contains_parent(&self, row: u32) -> bool { self.sites.iter().any(|(_, rows)| rows.contains(&row)) }

    pub fn allows(&self, name: &str, rows: &[u32]) -> bool {
        let [first, second] = rows else { return false };
        self.sites.contains(&(name.to_owned(), ordered([*first, *second])))
    }

    pub fn observed_parents(&self, name: &str, group: &AtomGroup) -> Option<Vec<u32>> {
        self.sites.iter().find_map(|(hydrogen, rows)| {
            (hydrogen == name && rows.iter().all(|row| group.atoms.iter().any(|atom| atom.src == *row)))
                .then(|| rows.to_vec())
        })
    }

    pub fn observed_name(&self, atom: &Atom, group: &AtomGroup) -> Option<String> {
        self.sites.iter().find_map(|(name, rows)| {
            (rows.iter().all(|row| group.atoms.iter().any(|parent| parent.src == *row))
                && (atom.name.trim() == name || (atom.elem() == "D"
                    && atom.name.trim() == name.replacen('H', "D", 1)))).then(|| name.clone())
        })
    }
}

fn source_paths(st: &Structure, rows: [u32; 2]) -> Result<[AtomPath; 2], String> {
    let mut paths = Vec::new();
    for row in rows {
        let found = st.atom_paths().into_iter().filter(|p| st.atom(*p).src == row).collect::<Vec<_>>();
        let [pth] = found[..] else { return Err("An allowed bridge metal row is absent or ambiguous in the selected model".into()) };
        let atom = st.atom(pth);
        if !metal(atom.elem()) || ![atom.xyz.x, atom.xyz.y, atom.xyz.z].iter().all(|value| value.is_finite()) {
            return Err("Allowed bridge hydrogen parents must have finite FE/NI coordinates".into());
        }
        paths.push(pth);
    }
    if !same_group(paths[0], paths[1]) { return Err("Allowed bridge hydrogen parents must share the same source atom group".into()); }
    Ok([paths[0], paths[1]])
}

fn validate_site(st: &Structure, ml: &MonLib, name: &str, rows: [u32; 2]) -> Result<(), String> {
    let paths = source_paths(st, rows)?;
    let group = st.atom_group(paths[0]);
    let isotopes = group.atoms.iter().filter(|atom| [name.to_owned(), name.replacen('H', "D", 1),
        name.replacen('H', "T", 1)].contains(&atom.name.trim().to_owned())).collect::<Vec<_>>();
    if isotopes.iter().any(|atom| atom.elem() == "T" || atom.name.trim() == name.replacen('H', "T", 1)) {
        return Err("Observed bridge tritium is not supported by hydrogen placement".into());
    }
    if isotopes.len() > 1 || isotopes.iter().any(|atom| !atom.is_hydrogen()
        || (atom.name.trim() != name && atom.elem() != "D")) {
        return Err("An observed bridge hydrogen has an ambiguous or mistyped isotope identity".into());
    }
    let names = group.atoms.iter().filter(|atom| !atom.is_hydrogen()).map(|atom| atom.name.to_string()).collect::<Vec<_>>();
    let comp = hplace::hydrogen_dictionary(ml, &group.resname, &names, false)
        .ok_or("An allowed bridge hydrogen has no selected dictionary")?;
    let parents = dictionary_parents(&comp, name).ok_or("An allowed bridge hydrogen lacks exactly two single-bond metal parents")?;
    if parents.iter().any(|name| group.atoms.iter().filter(|atom| !atom.is_hydrogen()
        && atom.name.trim() == name).count() != 1) {
        return Err("An allowed bridge dictionary metal name is absent or ambiguous in the source group".into());
    }
    let ids = mapping(&comp, group);
    let selected = group.atoms.iter().zip(ids).filter_map(|(atom, id)| {
        id.filter(|id| parents.contains(id) && comp.atom(id).is_some_and(|expected| expected.type_symbol == atom.elem()))
            .map(|_| atom.src)
    }).collect::<Vec<_>>();
    if selected.len() != 2 || ordered([selected[0], selected[1]]) != rows || !source_valid(&comp, name, &parents)? {
        return Err("An allowed bridge hydrogen lacks its exact negative-charge source metal bonds".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_source_metal_flags_negative_charge_and_exact_names_are_required() {
        let root = MonLib::locate(None).expect("Bridge tests require REDUCE3_CHEM_DATA");
        let library = MonLib::load(&root).unwrap();
        let comp = library.ccd_comp("EJ2", false).unwrap();
        let parents = dictionary_parents(&comp, "H").unwrap();
        let source = std::fs::read_to_string(&comp.source).unwrap();
        assert!(source_definition(&comp, "H", &parents, &source));
        for (charge, order, flag, hydrogen, second) in [
            ("0", "SING", "Y", "H", "NI"), ("?", "SING", "Y", "H", "NI"),
            ("-1", "DOUB", "Y", "H", "NI"), ("-1", "SING", "N", "H", "NI"),
            ("-1", "SING", "Y", "H2", "NI"), ("-1", "SING", "Y", "H", "C"),
        ] {
            let changed = format!(concat!(
                "data_test\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n",
                "_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nEJ2 {hydrogen} H {charge}\n",
                "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n",
                "_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_metal_coordination_flag\n",
                "EJ2 FE H SING Y\nEJ2 {second} H {order} {flag}\n",
            ), hydrogen = hydrogen, charge = charge, second = second, order = order, flag = flag);
            assert!(!source_definition(&comp, "H", &parents, &changed));
        }
    }
}
