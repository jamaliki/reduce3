//! Remove undefined H and map the existing optimizer graph exactly.

use super::{HydrogenUncertainty, HydrogenUncertaintyReason, Policy, mapping, observed_name};
use crate::hplace::Placed;
use crate::interp::{self, AtomDictInfo, EType};
use crate::model::{Atom, AtomGroup, AtomPath, Structure};
use crate::monlib::MonLib;
use crate::riding::RidingCoef;
use rustc_hash::{FxHashMap, FxHashSet};

fn group(st: &Structure, source: u32) -> Result<AtomPath, String> {
    st.atom_paths().into_iter().find(|p| st.atom(*p).src == source && !st.atom(*p).is_hydrogen())
        .ok_or_else(|| "A declared bridge parent disappeared during hydrogen placement".into())
}

fn mapped(map: &[u32], old: u32) -> Result<u32, String> {
    map.get(old as usize).copied().filter(|&index| index != u32::MAX)
        .ok_or_else(|| "A retained hydrogen restraint depends on a removed guessed hydride".into())
}

fn riding(mut coef: RidingCoef, map: &[u32]) -> Result<RidingCoef, String> {
    coef.ih = mapped(map, coef.ih)?;
    coef.a0 = mapped(map, coef.a0)?;
    coef.a1 = mapped(map, coef.a1)?;
    if coef.a2 >= 0 { coef.a2 = mapped(map, coef.a2 as u32)? as i64; }
    if coef.a3 >= 0 { coef.a3 = mapped(map, coef.a3 as u32)? as i64; }
    Ok(coef)
}

fn remove_sites(policy: &Policy, st: &mut Structure, placed: &Placed, unresolved: &FxHashSet<u32>,
    affected: &FxHashSet<u32>, ml: &MonLib) -> Result<(), String> {
    let mut remove = FxHashSet::default();
    let paths = st.atom_paths();
    for &(first, second, _) in &placed.bonds {
        for (hydrogen, parent) in [(first, second), (second, first)] {
            let parent = st.atom(paths[parent as usize]);
            if unresolved.contains(&parent.src) && st.atom(paths[hydrogen as usize]).is_hydrogen() {
                remove.insert(hydrogen as usize);
            }
        }
    }
    for (name, sources) in policy.bridges.sites() {
        let p = group(st, sources[0])?;
        let atoms = st.atom_group(p);
        let heavy = atoms.atoms.iter().filter(|atom| !atom.is_hydrogen()).map(|atom| atom.name.clone()).collect::<Vec<_>>();
        let comp = crate::hplace::hydrogen_dictionary(ml, &atoms.resname, &heavy, false)
            .ok_or("A declared bridge lost its selected dictionary")?;
        for atom in atoms.atoms.iter().filter(|atom| atom.is_hydrogen()) {
            let single = AtomGroup { atoms: vec![atom.clone()], resname: atoms.resname.clone(), altloc: atoms.altloc.clone() };
            if mapping(&comp, &single)[0].as_deref() == Some(name) { remove.insert(atom.i_seq); }
        }
    }
    for observed in policy.observed.iter().filter(|site| site.parent_sources.len() == 1
        && affected.contains(&site.parent_sources[0])) {
        let p = group(st, observed.parent_sources[0])?;
        let atoms = st.atom_group(p);
        for (atom, id) in atoms.atoms.iter().zip(mapping(&observed.comp, atoms)) {
            if atom.is_hydrogen() && (atom.src == observed.atom.src
                || observed_name(&observed.comp, atom, &id).as_deref() == Some(&observed.dictionary_name)) {
                remove.insert(atom.i_seq);
            }
        }
    }
    st.retain_atoms(|atom| !remove.contains(&atom.i_seq));
    Ok(())
}

fn restore_fixed(policy: &Policy, st: &mut Structure, unresolved: &FxHashSet<u32>) -> Result<(), String> {
    for observed in policy.observed.iter().filter(|observed| observed.parent_sources.len() == 2
        || observed.parent_sources.iter().any(|source| unresolved.contains(source))) {
        let p = group(st, observed.parent_sources[0])?;
        let mut atom = observed.atom.clone();
        atom.i_seq = usize::MAX;
        st.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
            .atom_groups[p.ag as usize].atoms.push(atom);
    }
    Ok(())
}

pub(super) fn prepare(policy: &Policy, st: &mut Structure, placed: &mut Placed,
    bonds: &[crate::interp::target_uncertainty::BondUncertainty], ml: &MonLib) -> Result<Vec<HydrogenUncertainty>, String> {
    let unresolved = bonds.iter().filter(|bond| policy.allow_bond_targets && bond.hydrogens_unknown)
        .map(|bond| bond.parent_source).collect::<FxHashSet<_>>();
    let affected = bonds.iter().filter(|_| policy.allow_bond_targets)
        .map(|bond| bond.parent_source).collect::<FxHashSet<_>>();
    if policy.bridges.is_empty() && affected.is_empty() { return Ok(Vec::new()); }
    st.reset_i_seq();
    let old_len = st.atoms_size();
    remove_sites(policy, st, placed, &unresolved, &affected, ml)?;
    restore_fixed(policy, st, &affected)?;
    project(policy, st, placed, old_len)?;
    let mut result = policy.bridges.sites().filter_map(|(name, sources)| {
        (!policy.observed.iter().any(|observed| observed.dictionary_name == name
            && observed.parent_sources == sources)).then(|| HydrogenUncertainty {
                parent_sources: sources.to_vec(), hydrogen_name: name.into(),
                reason: HydrogenUncertaintyReason::MissingObservedBridge,
            })
    }).collect::<Vec<_>>();
    result.sort_by(|first, second| first.parent_sources.cmp(&second.parent_sources).then_with(|| first.hydrogen_name.cmp(&second.hydrogen_name)));
    Ok(result)
}

fn project_arrays(policy: &Policy, st: &Structure, paths: &[AtomPath], placed: &mut Placed, map: &[u32]) -> Result<(), String> {
    let mut etype = Vec::with_capacity(paths.len());
    let mut dict = Vec::with_capacity(paths.len());
    let mut coefs = Vec::with_capacity(paths.len());
    for &p in paths {
        let old = st.atom(p).i_seq;
        if old == usize::MAX {
            let observed = policy.observed.iter().find(|observed| observed.atom.src == st.atom(p).src)
                .ok_or("A fixed hydrogen has no original source record")?;
            let energy = observed.comp.atom(&observed.dictionary_name).and_then(|atom| atom.type_energy.clone());
            etype.push(energy.map(EType::Typed).unwrap_or(EType::NoneType));
            dict.push(AtomDictInfo { comp: Some(std::sync::Arc::clone(&observed.comp)), dict_id: Some(observed.dictionary_name.clone()) });
            coefs.push(None);
        } else {
            etype.push(placed.etype[old].clone());
            dict.push(placed.dict[old].clone());
            coefs.push(placed.riding[old].map(|coef| riding(coef, map)).transpose()?);
        }
    }
    placed.etype = etype;
    placed.dict = dict;
    placed.riding = coefs;
    Ok(())
}

fn connect_observed(policy: &Policy, st: &Structure, paths: &[AtomPath], placed: &mut Placed) -> Result<(), String> {
    let source_indices = paths.iter().enumerate().filter_map(|(index, &p)| {
        (!st.atom(p).is_hydrogen()).then_some((st.atom(p).src, index as u32))
    }).collect::<FxHashMap<_, _>>();
    for (index, &p) in paths.iter().enumerate().filter(|(_, p)| st.atom(**p).i_seq == usize::MAX) {
        let atom = st.atom(p);
        if atom.src == Atom::NEW { return Err("A fixed hydrogen requires an observed original atom".into()); }
        let observed = policy.observed.iter().find(|observed| observed.atom.src == atom.src)
            .ok_or("A fixed hydrogen lacks its validated parent rows")?;
        let origin = if observed.parent_sources.len() == 2 { interp::ORIGIN_METAL } else { interp::ORIGIN_COVALENT };
        for source in &observed.parent_sources {
            placed.bonds.push((index as u32, *source_indices.get(source).ok_or("A fixed hydrogen parent is absent from its graph")?, origin));
        }
    }
    Ok(())
}

fn project(policy: &Policy, st: &mut Structure, placed: &mut Placed, old_len: usize) -> Result<(), String> {
    let paths = st.atom_paths();
    let mut old_to_new = vec![u32::MAX; old_len];
    for (index, &p) in paths.iter().enumerate() {
        let old = st.atom(p).i_seq;
        if old < old_len { old_to_new[old] = index as u32; }
    }
    project_arrays(policy, st, &paths, placed, &old_to_new)?;
    let mut bonds = Vec::with_capacity(placed.bonds.len());
    for &(first, second, origin) in &placed.bonds {
        let first = *old_to_new.get(first as usize).ok_or("A placed bond is outside the original atom map")?;
        let second = *old_to_new.get(second as usize).ok_or("A placed bond is outside the original atom map")?;
        if first != u32::MAX && second != u32::MAX { bonds.push((first, second, origin)); }
    }
    placed.bonds = bonds;
    connect_observed(policy, st, &paths, placed)?;
    placed.n_h_final = paths.iter().filter(|&&p| st.atom(p).is_hydrogen()).count();
    st.reset_i_seq();
    Ok(())
}
