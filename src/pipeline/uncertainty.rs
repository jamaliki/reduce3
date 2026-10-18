//! Retain observed H when a declared parent lacks placement restraints.

use super::{Approach, Params};
use crate::hplace;
use crate::model::{Atom, AtomGroup, AtomPath, Structure};
use crate::monlib::{Comp, MonLib};
use crate::{names, resclass};
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;

mod bridges;
mod projection;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HydrogenUncertaintyReason {
    MissingHeavyNeighbor,
    InsufficientRestraints,
    MissingObservedBridge,
}

/// One removed, unplaced H from the placement pass.
#[derive(Clone, Debug)]
pub struct UnplacedHydrogen {
    pub parent_sources: Vec<u32>,
    pub hydrogen_name: String,
    pub label: String,
    /// The existing strict policy checks this site, even without an allowance.
    pub strict_restraints: bool,
    pub reason: HydrogenUncertaintyReason,
}

/// An undefined H target on one retained original heavy atom.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HydrogenUncertainty {
    pub parent_sources: Vec<u32>,
    pub hydrogen_name: String,
    pub reason: HydrogenUncertaintyReason,
}

struct Observed {
    parent_sources: Vec<u32>,
    dictionary_name: String,
    atom: Atom,
    comp: Arc<Comp>,
}

pub(super) struct Policy {
    allowed: FxHashSet<u32>,
    observed: Vec<Observed>,
    bridges: bridges::Bridges,
    allow_bond_targets: bool,
}

fn mapping(comp: &Comp, group: &AtomGroup) -> Vec<Option<String>> {
    let names = group.atoms.iter().map(|atom| atom.name.to_string()).collect::<Vec<_>>();
    let is_na = comp.test_for_rna_dna().is_some() || resclass::get_class(&group.resname).is_rna_dna();
    names::map_atoms(comp, &names, None, is_na, None).ids
}

fn parent(comp: &Comp, hydrogen: &str) -> Option<String> {
    let mut parents = comp.bonds.iter().filter_map(|bond| {
        let other = if bond.a1 == hydrogen { &bond.a2 }
            else if bond.a2 == hydrogen { &bond.a1 } else { return None };
        (!comp.is_h(other)).then_some(other.clone())
    }).collect::<Vec<_>>();
    parents.sort();
    parents.dedup();
    (parents.len() == 1).then(|| parents.remove(0))
}

fn terminal_parent(comp: &Comp, name: &str) -> Option<String> {
    if comp.has_atom(name) { return None; }
    if comp.test_for_peptide() && matches!(name, "H1" | "H2" | "H3") { return Some("N".into()); }
    hplace::terminal::hydrogen_parents(comp).into_iter().find_map(|(hydrogen, parent)| (hydrogen == name).then_some(parent))
}

fn observed_name(comp: &Comp, atom: &Atom, id: &Option<String>) -> Option<String> {
    if let Some(id) = id { return Some(id.clone()); }
    let name = if atom.elem() == "D" { atom.name.trim().replacen('D', "H", 1) }
        else { atom.name.trim().to_owned() };
    terminal_parent(comp, &name).map(|_| name)
}

impl Policy {
    pub fn read(st: &Structure, ml: &MonLib, params: &Params) -> Result<Self, String> {
        let rows = &params.allowed_unplaced_hydrogen_parent_sources;
        let allowed = rows.iter().copied().collect::<FxHashSet<_>>();
        if !rows.is_empty() || !params.allowed_unplaced_bridge_hydrogens.is_empty() || params.allow_unresolved_bond_targets {
            if params.compat || params.approach != Approach::Add || !params.stop_on_any_missing_hydrogen {
                return Err("Allowed unplaced hydrogen parents need strict fixed-mode hydrogen addition".into());
            }
            let heavy = st.atom_paths().into_iter().filter_map(|pth| {
                let atom = st.atom(pth);
                (!atom.is_hydrogen() && atom.src != Atom::NEW).then_some(atom.src)
            }).collect::<FxHashSet<_>>();
            if allowed.len() != rows.len() || !allowed.is_subset(&heavy) {
                return Err("Allowed unplaced hydrogen parents must be unique heavy source rows in the selected model".into());
            }
        }
        let mut result = Self { allowed, observed: Vec::new(), bridges: bridges::Bridges::read(st, ml, params)?,
            allow_bond_targets: params.allow_unresolved_bond_targets };
        if !result.enabled() { return Ok(result); }
        for group in st.models.iter().flat_map(|model| &model.chains)
            .flat_map(|chain| &chain.residue_groups).flat_map(|residue| &residue.atom_groups) {
            let selected = result.allow_bond_targets || group.atoms.iter().any(|atom| !atom.is_hydrogen()
                && (result.allows(atom.src) || result.bridges.contains_parent(atom.src)));
            if !selected || !group.atoms.iter().any(Atom::is_hydrogen) { continue; }
            result.observe(group, ml)?;
        }
        Ok(result)
    }

    pub fn allows(&self, row: u32) -> bool { self.allowed.contains(&row) }

    pub fn enabled(&self) -> bool { !self.allowed.is_empty() || !self.bridges.is_empty() || self.allow_bond_targets }

    pub fn authorize_bond_targets(&mut self, bonds: &[crate::interp::target_uncertainty::BondUncertainty]) {
        if self.allow_bond_targets {
            self.allowed.extend(bonds.iter().filter(|bond| bond.hydrogens_unknown).map(|bond| bond.parent_source));
        }
    }

    pub fn annotate_observed_bond_targets(&self, st: &Structure, placed: &hplace::Placed,
        reports: &mut Vec<crate::interp::target_uncertainty::BondUncertainty>) {
        if !self.allow_bond_targets { return; }
        let observed = self.observed.iter().filter_map(|site| {
            let [parent] = site.parent_sources.as_slice() else { return None };
            (site.dictionary_name == "HG").then_some(*parent)
        }).collect();
        crate::interp::target_uncertainty::add_observed_sulfur_targets(st, &placed.bonds, &observed, reports);
    }

    pub fn allows_bridge(&self, name: &str, rows: &[u32]) -> bool { self.bridges.allows(name, rows) }

    pub fn prepare_targets(&self, st: &mut Structure, placed: &mut hplace::Placed,
        bonds: &[crate::interp::target_uncertainty::BondUncertainty], ml: &MonLib) -> Result<Vec<HydrogenUncertainty>, String> {
        projection::prepare(self, st, placed, bonds, ml)
    }

    fn observe(&mut self, group: &AtomGroup, ml: &MonLib) -> Result<(), String> {
        let heavy_names = group.atoms.iter().filter(|atom| !atom.is_hydrogen())
            .map(|atom| atom.name.to_string()).collect::<Vec<_>>();
        let Some(comp) = hplace::hydrogen_dictionary(ml, &group.resname, &heavy_names, false) else { return Ok(()) };
        let ids = mapping(&comp, group);
        let heavy = group.atoms.iter().zip(&ids).filter_map(|(atom, id)| {
            if atom.is_hydrogen() { None } else { id.clone().map(|id| (id, atom.src)) }
        }).collect::<FxHashMap<_, _>>();
        for (atom, id) in group.atoms.iter().zip(&ids).filter(|(atom, _)| atom.is_hydrogen()) {
            let Some((dictionary_name, sources)) = self.bridges.observed_name(atom, group)
                .or_else(|| observed_name(&comp, atom, id)).and_then(|id| {
                let sources = self.bridges.observed_parents(&id, group).or_else(|| {
                    let parent = parent(&comp, &id).or_else(|| terminal_parent(&comp, &id))?;
                    Some(vec![heavy.get(&parent).copied()?])
                })?;
                Some((id, sources))
            }) else {
                return Err(format!("Observed hydrogen {} has no unique retained dictionary parent", atom.name.trim()));
            };
            if self.allow_bond_targets || sources.iter().any(|&source| self.allows(source)) || self.allows_bridge(&dictionary_name, &sources) {
                self.observed.push(Observed { parent_sources: sources, dictionary_name, atom: atom.clone(), comp: Arc::clone(&comp) });
            }
        }
        Ok(())
    }

    pub fn restore_observed(&self, st: &mut Structure, uncertainties: &[HydrogenUncertainty],
        bonds: &[crate::interp::target_uncertainty::BondUncertainty], ml: &MonLib) -> Result<(), String> {
        let affected = uncertainties.iter().flat_map(|site| site.parent_sources.iter().copied())
            .chain(bonds.iter().filter(|_| self.allow_bond_targets)
                .map(|bond| bond.parent_source)).collect::<FxHashSet<_>>();
        if affected.is_empty() { return Ok(()); }
        let parents = st.atom_paths().into_iter().filter_map(|pth| {
            let atom = st.atom(pth);
            (!atom.is_hydrogen() && affected.contains(&atom.src)).then_some((atom.src, pth))
        }).collect::<FxHashMap<_, AtomPath>>();
        for observed in self.observed.iter().filter(|observed| observed.parent_sources.len() == 1
            && observed.parent_sources.iter().all(|row| affected.contains(row))) {
            let pth = *parents.get(&observed.parent_sources[0]).ok_or("An uncertain hydrogen parent disappeared during reduction")?;
            let group = st.atom_group(pth);
            let heavy_names = group.atoms.iter().filter(|atom| !atom.is_hydrogen())
                .map(|atom| atom.name.to_string()).collect::<Vec<_>>();
            let comp = hplace::hydrogen_dictionary(ml, &group.resname, &heavy_names, false)
                .ok_or("An uncertain hydrogen parent lost its dictionary during reduction")?;
            let ids = mapping(&comp, group);
            let remove = group.atoms.iter().zip(ids).map(|(atom, id)| {
                atom.is_hydrogen() && (atom.src == observed.atom.src
                    || observed_name(&comp, atom, &id).as_deref() == Some(&observed.dictionary_name))
            }).collect::<Vec<_>>();
            let group = &mut st.models[pth.model as usize].chains[pth.chain as usize]
                .residue_groups[pth.rg as usize].atom_groups[pth.ag as usize];
            let mut index = 0;
            group.atoms.retain(|_| { let keep = !remove[index]; index += 1; keep });
            group.atoms.push(observed.atom.clone());
        }
        Ok(())
    }
}
