//! Build weighted validation restraints without changing hydrogen preparation.

mod builder;
mod cdl;
mod charges;
mod chirality;
mod nomenclature;

pub use nomenclature::{SiteSwap, normalized_copy};

use super::*;

/// Parameters that affect chemical targets.
#[derive(Clone, Debug)]
pub struct ValidationOptions {
    pub nuclear: bool,
    pub use_cdl: bool,
    pub link_distance_cutoff: f64,
    pub chirality_sigma: f64,
    pub dictionary_links: Vec<DictionaryLink>,
}

/// Apply a declared link to the residues containing these flat atom indices.
#[derive(Clone, Debug)]
pub struct DictionaryLink {
    pub id: String,
    pub atom_indices: [u32; 2],
}

impl Default for ValidationOptions {
    fn default() -> Self {
        Self {
            nuclear: false,
            use_cdl: true,
            link_distance_cutoff: 3.0,
            chirality_sigma: 0.2,
            dictionary_links: Vec::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct WeightedBond {
    pub indices: [u32; 2],
    pub ideal: f64,
    pub weight: f64,
    pub origin: u16,
}

#[derive(Clone, Debug)]
pub struct WeightedSymmetryBond {
    pub indices: [u32; 2],
    pub image: crate::cell::SymOp,
    pub ideal: f64,
    pub weight: f64,
    pub origin: u16,
}

#[derive(Clone, Debug)]
pub struct WeightedAngle {
    pub indices: [u32; 3],
    pub ideal: f64,
    pub weight: f64,
    pub origin: u16,
}

#[derive(Clone, Debug)]
pub struct WeightedChirality {
    pub indices: [u32; 4],
    pub ideal: f64,
    pub weight: f64,
    pub both_signs: bool,
    pub pseudochiral: bool,
    pub origin: u16,
}

#[derive(Clone, Debug)]
pub struct WeightedPlane {
    pub indices: Vec<u32>,
    pub weights: Vec<f64>,
}

/// A check is complete only when no unsupported interpretation remains.
#[derive(Clone, Debug)]
pub struct Coverage {
    pub check: &'static str,
    pub unsupported: Vec<String>,
    pub unresolved: usize,
}

impl Coverage {
    pub fn complete(&self) -> bool {
        self.unsupported.is_empty() && self.unresolved == 0
    }
}

#[derive(Clone, Debug)]
pub struct DictionaryUse {
    pub residue: String,
    pub source: std::path::PathBuf,
    pub modifications: Vec<String>,
}

pub struct ValidationGraph {
    pub sites: Vec<Vec3>,
    pub atom_paths: Vec<AtomPath>,
    pub hydrogen_flags: Vec<bool>,
    /// Per-atom formal charges after reference scattering-label normalization.
    pub typed_charges: Vec<i32>,
    pub bonds: Vec<WeightedBond>,
    pub symmetry_bonds: Vec<WeightedSymmetryBond>,
    pub angles: Vec<WeightedAngle>,
    pub chiralities: Vec<WeightedChirality>,
    pub planes: Vec<WeightedPlane>,
    pub coverage: Vec<Coverage>,
    pub dictionaries: Vec<DictionaryUse>,
    pub topology: Interp,
    bond_keys: FxHashMap<[u32; 2], usize>,
    angle_keys: FxHashMap<[u32; 3], usize>,
    chiral_keys: FxHashMap<[u32; 4], usize>,
    plane_keys: FxHashMap<Vec<u32>, usize>,
}

fn bond_weight(ideal: f64, sigma: f64) -> Option<f64> {
    let weight = 1. / sigma.powi(2);
    (ideal.is_finite() && ideal > 0. && sigma.is_finite() && sigma > 0.
        && weight.is_finite() && weight > 0.).then_some(weight)
}

impl ValidationGraph {
    pub fn complete(&self) -> bool {
        self.coverage.iter().all(Coverage::complete)
    }

    fn image_bond(&mut self, mut indices: [u32;2], mut image: crate::cell::SymOp, ideal: f64, sigma: f64, origin: u16) {
        let Some(weight) = bond_weight(ideal,sigma) else {
            self.unsupported(&["bonds", "hydrogens"], "Missing or invalid symmetry bond target uncertainty");
            return;
        };
        let Some(reverse) = image.inverse() else {
            self.unsupported(builder::ALL, "Symmetry bond operator is singular");
            return;
        };
        if indices[0] > indices[1] { indices.swap(0,1); image = reverse; }
        self.symmetry_bonds.push(WeightedSymmetryBond { indices, image, ideal, weight, origin });
        // MCL inserts one directed ASU proxy. Disulfide discovery supplies its
        // two directed images separately.
        self.topology.sym_bonds.push(BondProxy { i:indices[0], j:indices[1], ideal, origin });
    }

    fn unsupported(&mut self, checks: &[&str], reason: impl Into<String>) {
        let reason = reason.into();
        for coverage in &mut self.coverage {
            if checks.contains(&coverage.check) && !coverage.unsupported.contains(&reason) {
                coverage.unsupported.push(reason.clone());
            }
        }
    }

    fn unresolved(&mut self, check: &str) {
        self.coverage
            .iter_mut()
            .find(|c| c.check == check)
            .unwrap()
            .unresolved += 1;
    }

    fn bond(&mut self, mut indices: [u32; 2], ideal: f64, sigma: f64, origin: u16) {
        indices.sort_unstable();
        let Some(weight) = bond_weight(ideal,sigma) else {
            self.unsupported(
                &["bonds", "hydrogens"],
                "Missing or invalid bond target uncertainty",
            );
            return;
        };
        if let Some(&k) = self.bond_keys.get(&indices) {
            if self.bonds[k].ideal != ideal || self.bonds[k].weight != weight {
                self.unsupported(
                    &["bonds"],
                    "Conflicting bond proxies require reference registry resolution",
                );
            }
            return;
        }
        self.bond_keys.insert(indices, self.bonds.len());
        self.bonds.push(WeightedBond {
            indices,
            ideal,
            weight,
            origin,
        });
        self.topology
            .upsert_bond(indices[0], indices[1], ideal, origin);
    }

    fn angle(&mut self, mut indices: [u32; 3], ideal: f64, sigma: f64) {
        if indices[0] > indices[2] {
            indices.swap(0, 2);
        }
        if !ideal.is_finite() || !sigma.is_finite() || sigma <= 0. {
            self.unsupported(&["angles"], "Missing or invalid angle target uncertainty");
            return;
        }
        let weight = 1. / sigma.powi(2);
        if let Some(&k) = self.angle_keys.get(&indices) {
            if self.angles[k].ideal != ideal || self.angles[k].weight != weight {
                self.unsupported(
                    &["angles"],
                    "Conflicting angle proxies require reference registry resolution",
                );
            }
            return;
        }
        self.angle_keys.insert(indices, self.angles.len());
        self.angles.push(WeightedAngle {
            indices,
            ideal,
            weight,
            origin: ORIGIN_COVALENT,
        });
    }

    fn angle_with_origin(&mut self, indices: [u32; 3], ideal: f64, sigma: f64, origin: u16) {
        let before = self.angles.len();
        self.angle(indices, ideal, sigma);
        for angle in &mut self.angles[before..] {
            angle.origin = origin;
        }
    }

    fn plane(&mut self, mut atoms: Vec<(u32, f64)>) {
        if atoms.len() < 4 {
            return;
        }
        atoms.sort_by_key(|a| a.0);
        let indices: Vec<u32> = atoms.iter().map(|a| a.0).collect();
        let weights: Vec<f64> = atoms.iter().map(|a| 1. / a.1.powi(2)).collect();
        if let Some(&k) = self.plane_keys.get(&indices) {
            if self.planes[k].weights != weights {
                self.unsupported(
                    &["planarities"],
                    "Conflicting plane weights need reference registry resolution",
                );
            }
            return;
        }
        self.plane_keys.insert(indices.clone(), self.planes.len());
        self.planes.push(WeightedPlane { indices, weights });
    }

    fn chirality(
        &mut self,
        mut indices: [u32; 4],
        mut ideal: f64,
        weight: f64,
        both_signs: bool,
        pseudochiral: bool,
    ) {
        for i in 1..4 {
            for j in i + 1..4 {
                if indices[i] > indices[j] {
                    indices.swap(i, j);
                    ideal = -ideal;
                }
            }
        }
        if !ideal.is_finite() || !weight.is_finite() || weight <= 0. {
            self.unsupported(&["chiralities"], "Invalid chiral target or uncertainty");
            return;
        }
        if let Some(&index) = self.chiral_keys.get(&indices) {
            let previous = &self.chiralities[index];
            if previous.ideal != ideal
                || previous.weight != weight
                || previous.both_signs != both_signs
                || previous.pseudochiral != pseudochiral
            {
                self.unsupported(
                    &["chiralities"],
                    "Conflicting chiral proxies need reference registry resolution",
                );
            }
            return;
        }
        self.chiral_keys.insert(indices, self.chiralities.len());
        self.chiralities.push(WeightedChirality {
            indices,
            ideal,
            weight,
            both_signs,
            pseudochiral,
            origin: ORIGIN_COVALENT,
        });
    }
}

/// Interpret fresh input sites and retain explicit coverage for every check.
pub fn build(
    st: &Structure,
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) -> Result<ValidationGraph, String> {
    builder::build(st, flat, ml, options)
}

#[cfg(test)]
mod target_tests {
    use super::bond_weight;

    #[test]
    fn bond_targets_need_finite_positive_weights() {
        assert_eq!(bond_weight(2.,0.5),Some(4.));
        for sigma in [0.,-0.1, f64::NAN, f64::INFINITY, 1e-250, 1e250] {
            assert!(bond_weight(2.,sigma).is_none());
        }
    }
}
