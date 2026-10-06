//! Restraints for residues that only the wwPDB chemical component dictionary
//! (CCD) describes: port of Reduce2's fallback in `mon_lib_query`,
//! `mmtbx.hydrogens.reduce_hydrogen.get_h_restraints(resname, strict=False)`.
//!
//! Reduce2 first builds an RDKit molecule from the CCD entry
//! (`mmtbx.ligands.rdkit_utils.read_chemical_component_filename`). That fails,
//! and the residue is left without restraints, for bond orders other than
//! SING/DOUB/TRIP, elements or charges RDKit cannot read, misaligned or
//! missing coordinates, and valences RDKit's sanitization rejects
//! ([`crate::rdkit_valence`] holds RDKit's own verdicts). The dictionary then
//! has every CCD atom without an energy type, every bond at 0.9 times its
//! length in the ideal coordinates, and every angle and torsion at its value
//! there. Angles and torsions come out in the iteration order of the Python
//! sets RDKit's neighbours are collected in; riding H placement depends on the
//! torsion order, so [`PySet`] reproduces CPython's set layout.

use crate::cif;
use crate::monlib::{Comp, CompAngle, CompAtom, CompBond, CompTor};
use crate::rdkit_valence::{ACCEPTED, MIN_CHARGE};
use std::path::Path;

type P3 = [f64; 3];

fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
// RDKit's Point3D arithmetic. The arm64 builds that Reduce2 runs with fuse
// multiply-adds as below (read from the disassembly of
// MolTransforms::getDihedralRad); the sign of an exactly planar torsion
// (+-180) depends on it.
const FUSED: bool = cfg!(target_arch = "aarch64");

fn dot(a: P3, b: P3) -> f64 {
    if FUSED {
        a[2].mul_add(b[2], a[0].mul_add(b[0], a[1] * b[1]))
    } else {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
}
fn cross(a: P3, b: P3) -> P3 {
    if FUSED {
        [
            a[1].mul_add(b[2], -(a[2] * b[1])),
            (-a[0]).mul_add(b[2], a[2] * b[0]),
            a[0].mul_add(b[1], -(a[1] * b[0])),
        ]
    } else {
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    }
}

/// `rdMolTransforms.GetBondLength`.
fn bond_length(p: &[P3], i: usize, j: usize) -> f64 {
    let d = sub(p[i], p[j]);
    dot(d, d).sqrt()
}

/// `rdMolTransforms.GetAngleDeg`; `None` where RDKit raises (coincident atoms).
fn angle_deg(p: &[P3], i: usize, j: usize, k: usize) -> Option<f64> {
    let rji = sub(p[i], p[j]);
    let rjk = sub(p[k], p[j]);
    if dot(rji, rji) <= 1e-16 || dot(rjk, rjk) <= 1e-16 {
        return None;
    }
    let lsq = dot(rji, rji) * dot(rjk, rjk);
    let d = dot(rji, rjk) / lsq.sqrt();
    let rad = if d <= -1.0 {
        std::f64::consts::PI
    } else if d >= 1.0 {
        0.0
    } else {
        d.acos()
    };
    Some(180.0 * rad / std::f64::consts::PI)
}

/// `rdMolTransforms.GetDihedralDeg`; `None` where RDKit raises.
fn dihedral_deg(p: &[P3], i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
    let rij = sub(p[j], p[i]);
    let rjk = sub(p[k], p[j]);
    let rkl = sub(p[l], p[k]);
    if dot(rij, rij) <= 1e-16 || dot(rjk, rjk) <= 1e-16 || dot(rkl, rkl) <= 1e-16 {
        return None;
    }
    let nijk = cross(rij, rjk);
    let njkl = cross(rjk, rkl);
    let m = cross(nijk, rjk);
    let rad = -(dot(m, njkl) / (dot(njkl, njkl) * dot(m, m)).sqrt())
        .atan2(dot(nijk, njkl) / (dot(nijk, nijk) * dot(njkl, njkl)).sqrt());
    Some(180.0 * rad / std::f64::consts::PI)
}

/// A CPython `set`'s slot layout, for reproducing its iteration order
/// (Objects/setobject.c: linear probing over 9 slots, then perturbation;
/// resize to 4 times the used count once 3/5 full).
pub struct PySet<K> {
    table: Vec<Option<(u64, K)>>,
    fill: usize,
}

const LINEAR_PROBES: usize = 9;

impl<K: PartialEq> PySet<K> {
    pub fn new() -> PySet<K> {
        PySet { table: (0..8).map(|_| None).collect(), fill: 0 }
    }

    pub fn add(&mut self, key: K, hash: u64) {
        let mask = self.table.len() - 1;
        let mut i = hash as usize & mask;
        let mut perturb = hash as usize;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for e in i..=i + probes {
                match &self.table[e] {
                    None => {
                        self.table[e] = Some((hash, key));
                        self.fill += 1;
                        if self.fill * 5 >= mask * 3 {
                            let used = self.fill;
                            self.resize(if used > 50000 { used * 2 } else { used * 4 });
                        }
                        return;
                    }
                    Some((h, k)) if *h == hash && *k == key => return,
                    Some(_) => {}
                }
            }
            perturb >>= 5;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb)) & mask;
        }
    }

    fn resize(&mut self, min_used: usize) {
        let mut size = 8;
        while size <= min_used {
            size <<= 1;
        }
        let old = std::mem::replace(&mut self.table, (0..size).map(|_| None).collect());
        let mask = size - 1;
        for (hash, key) in old.into_iter().flatten() {
            let mut i = hash as usize & mask;
            let mut perturb = hash as usize;
            'probe: loop {
                if self.table[i].is_none() {
                    self.table[i] = Some((hash, key));
                    break;
                }
                if i + LINEAR_PROBES <= mask {
                    for e in i + 1..=i + LINEAR_PROBES {
                        if self.table[e].is_none() {
                            self.table[e] = Some((hash, key));
                            break 'probe;
                        }
                    }
                }
                perturb >>= 5;
                i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb)) & mask;
            }
        }
    }

    /// Keys in iteration order.
    pub fn into_keys(self) -> Vec<K> {
        self.table.into_iter().flatten().map(|(_, k)| k).collect()
    }
}

/// CPython's `hash()` of a tuple of small non-negative ints (xxHash-based,
/// Objects/tupleobject.c, 64-bit).
pub fn py_tuple_hash(items: &[usize]) -> u64 {
    const P1: u64 = 11400714785074694791;
    const P2: u64 = 14029467366897019727;
    const P5: u64 = 2870177450012600261;
    let mut acc = P5;
    for &x in items {
        acc = acc.wrapping_add((x as u64).wrapping_mul(P2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(P1);
    }
    acc = acc.wrapping_add(items.len() as u64 ^ (P5 ^ 3527539));
    if acc == u64::MAX { 1546275796 } else { acc }
}

/// `enumerate_angles` in set iteration order.
fn enumerate_angles(nbrs: &[Vec<usize>]) -> Vec<[usize; 3]> {
    let mut set = PySet::new();
    for a in 0..nbrs.len() {
        for &b in &nbrs[a] {
            for &c in &nbrs[b] {
                if a == b || b == c || a == c {
                    continue;
                }
                let (x, z) = if a > c { (c, a) } else { (a, c) };
                let key = [x, b, z];
                set.add(key, py_tuple_hash(&key));
            }
        }
    }
    set.into_keys()
}

/// `enumerate_torsions` in set iteration order.
fn enumerate_torsions(nbrs: &[Vec<usize>]) -> Vec<[usize; 4]> {
    let mut set = PySet::new();
    for i0 in 0..nbrs.len() {
        for &i1 in &nbrs[i0] {
            for &i2 in &nbrs[i1] {
                if i2 == i0 {
                    continue;
                }
                for &i3 in &nbrs[i2] {
                    if i3 == i1 || i3 == i0 {
                        continue;
                    }
                    let key = if i0 < i3 { [i0, i1, i2, i3] } else { [i3, i2, i1, i0] };
                    set.add(key, py_tuple_hash(&key));
                }
            }
        }
    }
    set.into_keys()
}

/// Atomic number of an upper-case element symbol (`*` is 0).
fn atomic_number(element: &str) -> Option<usize> {
    if element == "*" {
        return Some(0);
    }
    ACCEPTED.iter().position(|(e, _)| *e == element).map(|k| k + 1)
}

/// Whether RDKit's sanitization accepts an atom with this atomic number,
/// formal charge and valence.
fn valence_accepted(z: usize, charge: i32, valence: usize) -> bool {
    if z == 0 {
        return true;
    }
    let masks = &ACCEPTED[z - 1].1;
    let k = charge - MIN_CHARGE;
    if k < 0 || k as usize >= masks.len() || valence > 15 {
        return false;
    }
    masks[k as usize] & (1 << valence) != 0
}

/// The largest valence RDKit allows a neutral atom (`getValenceList().back()`),
/// or `None` where RDKit skips the hypervalence test (any valence, or none).
fn max_valence(z: usize) -> Option<usize> {
    let mask = ACCEPTED.get(z.checked_sub(1)?)?.1[(-MIN_CHARGE) as usize];
    let max = 15 - mask.leading_zeros() as usize;
    if mask == 0xffff || mask == 0 || max == 0 { None } else { Some(max) }
}

/// `QueryOps::isMetal`: everything except the Marvin non-metals.
fn is_metal(z: usize) -> bool {
    ![0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86].contains(&z)
}

/// A bond of the RDKit molecule: atoms and order (0 for dative, which counts
/// for the acceptor only).
#[derive(Clone, Copy)]
struct RdBond {
    a: usize,
    b: usize,
    order: usize,
}

/// `Chem.SanitizeMol` as far as it can fail for these molecules: the clean-up
/// steps (nitro and azide nitrogens, phosphorus and halogen oxides, bonds to
/// metals made dative), then RDKit's valence check.
fn rdkit_sanitizes(z: &[usize], charge: &mut [i32], bonds: &mut [RdBond]) -> bool {
    let n = z.len();
    // bonds of each atom in insertion order (RDKit's neighbour order)
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (k, b) in bonds.iter().enumerate() {
        incident[b.a].push(k);
        incident[b.b].push(k);
    }
    let other = |b: &RdBond, a: usize| if b.a == a { b.b } else { b.a };
    let valence = |bonds: &[RdBond], a: usize| -> usize {
        incident[a].iter().map(|&k| if bonds[k].order == 0 { usize::from(bonds[k].b == a) } else { bonds[k].order }).sum()
    };
    // nitrogensCleanup
    let mut considered = Vec::new();
    for a in 0..n {
        if z[a] != 7 || charge[a] != 0 || valence(bonds, a) != 5 {
            continue;
        }
        considered.push(a);
        for &k in &incident[a] {
            let o = other(&bonds[k], a);
            if z[o] == 8 && charge[o] == 0 && bonds[k].order == 2 {
                bonds[k].order = 1;
                charge[a] = 1;
                charge[o] = -1;
                break;
            }
        }
    }
    for &a in &considered {
        for &k in &incident[a] {
            let o = other(&bonds[k], a);
            if z[o] == 7 && charge[o] == 0 && bonds[k].order == 3 {
                bonds[k].order = 2;
                charge[a] = 1;
                charge[o] = -1;
                break;
            }
        }
    }
    for a in 0..n {
        match z[a] {
            15 if charge[a] == 0 && valence(bonds, a) == 5 && incident[a].len() == 3 => {
                let mut double_to_o = None;
                let mut double_to_c_or_n = false;
                for &k in &incident[a] {
                    let o = other(&bonds[k], a);
                    if z[o] == 8 && charge[o] == 0 && bonds[k].order == 2 {
                        double_to_o = Some((k, o));
                    } else if (z[o] == 6 || z[o] == 7) && incident[o].len() >= 2 && bonds[k].order == 2 {
                        double_to_c_or_n = true;
                    }
                }
                if let (true, Some((k, o))) = (double_to_c_or_n, double_to_o) {
                    charge[o] = -1;
                    bonds[k].order = 1;
                    charge[a] = 1;
                }
            }
            17 | 35 | 53 if charge[a] == 0 && matches!(valence(bonds, a), 3 | 5 | 7) => {
                if incident[a].iter().all(|&k| z[other(&bonds[k], a)] == 8) {
                    let mut c = 0;
                    for &k in &incident[a] {
                        if bonds[k].order == 2 {
                            bonds[k].order = 1;
                            c += 1;
                            charge[other(&bonds[k], a)] = -1;
                        }
                    }
                    charge[a] = c;
                }
            }
            _ => {}
        }
    }
    // cleanUpOrganometallics: a hypervalent non-metal donates one single bond
    // to a metal as a dative bond
    for a in 0..n {
        if is_metal(z[a]) || [1, 2, 9, 10].contains(&z[a]) {
            continue;
        }
        let eff = z[a] as i64 - charge[a] as i64;
        if eff <= 0 {
            continue;
        }
        let Some(max) = max_valence(eff as usize) else { continue };
        if valence(bonds, a) <= max {
            continue;
        }
        if let Some(&k) = incident[a].iter().find(|&&k| bonds[k].order == 1 && is_metal(z[other(&bonds[k], a)])) {
            let m = other(&bonds[k], a);
            bonds[k] = RdBond { a, b: m, order: 0 };
        }
    }
    (0..n).all(|a| valence_accepted(z[a], charge[a], valence(bonds, a)))
}

/// Format like Python `'%0.<n>f' % v` and read the result back, as the
/// dictionary goes through CIF text.
fn rounded(v: f64, decimals: usize) -> f64 {
    format!("{:.*}", decimals, v).parse().unwrap_or(v)
}

/// Build the dictionary from a CCD file's text, or `None` where Reduce2's
/// construction fails.
pub fn comp_from_ccd(text: &str, source: &Path) -> Option<Comp> {
    let doc = cif::parse(text);
    let block = doc.blocks.first()?;
    let cc = block.category("_chem_comp")?;
    let atoms = block.category("_chem_comp_atom")?;
    let col = |tag: &str| atoms.col(tag);
    let (c_id, c_type, c_charge) = (col("atom_id")?, col("type_symbol")?, col("charge"));
    let n = atoms.nrows();
    let get = |r: usize, c: Option<usize>| c.map_or("?", |c| atoms.get(r, c));

    // coordinates: the ideal set without its all-unknown rows, else the model set
    let coords = |names: [&str; 3]| -> Vec<[String; 3]> {
        let cols = names.map(|t| col(t));
        (0..n)
            .map(|r| cols.map(|c| get(r, c).to_string()))
            .filter(|xyz| xyz.iter().any(|v| v != "?"))
            .collect()
    };
    let mut xyzs = coords(["pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal"]);
    if xyzs.is_empty() {
        xyzs = coords(["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]);
    }
    if xyzs.is_empty() {
        return None;
    }

    // the RDKit molecule: atoms, positions, bonds
    let mut ids: Vec<String> = Vec::with_capacity(n);
    let mut elements: Vec<String> = Vec::with_capacity(n);
    let mut numbers: Vec<usize> = Vec::with_capacity(n);
    let mut charges: Vec<i32> = Vec::with_capacity(n);
    let mut pos: Vec<P3> = vec![[0.0; 3]; n];
    for r in 0..n {
        let element = get(r, Some(c_type)).to_ascii_uppercase();
        numbers.push(atomic_number(&element)?);
        charges.push(get(r, c_charge).trim().parse::<i32>().ok()?);
        let xyz = xyzs.get(r)?;
        if xyz[0] != "?" {
            let mut p = [0.0; 3];
            for k in 0..3 {
                p[k] = xyz[k].trim().parse::<f64>().ok()?;
            }
            pos[r] = p;
        }
        ids.push(get(r, Some(c_id)).to_string());
        elements.push(element);
    }
    let index_of = |name: &str| -> Option<usize> { ids.iter().rposition(|x| x == name) };
    let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut rd_bonds: Vec<RdBond> = Vec::new();
    let mut ccd_bonds: Vec<(String, String, String)> = Vec::new();
    if let Some(bonds) = block.category("_chem_comp_bond") {
        let (b1, b2, bo) = (bonds.col("atom_id_1"), bonds.col("atom_id_2"), bonds.col("value_order"));
        for r in 0..bonds.nrows() {
            let field = |c: Option<usize>| c.map_or("?", |c| bonds.get(r, c));
            let (a1, a2, order) = (field(b1), field(b2), field(bo));
            let i = index_of(a1)?;
            let j = index_of(a2)?;
            let order_n = match order {
                "SING" => 1,
                "DOUB" => 2,
                "TRIP" => 3,
                _ => return None,
            };
            if i == j || nbrs[i].contains(&j) {
                return None;
            }
            nbrs[i].push(j);
            nbrs[j].push(i);
            rd_bonds.push(RdBond { a: i, b: j, order: order_n });
            ccd_bonds.push((a1.to_string(), a2.to_string(), order.to_string()));
        }
    }
    if !rdkit_sanitizes(&numbers, &mut charges.clone(), &mut rd_bonds) {
        return None;
    }

    // the dictionary
    let mut comp = Comp {
        id: cc.get_tag(0, "id").unwrap_or("").to_string(),
        group: cc.get_tag(0, "type").unwrap_or("").to_string(),
        source: source.to_path_buf(),
        from_ccd: true,
        ..Default::default()
    };
    for k in 0..n {
        comp.atoms.push(CompAtom { id: ids[k].clone(), type_symbol: elements[k].clone(), type_energy: None });
    }
    for (a1, a2, order) in ccd_bonds {
        let (i, j) = (index_of(&a1)?, index_of(&a2)?);
        comp.bonds.push(CompBond {
            a1,
            a2,
            type_: order,
            value_dist: Some(rounded(bond_length(&pos, i, j) * 0.9, 3)),
            esd: Some(0.1),
            value_dist_neutron: None,
        });
    }
    for [a, b, c] in enumerate_angles(&nbrs) {
        comp.angles.push(CompAngle {
            a1: ids[a].clone(),
            a2: ids[b].clone(),
            a3: ids[c].clone(),
            value: Some(rounded(angle_deg(&pos, a, b, c)?, 1)),
            esd: Some(1.0),
        });
    }
    for (k, [a, b, c, d]) in enumerate_torsions(&nbrs).into_iter().enumerate() {
        comp.tors.push(CompTor {
            id: format!("Var_{:03}", k),
            a: [ids[a].clone(), ids[b].clone(), ids[c].clone(), ids[d].clone()],
            value: Some(rounded(dihedral_deg(&pos, a, b, c, d)?, 1)),
            esd: Some(1.0),
            period: 1,
            alt_values: None,
        });
    }
    Some(comp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ccd(atoms: &str, bonds: &str) -> String {
        format!(
            "data_TST\n_chem_comp.id TST\n_chem_comp.type NON-POLYMER\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
             _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_model_Cartn_x_ideal\n\
             _chem_comp_atom.pdbx_model_Cartn_y_ideal\n_chem_comp_atom.pdbx_model_Cartn_z_ideal\n{}\nloop_\n\
             _chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n{}\n",
            atoms, bonds
        )
    }

    #[test]
    fn builds_restraints_like_reduce2() {
        // water-like O with two H: bonds at 0.9 times the ideal length
        let text = ccd(
            "TST O O 0 0.000 0.000 0.000\nTST H1 H 0 1.000 0.000 0.000\nTST H2 H 0 0.000 1.000 0.000",
            "TST O H1 SING\nTST O H2 SING",
        );
        let c = comp_from_ccd(&text, Path::new("x")).unwrap();
        assert!(c.from_ccd && c.group == "NON-POLYMER" && c.atoms.iter().all(|a| a.type_energy.is_none()));
        assert_eq!(c.bonds.iter().map(|b| b.value_dist.unwrap()).collect::<Vec<_>>(), [0.9, 0.9]);
        assert_eq!(c.angles.len(), 1);
        assert_eq!(c.angles[0].value, Some(90.0));
        assert!(c.tors.is_empty());
    }

    #[test]
    fn rejects_what_rdkit_rejects() {
        let atoms = "TST C1 C 0 0 0 0\nTST C2 C 0 1.5 0 0";
        assert!(comp_from_ccd(&ccd(atoms, "TST C1 C2 AROM"), Path::new("x")).is_none(), "bond order");
        assert!(comp_from_ccd(&ccd("TST C1 C ? 0 0 0\nTST C2 C 0 1.5 0 0", "TST C1 C2 SING"), Path::new("x")).is_none(), "charge");
        assert!(comp_from_ccd(&ccd("TST C1 X 0 0 0 0\nTST C2 C 0 1.5 0 0", "TST C1 C2 SING"), Path::new("x")).is_none(), "element");
        // three bonds on a neutral O: too many, unless one goes to a metal
        let o3 = "TST O O 0 0 0 0\nTST C1 C 0 1 0 0\nTST C2 C 0 0 1 0\nTST C3 C 0 0 0 1";
        assert!(comp_from_ccd(&ccd(o3, "TST O C1 SING\nTST O C2 SING\nTST O C3 SING"), Path::new("x")).is_none());
        let o2m = "TST O O 0 0 0 0\nTST C1 C 0 1 0 0\nTST C2 C 0 0 1 0\nTST FE FE 0 0 0 2";
        assert!(comp_from_ccd(&ccd(o2m, "TST O C1 SING\nTST O C2 SING\nTST O FE SING"), Path::new("x")).is_some());
    }

    #[test]
    fn tuple_hashes_match_cpython() {
        assert_eq!(py_tuple_hash(&[]), 5740354900026072187);
        assert_eq!(py_tuple_hash(&[1, 2]) as i64, -3550055125485641917);
        assert_eq!(py_tuple_hash(&[0, 1, 2, 3]) as i64, -8281178343658874797);
    }

    /// 300 random 3-tuples added to a CPython 3.12 set, and the set's order.
    #[test]
    fn set_order_matches_cpython() {
        let text = include_str!("../tests/data/pyset_order.txt");
        let mut lines = text.lines();
        let parse = |l: &str| -> Vec<[usize; 3]> {
            l.split(';')
                .map(|t| {
                    let v: Vec<usize> = t.split(',').map(|x| x.parse().unwrap()).collect();
                    [v[0], v[1], v[2]]
                })
                .collect()
        };
        let keys = parse(lines.next().unwrap());
        let expected = parse(lines.next().unwrap());
        let mut set = PySet::new();
        for k in keys.iter().chain(keys.iter().take(20)) {
            set.add(*k, py_tuple_hash(k));
        }
        assert_eq!(set.into_keys(), expected);
    }
}
