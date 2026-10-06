//! Flat, index-based view of a structure used by the optimizer and Movers.

use crate::geom::Vec3;
use crate::model::Structure;
use crate::movers::RidingRef;
use crate::probe::AtomInfo;
use crate::resclass::{self, ResClass};
use std::sync::Arc;

/// Static labels of an atom (hierarchy context), mirroring what Reduce2
/// reads through `atom.parent()...`.
#[derive(Clone, Debug, Default)]
pub struct AtomLabels {
    /// Stripped, upper-case atom name.
    pub name: String,
    /// Raw (padded) atom name.
    pub raw_name: String,
    /// Stripped, upper-case element.
    pub element: String,
    // Residue-level labels are shared by the atoms of an atom group.
    pub resname: Arc<str>,
    /// Raw (padded) residue name as stored by iotbx.
    pub resname_raw: Arc<str>,
    pub chain: Arc<str>,
    pub resseq: i32,
    pub resseq_raw: Arc<str>,
    pub icode: Arc<str>,
    /// Raw atom-group altloc ("" for the main conformation).
    pub altloc: Arc<str>,
    pub model_index: u32,
    pub model_id: Arc<str>,
    /// Unique id of the atom group (parent) and residue group.
    pub ag: u32,
    pub rg: u32,
    pub hetero: bool,
    pub is_water: bool,
    pub class: Option<ResClass>,
}

pub struct World {
    pub labels: Vec<AtomLabels>,
    pub pos: Vec<Vec3>,
    pub occ: Vec<f64>,
    pub b: Vec<f64>,
    pub info: Vec<AtomInfo>,
    pub bonded: Vec<Vec<u32>>,
    pub riding: Vec<Option<RidingRef>>,
    /// Number of real atoms; phantom hydrogens are appended after these.
    pub n_real: usize,
}

impl World {
    #[inline]
    pub fn len(&self) -> usize {
        self.pos.len()
    }
    #[inline]
    pub fn elem(&self, a: u32) -> &str {
        &self.labels[a as usize].element
    }
    #[inline]
    pub fn name(&self, a: u32) -> &str {
        &self.labels[a as usize].name
    }
    #[inline]
    pub fn resname(&self, a: u32) -> &str {
        &self.labels[a as usize].resname
    }
    #[inline]
    pub fn is_h(&self, a: u32) -> bool {
        let e = self.elem(a);
        e == "H" || e == "D"
    }
    #[inline]
    pub fn altloc(&self, a: u32) -> &str {
        &self.labels[a as usize].altloc
    }
    #[inline]
    pub fn set_pos(&mut self, a: u32, p: Vec3) {
        self.pos[a as usize] = p;
    }
    pub fn is_positive_ion(&self, a: u32) -> bool {
        resclass::element_is_positive_ion(self.elem(a))
    }
    /// `compatibleConformations` from probe Helpers (hierarchy altlocs).
    pub fn compatible(&self, a: u32, b: u32) -> bool {
        let x = self.altloc(a);
        let y = self.altloc(b);
        x.is_empty() || x == " " || y.is_empty() || y == " " || x == y
    }

    /// Sorting key used by MoverSingleHydrogenRotator for acceptors:
    /// chain + resname + resseq + name + altloc, all raw.
    pub fn atom_id_string(&self, a: u32) -> String {
        let l = &self.labels[a as usize];
        format!("{}{}{}{}{}", l.chain, l.resname_raw, l.resseq, l.raw_name, l.altloc)
    }

    /// `_ResNameAndID(a)`: "chain A THR 1" style residue description.
    pub fn res_name_and_id(&self, a: u32) -> String {
        let l = &self.labels[a as usize];
        format!("chain {} {}{} {}{}", l.chain, l.altloc, l.resname, l.resseq, l.icode)
    }

    /// Build a world from a structure. Bonds, atom info and riding references
    /// are filled in by the caller.
    pub fn from_structure(st: &Structure) -> World {
        let n = st.atoms_size();
        let mut labels = Vec::with_capacity(n);
        let mut pos = Vec::with_capacity(n);
        let mut occ = Vec::with_capacity(n);
        let mut bfac = Vec::with_capacity(n);
        let mut ag_id = 0u32;
        let mut rg_id = 0u32;
        for (mi, m) in st.models.iter().enumerate() {
            let model_id: Arc<str> = Arc::from(m.id.as_str());
            for c in &m.chains {
                let chain: Arc<str> = Arc::from(c.id.as_str());
                for rg in &c.residue_groups {
                    let resseq_raw: Arc<str> = Arc::from(rg.resseq.as_str());
                    let icode: Arc<str> = Arc::from(rg.icode.trim());
                    for ag in &rg.atom_groups {
                        let resname: Arc<str> = Arc::from(ag.resname.trim().to_ascii_uppercase());
                        let resname_raw: Arc<str> = Arc::from(ag.resname.as_str());
                        let altloc: Arc<str> = Arc::from(ag.altloc.as_str());
                        let class = resclass::get_class(&ag.resname);
                        for a in &ag.atoms {
                            labels.push(AtomLabels {
                                name: a.name.trim().to_ascii_uppercase(),
                                raw_name: a.name.clone(),
                                element: a.elem().to_ascii_uppercase(),
                                resname: resname.clone(),
                                resname_raw: resname_raw.clone(),
                                chain: chain.clone(),
                                resseq: rg.resseq_as_int(),
                                resseq_raw: resseq_raw.clone(),
                                icode: icode.clone(),
                                altloc: altloc.clone(),
                                model_index: mi as u32,
                                model_id: model_id.clone(),
                                ag: ag_id,
                                rg: rg_id,
                                hetero: a.hetero,
                                is_water: class == ResClass::CommonWater,
                                class: Some(class),
                            });
                            pos.push(a.xyz);
                            occ.push(a.occ);
                            bfac.push(a.b);
                        }
                        ag_id += 1;
                    }
                    rg_id += 1;
                }
            }
        }
        World {
            labels,
            pos,
            occ,
            b: bfac,
            info: vec![AtomInfo::default(); n],
            bonded: vec![Vec::new(); n],
            riding: vec![None; n],
            n_real: n,
        }
    }

    /// Append a phantom hydrogen on a water oxygen; returns its index.
    pub fn add_phantom(&mut self, parent: u32, p: Vec3, info: AtomInfo) -> u32 {
        let mut l = self.labels[parent as usize].clone();
        l.name = "H?".into();
        l.raw_name = " H?".into();
        l.element = "H".into();
        let idx = self.pos.len() as u32;
        self.labels.push(l);
        self.pos.push(p);
        self.occ.push(self.occ[parent as usize]);
        self.b.push(self.b[parent as usize]);
        self.info.push(info);
        self.bonded.push(vec![parent]);
        self.riding.push(None);
        idx
    }

    /// Drop all phantom hydrogens.
    pub fn clear_phantoms(&mut self) {
        let n = self.n_real;
        self.labels.truncate(n);
        self.pos.truncate(n);
        self.occ.truncate(n);
        self.b.truncate(n);
        self.info.truncate(n);
        self.bonded.truncate(n);
        self.riding.truncate(n);
    }
}
