//! CCP4/GeoStd monomer library access (port of the parts of
//! mmtbx.monomer_library.server that Reduce2 uses), plus the wwPDB chemical
//! component dictionary (CCD).

use crate::cif;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug, Default)]
pub struct CompAtom {
    pub id: String,
    pub type_symbol: String,
    /// `None` for dictionaries without energy types.
    pub type_energy: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct CompBond {
    pub a1: String,
    pub a2: String,
    pub type_: String,
    pub value_dist: Option<f64>,
    pub esd: Option<f64>,
    pub value_dist_neutron: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct CompAngle {
    pub a1: String,
    pub a2: String,
    pub a3: String,
    pub value: Option<f64>,
    pub esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct CompTor {
    pub id: String,
    pub a: [String; 4],
    pub value: Option<f64>,
    pub esd: Option<f64>,
    pub period: i32,
    pub alt_values: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct CompChir {
    pub id: String,
    pub centre: String,
    pub a: [String; 3],
    pub volume_sign: String,
}

#[derive(Clone, Debug, Default)]
pub struct CompPlaneAtom {
    pub plane_id: String,
    pub atom: String,
    pub esd: Option<f64>,
}

/// A residue dictionary (`comp_comp_id`), possibly with modifications applied.
#[derive(Clone, Debug, Default)]
pub struct Comp {
    pub id: String,
    pub group: String,
    pub atoms: Vec<CompAtom>,
    pub bonds: Vec<CompBond>,
    pub angles: Vec<CompAngle>,
    pub tors: Vec<CompTor>,
    pub chirs: Vec<CompChir>,
    pub planes: Vec<CompPlaneAtom>,
    pub source: PathBuf,
    pub applied_mods: Vec<String>,
    pub is_terminus: bool,
    /// Built from the CCD for a residue the libraries lack (Reduce2's
    /// `auto_<resname>` dictionary); its atoms have no energy types.
    pub from_ccd: bool,
}

impl Comp {
    pub fn atom(&self, id: &str) -> Option<&CompAtom> {
        self.atoms.iter().find(|a| a.id == id)
    }
    pub fn has_atom(&self, id: &str) -> bool {
        self.atoms.iter().any(|a| a.id == id)
    }
    pub fn is_h(&self, id: &str) -> bool {
        self.atom(id).map(|a| a.type_symbol == "H" || a.type_symbol == "D").unwrap_or(false)
    }

    /// `comp_comp_id.test_for_peptide`.
    pub fn test_for_peptide(&self) -> bool {
        for n in ["N", "CA", "C", "O"] {
            if !self.has_atom(n) {
                return false;
            }
        }
        let has = |x: &str, y: &str, check: &dyn Fn(&CompBond) -> bool| {
            self.bonds.iter().any(|b| ((b.a1 == x && b.a2 == y) || (b.a1 == y && b.a2 == x)) && check(b))
        };
        has("CA", "N", &|_| true) && has("C", "CA", &|_| true) && has("C", "O", &|b| b.type_ != "single")
    }

    /// RNA/DNA classification of the dictionary (`test_for_rna_dna`), returning
    /// Some(true) for RNA, Some(false) for DNA.
    pub fn test_for_rna_dna(&self) -> Option<bool> {
        let pick = |alts: &[&str]| -> Option<String> {
            alts.iter().find(|a| self.has_atom(a)).map(|s| s.to_string())
        };
        let p = pick(&["P"])?;
        let op1 = pick(&["OP1", "O1P"])?;
        let op2 = pick(&["OP2", "O2P"])?;
        let star = self.has_atom("O5*");
        let n = |base: &str| -> String { if star { base.replace('\'', "*") } else { base.to_string() } };
        let o5 = pick(&["O5'", "O5*"])?;
        for a in ["C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"] {
            if !self.has_atom(&n(a)) && !self.has_atom(a) {
                return None;
            }
        }
        let bonded = |x: &str, y: &str| {
            let xs = [x.to_string(), n(x)];
            let ys = [y.to_string(), n(y)];
            self.bonds.iter().any(|b| {
                (xs.contains(&b.a1) && ys.contains(&b.a2)) || (ys.contains(&b.a1) && xs.contains(&b.a2))
            })
        };
        let req = [
            (op1.as_str(), p.as_str()),
            (op2.as_str(), p.as_str()),
            (o5.as_str(), p.as_str()),
            ("C1'", "C2'"),
            ("C2'", "C3'"),
            ("C3'", "C4'"),
            ("C3'", "O3'"),
            ("C4'", "C5'"),
            ("C4'", "O4'"),
            ("C1'", "O4'"),
            ("C5'", "O5'"),
        ];
        for (x, y) in req {
            if !bonded(x, y) {
                return None;
            }
        }
        Some(bonded("C2'", "O2'"))
    }

    /// `hydrogen_deuterium_aliases`: "D"+id[1:] -> id for H atoms starting with H.
    pub fn hd_alias(&self, name: &str) -> Option<String> {
        if !name.starts_with('D') {
            return None;
        }
        let cand = format!("H{}", &name[1..]);
        if self.atoms.iter().any(|a| a.id == cand && (a.type_symbol == "H")) { Some(cand) } else { None }
    }
}

// ----------------------------------------------------------------------------
// Links and modifications

#[derive(Clone, Debug, Default)]
pub struct LinkBond {
    pub c1: u8,
    pub a1: String,
    pub c2: u8,
    pub a2: String,
    pub value_dist: Option<f64>,
    pub esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct LinkAngle {
    pub c: [u8; 3],
    pub a: [String; 3],
    pub value: Option<f64>,
    pub esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct LinkTor {
    pub id: String,
    pub c: [u8; 4],
    pub a: [String; 4],
    pub value: Option<f64>,
    pub esd: Option<f64>,
    pub period: i32,
}

#[derive(Clone, Debug, Default)]
pub struct LinkChir {
    pub c: [u8; 4],
    pub a: [String; 4],
    pub volume_sign: String,
}

#[derive(Clone, Debug, Default)]
pub struct LinkPlaneAtom {
    pub plane_id: String,
    pub c: u8,
    pub atom: String,
    pub esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct ChemLink {
    pub id: String,
    pub comp_id_1: String,
    pub mod_id_1: String,
    pub group_comp_1: String,
    pub comp_id_2: String,
    pub mod_id_2: String,
    pub group_comp_2: String,
    pub name: String,
    /// Whether a `_chem_link` table row describes the link; cctbx leaves the
    /// fields of the others None rather than empty.
    pub listed: bool,
    pub bonds: Vec<LinkBond>,
    pub angles: Vec<LinkAngle>,
    pub tors: Vec<LinkTor>,
    pub chirs: Vec<LinkChir>,
    pub planes: Vec<LinkPlaneAtom>,
}

#[derive(Clone, Debug, Default)]
pub struct ModAtom {
    pub function: String,
    pub atom_id: String,
    pub new_atom_id: String,
    pub new_type_symbol: String,
    pub new_type_energy: String,
}

#[derive(Clone, Debug, Default)]
pub struct ModBond {
    pub function: String,
    pub a1: String,
    pub a2: String,
    pub new_type: String,
    pub new_value_dist: Option<f64>,
    pub new_esd: Option<f64>,
    pub new_value_dist_neutron: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct ModAngle {
    pub function: String,
    pub a: [String; 3],
    pub new_value: Option<f64>,
    pub new_esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct ModTor {
    pub function: String,
    pub id: String,
    pub a: [String; 4],
    pub new_value: Option<f64>,
    pub new_esd: Option<f64>,
    pub new_period: Option<i32>,
    pub new_alt: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct ModChir {
    pub function: String,
    pub id: String,
    pub centre: String,
    pub a: [String; 3],
    pub new_volume_sign: String,
}

#[derive(Clone, Debug, Default)]
pub struct ModPlane {
    pub function: String,
    pub plane_id: String,
    pub atom: String,
    pub new_esd: Option<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct ChemMod {
    pub id: String,
    pub name: String,
    pub atoms: Vec<ModAtom>,
    pub bonds: Vec<ModBond>,
    pub angles: Vec<ModAngle>,
    pub tors: Vec<ModTor>,
    pub chirs: Vec<ModChir>,
    pub planes: Vec<ModPlane>,
}

#[derive(Clone, Debug, Default)]
pub struct EnerAtom {
    pub type_: String,
    pub hb_type: String,
    pub vdw_radius: Option<f64>,
    pub vdwh_radius: Option<f64>,
    pub ion_radius: Option<f64>,
    pub element: String,
    pub vdw_radius_neutron: Option<f64>,
}

/// CCD entry: atoms with model/ideal coordinates and bonds with orders.
#[derive(Clone, Debug, Default)]
pub struct CcdEntry {
    pub id: String,
    pub type_: String,
    /// (atom_id, element, model xyz, ideal xyz)
    pub atoms: Vec<(String, String, Option<[f64; 3]>, Option<[f64; 3]>)>,
    /// (atom1, atom2, order)
    pub bonds: Vec<(String, String, String)>,
}

pub struct MonLib {
    pub root: PathBuf,
    pub comp_synonyms: FxHashMap<String, String>,
    pub atom_synonyms: FxHashMap<String, FxHashMap<String, String>>,
    pub links: Vec<ChemLink>,
    pub link_index: FxHashMap<String, usize>,
    /// Links in cctbx's `link_link_id_list` order: the link blocks of the list
    /// files that hold restraints, then the RNA/DNA chain links.
    pub link_list_order: Vec<usize>,
    pub mods: FxHashMap<String, ChemMod>,
    pub ener: FxHashMap<String, EnerAtom>,
    comp_cache: RwLock<FxHashMap<String, Option<Arc<Comp>>>>,
    ccd_cache: RwLock<FxHashMap<String, Option<Arc<CcdEntry>>>>,
    ccd_comp_cache: RwLock<FxHashMap<String, Option<Arc<Comp>>>>,
    variant_cache: RwLock<FxHashMap<String, Option<Arc<Comp>>>>,
    sites_cache: RwLock<FxHashMap<(PathBuf, String), Option<Arc<FxHashMap<String, crate::geom::Vec3>>>>>,
    /// User-supplied dictionaries (restraint CIF files), which win over the library.
    user_comps: FxHashMap<String, Arc<Comp>>,
}

fn s(v: &str) -> String {
    if cif::is_null(v) { String::new() } else { v.to_string() }
}
fn opt_f(v: Option<&str>) -> Option<f64> {
    v.and_then(cif::parse_f64)
}

fn get<'a>(c: &cif::Category<'a>, row: usize, col: Option<usize>) -> &'a str {
    match col {
        Some(k) => c.get(row, k),
        None => ".",
    }
}

fn comp_digit(v: &str) -> u8 {
    v.trim().parse::<u8>().unwrap_or(1)
}

impl MonLib {
    /// Locate `chem_data` from (in order) an explicit path, `REDUCE3_CHEM_DATA`,
    /// `MMTBX_CCP4_MONOMER_LIB`'s parent, `CLIBD_MON`'s parent, or a few
    /// conventional install locations.
    pub fn locate(explicit: Option<&str>) -> Option<PathBuf> {
        let mut cands: Vec<PathBuf> = Vec::new();
        if let Some(e) = explicit {
            cands.push(PathBuf::from(e));
        }
        for var in ["REDUCE3_CHEM_DATA", "CHEM_DATA"] {
            if let Ok(v) = std::env::var(var) {
                cands.push(PathBuf::from(v));
            }
        }
        for var in ["MMTBX_CCP4_MONOMER_LIB", "CLIBD_MON"] {
            if let Ok(v) = std::env::var(var) {
                if let Some(p) = Path::new(&v).parent() {
                    cands.push(p.to_path_buf());
                }
            }
        }
        if let Ok(v) = std::env::var("CONDA_PREFIX") {
            if let Ok(rd) = std::fs::read_dir(Path::new(&v).join("lib")) {
                for e in rd.flatten() {
                    cands.push(e.path().join("site-packages/chem_data"));
                }
            }
        }
        cands.into_iter().find(|c| c.join("geostd").is_dir() && c.join("mon_lib").is_dir())
    }

    pub fn load(root: &Path) -> Result<MonLib, String> {
        let mut ml = MonLib {
            root: root.to_path_buf(),
            comp_synonyms: FxHashMap::default(),
            atom_synonyms: FxHashMap::default(),
            links: Vec::new(),
            link_index: FxHashMap::default(),
            link_list_order: Vec::new(),
            mods: FxHashMap::default(),
            ener: FxHashMap::default(),
            comp_cache: RwLock::new(FxHashMap::default()),
            ccd_cache: RwLock::new(FxHashMap::default()),
            ccd_comp_cache: RwLock::new(FxHashMap::default()),
            variant_cache: RwLock::new(FxHashMap::default()),
            sites_cache: RwLock::new(FxHashMap::default()),
            user_comps: FxHashMap::default(),
        };
        let read = |p: PathBuf| -> Result<String, String> {
            std::fs::read_to_string(&p).map_err(|e| format!("cannot read {}: {}", p.display(), e))
        };
        // mon_lib_list.cif with the geostd list merged on top; cctbx finds the
        // copy under geostd first (`mon_lib_list_cif`)
        let geostd_mlist = root.join("geostd/list/mon_lib_list.cif");
        let mlist = if geostd_mlist.exists() { read(geostd_mlist)? } else { read(root.join("mon_lib/list/mon_lib_list.cif"))? };
        let glist = read(root.join("geostd/list/geostd_list.cif")).unwrap_or_default();
        let mdoc = cif::parse(&mlist);
        let gdoc = cif::parse(&glist);
        // list tables: mon_lib rows then geostd rows (later rows win)
        let mut link_rows: Vec<ChemLink> = Vec::new();
        let mut mod_names: FxHashMap<String, String> = FxHashMap::default();
        for doc in [&mdoc, &gdoc] {
            for b in &doc.blocks {
                if let Some(c) = b.category("_chem_comp_synonym") {
                    let (ci, ai) = (c.col("comp_id"), c.col("comp_alternative_id"));
                    for r in 0..c.nrows() {
                        let id = s(get(c, r, ci));
                        let alt = s(get(c, r, ai));
                        if !id.is_empty() && !alt.is_empty() {
                            ml.comp_synonyms.insert(alt.to_ascii_uppercase(), id);
                        }
                    }
                }
                if let Some(c) = b.category("_chem_comp_synonym_atom") {
                    let (ci, ai, aa) = (c.col("comp_id"), c.col("atom_id"), c.col("atom_alternative_id"));
                    for r in 0..c.nrows() {
                        let comp = s(get(c, r, ci));
                        ml.atom_synonyms.entry(comp).or_default().insert(s(get(c, r, aa)), s(get(c, r, ai)));
                    }
                }
                if let Some(c) = b.category("_chem_link") {
                    let cols: Vec<Option<usize>> = ["id", "comp_id_1", "mod_id_1", "group_comp_1", "comp_id_2", "mod_id_2", "group_comp_2", "name"]
                        .iter()
                        .map(|t| c.col(t))
                        .collect();
                    for r in 0..c.nrows() {
                        let g = |k: usize| s(get(c, r, cols[k]));
                        link_rows.push(ChemLink {
                            id: g(0),
                            comp_id_1: g(1),
                            mod_id_1: g(2),
                            group_comp_1: g(3),
                            comp_id_2: g(4),
                            mod_id_2: g(5),
                            group_comp_2: g(6),
                            name: g(7),
                            listed: true,
                            ..Default::default()
                        });
                    }
                }
                if let Some(c) = b.category("_chem_mod") {
                    let (ii, ni) = (c.col("id"), c.col("name"));
                    for r in 0..c.nrows() {
                        mod_names.insert(s(get(c, r, ii)), s(get(c, r, ni)));
                    }
                }
            }
        }
        // later rows win for the same id, but keep the first position in the list
        for l in link_rows {
            match ml.link_index.get(&l.id) {
                Some(&i) => ml.links[i] = l,
                None => {
                    ml.link_index.insert(l.id.clone(), ml.links.len());
                    ml.links.push(l);
                }
            }
        }
        // link and mod definition blocks: geostd blocks replace mon_lib ones
        let mut link_blocks: FxHashMap<String, &cif::Block> = FxHashMap::default();
        let mut mod_blocks: FxHashMap<String, &cif::Block> = FxHashMap::default();
        for doc in [&mdoc, &gdoc] {
            for b in &doc.blocks {
                if let Some(id) = b.name.strip_prefix("link_") {
                    if !id.ends_with("_list") {
                        link_blocks.insert(id.to_string(), b);
                    }
                } else if let Some(id) = b.name.strip_prefix("mod_") {
                    if !id.ends_with("_list") {
                        mod_blocks.insert(id.to_string(), b);
                    }
                }
            }
        }
        // RNA/DNA chain links and pucker mods from geostd/rna_dna
        let mut extra_docs: Vec<String> = Vec::new();
        for f in [
            "chain_link_rna2p", "chain_link_rna3p", "mod_rna2p", "mod_rna3p", "mod_rna2p_pur", "mod_rna3p_pur",
            "mod_rna2p_pyr", "mod_rna3p_pyr",
        ] {
            if let Ok(t) = read(root.join(format!("geostd/rna_dna/{}.cif", f))) {
                extra_docs.push(t);
            }
        }
        let extra_parsed: Vec<cif::Document> = extra_docs.iter().map(|t| cif::parse(t)).collect();
        for doc in &extra_parsed {
            for b in &doc.blocks {
                if let Some(c) = b.category("_chem_link") {
                    let cols: Vec<Option<usize>> = ["id", "comp_id_1", "mod_id_1", "group_comp_1", "comp_id_2", "mod_id_2", "group_comp_2", "name"]
                        .iter()
                        .map(|t| c.col(t))
                        .collect();
                    for r in 0..c.nrows() {
                        let g = |k: usize| s(get(c, r, cols[k]));
                        let l = ChemLink {
                            id: g(0),
                            comp_id_1: g(1),
                            mod_id_1: g(2),
                            group_comp_1: g(3),
                            comp_id_2: g(4),
                            mod_id_2: g(5),
                            group_comp_2: g(6),
                            name: g(7),
                            listed: true,
                            ..Default::default()
                        };
                        if !ml.link_index.contains_key(&l.id) {
                            ml.link_index.insert(l.id.clone(), ml.links.len());
                            ml.links.push(l);
                        }
                    }
                }
                if let Some(c) = b.category("_chem_mod") {
                    let (ii, ni) = (c.col("id"), c.col("name"));
                    for r in 0..c.nrows() {
                        mod_names.insert(s(get(c, r, ii)), s(get(c, r, ni)));
                    }
                }
                if let Some(id) = b.name.strip_prefix("link_") {
                    link_blocks.insert(id.to_string(), b);
                } else if let Some(id) = b.name.strip_prefix("mod_") {
                    mod_blocks.insert(id.to_string(), b);
                }
            }
        }
        // `convert_list_block`: one entry per link block with restraint loops,
        // in block order of the merged list (geostd-only blocks after the rest)
        let has_restraints = |b: &cif::Block| {
            ["_chem_link_bond", "_chem_link_angle", "_chem_link_tor", "_chem_link_chir", "_chem_link_plane"]
                .iter()
                .any(|c| b.category(c).is_some())
        };
        let mut order_ids: Vec<String> = Vec::new();
        let geostd_restraints: FxHashSet<&str> =
            gdoc.blocks.iter().filter(|b| has_restraints(b)).filter_map(|b| b.name.strip_prefix("link_")).collect();
        for (doc, base) in [(&mdoc, true), (&gdoc, false)] {
            for b in &doc.blocks {
                let Some(id) = b.name.strip_prefix("link_") else { continue };
                let restrained = has_restraints(b) || (base && geostd_restraints.contains(id));
                if id != "list" && restrained && !order_ids.iter().any(|x| x == id) {
                    order_ids.push(id.to_string());
                }
            }
        }
        for doc in &extra_parsed {
            for b in &doc.blocks {
                if let Some(id) = b.name.strip_prefix("link_") {
                    if has_restraints(b) && !order_ids.iter().any(|x| x == id) {
                        order_ids.push(id.to_string());
                    }
                }
            }
        }
        for (id, b) in &link_blocks {
            let idx = match ml.link_index.get(id) {
                Some(&i) => i,
                None => {
                    ml.link_index.insert(id.clone(), ml.links.len());
                    ml.links.push(ChemLink { id: id.clone(), ..Default::default() });
                    ml.links.len() - 1
                }
            };
            parse_link_block(b, &mut ml.links[idx]);
        }
        ml.link_list_order = order_ids.iter().filter_map(|id| ml.link_index.get(id).copied()).collect();
        for (id, b) in &mod_blocks {
            let mut m = ChemMod { id: id.clone(), name: mod_names.get(id).cloned().unwrap_or_default(), ..Default::default() };
            parse_mod_block(b, &mut m);
            ml.mods.insert(id.clone(), m);
        }
        // energy library (geostd)
        let etext = read(root.join("geostd/ener_lib.cif"))?;
        let edoc = cif::parse(&etext);
        for b in &edoc.blocks {
            if let Some(c) = b.category("_lib_atom") {
                let cols: Vec<Option<usize>> = ["type", "hb_type", "vdw_radius", "vdwh_radius", "ion_radius", "element", "vdw_radius_neutron"]
                    .iter()
                    .map(|t| c.col(t))
                    .collect();
                for r in 0..c.nrows() {
                    let e = EnerAtom {
                        type_: s(get(c, r, cols[0])),
                        hb_type: s(get(c, r, cols[1])),
                        vdw_radius: cif::parse_f64(get(c, r, cols[2])),
                        vdwh_radius: cif::parse_f64(get(c, r, cols[3])),
                        ion_radius: cif::parse_f64(get(c, r, cols[4])),
                        element: s(get(c, r, cols[5])),
                        vdw_radius_neutron: cif::parse_f64(get(c, r, cols[6])),
                    };
                    ml.ener.insert(e.type_.clone(), e);
                }
            }
        }
        Ok(ml)
    }

    pub fn link(&self, id: &str) -> Option<&ChemLink> {
        self.link_index.get(id).map(|&i| &self.links[i])
    }

    /// Register dictionaries from a user-supplied restraints CIF.
    pub fn add_user_cif(&mut self, text: &str, path: &Path) {
        let doc = cif::parse(text);
        for comp in parse_comp_doc(&doc, path) {
            self.user_comps.insert(comp.id.to_ascii_uppercase(), Arc::new(comp));
        }
    }

    /// `get_comp_comp_id_direct`: geostd(std) > mon_lib(std) > geostd(id) > mon_lib(id).
    pub fn comp(&self, comp_id: &str) -> Option<Arc<Comp>> {
        let id = comp_id.trim().to_ascii_uppercase();
        if id.is_empty() {
            return None;
        }
        if let Some(c) = self.user_comps.get(&id) {
            return Some(c.clone());
        }
        if let Some(c) = self.comp_cache.read().unwrap().get(&id) {
            return c.clone();
        }
        let std = self.comp_synonyms.get(&id).cloned().unwrap_or_default();
        let mut found: Option<Arc<Comp>> = None;
        'outer: for case_insensitive in [false, true] {
            for trial in [std.as_str(), id.as_str()] {
                if trial.is_empty() {
                    continue;
                }
                let first = trial.chars().next().unwrap().to_ascii_lowercase().to_string();
                let cands = [
                    self.root.join("geostd").join(&first).join(format!("data_{}.cif", trial)),
                    self.root.join("mon_lib").join(&first).join(format!("{}.cif", mon_lib_file_name(trial))),
                ];
                for (k, p) in cands.iter().enumerate() {
                    let path = if case_insensitive { find_case_insensitive(p) } else if p.is_file() { Some(p.clone()) } else { None };
                    let Some(path) = path else { continue };
                    let _ = k;
                    if let Ok(t) = std::fs::read_to_string(&path) {
                        let doc = cif::parse(&t);
                        let comps = parse_comp_doc(&doc, &path);
                        let want_std = std.to_ascii_uppercase();
                        let pick = comps
                            .iter()
                            .find(|c| !want_std.is_empty() && c.id.to_ascii_uppercase() == want_std)
                            .or_else(|| comps.iter().find(|c| c.id.to_ascii_uppercase() == id))
                            .or_else(|| comps.first());
                        if let Some(c) = pick {
                            found = Some(Arc::new(c.clone()));
                            break 'outer;
                        }
                    }
                }
            }
        }
        self.comp_cache.write().unwrap().insert(id, found.clone());
        found
    }

    /// Restraints built from the CCD for a residue that neither library
    /// describes (Reduce2's fallback with `strict`, the fixed-mode one without;
    /// see [`crate::ccdrestraints::comp_from_ccd`]).
    pub fn ccd_comp(&self, comp_id: &str, strict: bool) -> Option<Arc<Comp>> {
        let id = comp_id.trim().to_ascii_uppercase();
        if id.is_empty() || id == "UNL" {
            return None;
        }
        let key = format!("{}{}", if strict { "strict:" } else { "" }, id);
        if let Some(c) = self.ccd_comp_cache.read().unwrap().get(&key) {
            return c.clone();
        }
        let first = id.chars().next().unwrap().to_ascii_lowercase().to_string();
        let p = self.root.join("chemical_components").join(first).join(format!("data_{}.cif", id));
        let comp = std::fs::read_to_string(&p)
            .ok()
            .and_then(|t| crate::ccdrestraints::comp_from_ccd(&t, &p, strict))
            .map(Arc::new);
        self.ccd_comp_cache.write().unwrap().insert(key, comp.clone());
        comp
    }

    /// The GeoStd neutron or low-pH variant of a residue dictionary
    /// (`get_comp_comp_id_direct(resname, pH_range=...)`, neutron first), which
    /// cctbx applies to residues whose hydrogens have no energy type.
    /// Ideal sites of `resname` in a restraint file that carries coordinates
    /// (`_dictionary_sites`), read once per file and residue name.
    pub fn dictionary_sites(&self, source: &Path, resname: &str) -> Option<Arc<FxHashMap<String, crate::geom::Vec3>>> {
        let key = (source.to_path_buf(), resname.to_string());
        if let Some(s) = self.sites_cache.read().unwrap().get(&key) {
            return s.clone();
        }
        let sites = read_dictionary_sites(source, resname).map(Arc::new);
        self.sites_cache.write().unwrap().insert(key, sites.clone());
        sites
    }

    pub fn ph_variant(&self, comp_id: &str) -> Option<Arc<Comp>> {
        let id = comp_id.trim().to_ascii_uppercase();
        if id.is_empty() {
            return None;
        }
        if let Some(c) = self.variant_cache.read().unwrap().get(&id) {
            return c.clone();
        }
        let first = id.chars().next().unwrap().to_ascii_lowercase().to_string();
        let dir = self.root.join("geostd").join(first);
        let mut found = None;
        for suffix in ["_neutron", "_pH_low"] {
            let p = dir.join(format!("data_{}{}.cif", id, suffix));
            let Ok(t) = std::fs::read_to_string(&p) else { continue };
            let doc = cif::parse(&t);
            found = parse_comp_doc(&doc, &p).into_iter().find(|c| c.id.eq_ignore_ascii_case(&id)).map(Arc::new);
            if found.is_some() {
                break;
            }
        }
        self.variant_cache.write().unwrap().insert(id, found.clone());
        found
    }

    /// The CCD entry for a residue name.
    pub fn ccd(&self, comp_id: &str) -> Option<Arc<CcdEntry>> {
        let id = comp_id.trim().to_ascii_uppercase();
        if id.is_empty() {
            return None;
        }
        if let Some(c) = self.ccd_cache.read().unwrap().get(&id) {
            return c.clone();
        }
        let first = id.chars().next().unwrap().to_ascii_lowercase().to_string();
        let p = self.root.join("chemical_components").join(first).join(format!("data_{}.cif", id));
        let entry = std::fs::read_to_string(&p).ok().and_then(|t| {
            let doc = cif::parse(&t);
            let b = doc.blocks.first()?;
            let mut e = CcdEntry { id: id.clone(), ..Default::default() };
            if let Some(c) = b.category("_chem_comp") {
                e.type_ = c.get_tag(0, "type").map(s).unwrap_or_default();
            }
            if let Some(c) = b.category("_chem_comp_atom") {
                let ai = c.col("atom_id");
                let ti = c.col("type_symbol");
                let mx = [c.col("model_Cartn_x"), c.col("model_Cartn_y"), c.col("model_Cartn_z")];
                let ix = [
                    c.col("pdbx_model_Cartn_x_ideal"),
                    c.col("pdbx_model_Cartn_y_ideal"),
                    c.col("pdbx_model_Cartn_z_ideal"),
                ];
                for r in 0..c.nrows() {
                    let xyz = |cols: &[Option<usize>; 3]| -> Option<[f64; 3]> {
                        let v: Vec<Option<f64>> = cols.iter().map(|k| cif::parse_f64(get(c, r, *k))).collect();
                        Some([v[0]?, v[1]?, v[2]?])
                    };
                    e.atoms.push((s(get(c, r, ai)), s(get(c, r, ti)).to_ascii_uppercase(), xyz(&mx), xyz(&ix)));
                }
            }
            if let Some(c) = b.category("_chem_comp_bond") {
                let (a1, a2, o) = (c.col("atom_id_1"), c.col("atom_id_2"), c.col("value_order"));
                for r in 0..c.nrows() {
                    e.bonds.push((s(get(c, r, a1)), s(get(c, r, a2)), s(get(c, r, o))));
                }
            }
            Some(Arc::new(e))
        });
        self.ccd_cache.write().unwrap().insert(id, entry.clone());
        entry
    }

    /// Apply a modification to a dictionary copy (`comp_comp_id.apply_mod`).
    pub fn apply_mod(&self, comp: &Comp, mod_id: &str) -> Option<Comp> {
        let m = self.mods.get(mod_id)?;
        if comp.applied_mods.iter().any(|x| x == mod_id) {
            return None;
        }
        let mut c = comp.clone();
        for ma in &m.atoms {
            match ma.function.as_str() {
                "add" => c.atoms.push(CompAtom {
                    id: ma.new_atom_id.clone(),
                    type_symbol: ma.new_type_symbol.clone(),
                    type_energy: if ma.new_type_energy.is_empty() { None } else { Some(ma.new_type_energy.clone()) },
                }),
                "delete" => {
                    let id = ma.atom_id.clone();
                    if c.has_atom(&id) {
                        c.atoms.retain(|a| a.id != id);
                        c.bonds.retain(|b| b.a1 != id && b.a2 != id);
                        c.angles.retain(|a| a.a1 != id && a.a2 != id && a.a3 != id);
                        c.tors.retain(|t| !t.a.contains(&id));
                        c.chirs.retain(|t| t.centre != id && !t.a.contains(&id));
                        c.planes.retain(|p| p.atom != id);
                    }
                }
                "change" => {
                    let id = ma.atom_id.clone();
                    let new_id = if ma.new_atom_id.is_empty() { id.clone() } else { ma.new_atom_id.clone() };
                    if let Some(a) = c.atoms.iter_mut().find(|a| a.id == id) {
                        if !ma.new_type_symbol.is_empty() {
                            a.type_symbol = ma.new_type_symbol.clone();
                        }
                        if !ma.new_type_energy.is_empty() {
                            a.type_energy = Some(ma.new_type_energy.clone());
                        }
                        a.id = new_id.clone();
                    }
                    if new_id != id {
                        let ren = |x: &mut String| {
                            if *x == id {
                                *x = new_id.clone();
                            }
                        };
                        for b in c.bonds.iter_mut() {
                            ren(&mut b.a1);
                            ren(&mut b.a2);
                        }
                        for a in c.angles.iter_mut() {
                            ren(&mut a.a1);
                            ren(&mut a.a2);
                            ren(&mut a.a3);
                        }
                        for t in c.tors.iter_mut() {
                            for x in t.a.iter_mut() {
                                ren(x);
                            }
                        }
                        for t in c.chirs.iter_mut() {
                            ren(&mut t.centre);
                            for x in t.a.iter_mut() {
                                ren(x);
                            }
                        }
                        for p in c.planes.iter_mut() {
                            ren(&mut p.atom);
                        }
                    }
                }
                _ => {}
            }
        }
        for mb in &m.bonds {
            let same = |b: &CompBond| (b.a1 == mb.a1 && b.a2 == mb.a2) || (b.a1 == mb.a2 && b.a2 == mb.a1);
            match mb.function.as_str() {
                "add" => c.bonds.push(CompBond {
                    a1: mb.a1.clone(),
                    a2: mb.a2.clone(),
                    type_: mb.new_type.clone(),
                    value_dist: mb.new_value_dist,
                    esd: mb.new_esd,
                    value_dist_neutron: mb.new_value_dist_neutron,
                }),
                "delete" => c.bonds.retain(|b| !same(b)),
                "change" => {
                    for b in c.bonds.iter_mut().filter(|b| same(b)) {
                        if !mb.new_type.is_empty() {
                            b.type_ = mb.new_type.clone();
                        }
                        if mb.new_value_dist.is_some() {
                            b.value_dist = mb.new_value_dist;
                        }
                        if mb.new_esd.is_some() {
                            b.esd = mb.new_esd;
                        }
                        if mb.new_value_dist_neutron.is_some() {
                            b.value_dist_neutron = mb.new_value_dist_neutron;
                        }
                    }
                }
                _ => {}
            }
        }
        for ma in &m.angles {
            let same = |a: &CompAngle| {
                a.a2 == ma.a[1] && ((a.a1 == ma.a[0] && a.a3 == ma.a[2]) || (a.a1 == ma.a[2] && a.a3 == ma.a[0]))
            };
            match ma.function.as_str() {
                "add" => c.angles.push(CompAngle {
                    a1: ma.a[0].clone(),
                    a2: ma.a[1].clone(),
                    a3: ma.a[2].clone(),
                    value: ma.new_value,
                    esd: ma.new_esd,
                }),
                "delete" => c.angles.retain(|a| !same(a)),
                "change" => {
                    for a in c.angles.iter_mut().filter(|a| same(a)) {
                        if ma.new_value.is_some() {
                            a.value = ma.new_value;
                        }
                        if ma.new_esd.is_some() {
                            a.esd = ma.new_esd;
                        }
                    }
                }
                _ => {}
            }
        }
        for mt in &m.tors {
            let same = |t: &CompTor| t.a == mt.a;
            match mt.function.as_str() {
                "add" => c.tors.push(CompTor {
                    id: mt.id.clone(),
                    a: mt.a.clone(),
                    value: mt.new_value,
                    esd: mt.new_esd,
                    period: mt.new_period.unwrap_or(0),
                    alt_values: mt.new_alt.clone(),
                }),
                "delete" => c.tors.retain(|t| !same(t)),
                "change" => {
                    for t in c.tors.iter_mut().filter(|t| same(t)) {
                        if mt.new_value.is_some() {
                            t.value = mt.new_value;
                        }
                        if mt.new_esd.is_some() {
                            t.esd = mt.new_esd;
                        }
                        if let Some(p) = mt.new_period {
                            t.period = p;
                        }
                        if mt.new_alt.is_some() {
                            t.alt_values = mt.new_alt.clone();
                        }
                    }
                }
                _ => {}
            }
        }
        for mc in &m.chirs {
            match mc.function.as_str() {
                "add" => c.chirs.push(CompChir {
                    id: mc.id.clone(),
                    centre: mc.centre.clone(),
                    a: mc.a.clone(),
                    volume_sign: mc.new_volume_sign.clone(),
                }),
                "delete" => c.chirs.retain(|x| !(x.centre == mc.centre && x.a == mc.a) && !(mc.id.len() > 0 && x.id == mc.id)),
                "change" => {
                    for x in c.chirs.iter_mut().filter(|x| x.centre == mc.centre || (!mc.id.is_empty() && x.id == mc.id)) {
                        if !mc.new_volume_sign.is_empty() {
                            x.volume_sign = mc.new_volume_sign.clone();
                        }
                    }
                }
                _ => {}
            }
        }
        for mp in &m.planes {
            match mp.function.as_str() {
                "add" => c.planes.push(CompPlaneAtom { plane_id: mp.plane_id.clone(), atom: mp.atom.clone(), esd: mp.new_esd }),
                "delete" => c.planes.retain(|p| !(p.plane_id == mp.plane_id && p.atom == mp.atom)),
                "change" => {
                    for p in c.planes.iter_mut().filter(|p| p.plane_id == mp.plane_id && p.atom == mp.atom) {
                        if mp.new_esd.is_some() {
                            p.esd = mp.new_esd;
                        }
                    }
                }
                _ => {}
            }
        }
        c.applied_mods.push(mod_id.to_string());
        if m.name.contains("terminus") {
            c.is_terminus = true;
        }
        Some(c)
    }
}

/// mon_lib file name for Windows device names (CON -> CON_CON etc.).
fn mon_lib_file_name(id: &str) -> String {
    let dev = ["CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4", "LPT1", "LPT2", "LPT3"];
    if dev.contains(&id) { format!("{}_{}", id, id) } else { id.to_string() }
}

fn find_case_insensitive(p: &Path) -> Option<PathBuf> {
    let dir = p.parent()?;
    let want = p.file_name()?.to_string_lossy().to_ascii_lowercase();
    for e in std::fs::read_dir(dir).ok()?.flatten() {
        if e.file_name().to_string_lossy().to_ascii_lowercase() == want {
            return Some(e.path());
        }
    }
    None
}

fn alt_list(v: &str) -> Option<String> {
    if cif::is_null(v) { None } else { Some(v.to_string()) }
}

/// Parse all `comp_*` blocks of a dictionary file.
pub fn parse_comp_doc(doc: &cif::Document, path: &Path) -> Vec<Comp> {
    let mut groups: FxHashMap<String, String> = FxHashMap::default();
    for b in &doc.blocks {
        if let Some(c) = b.category("_chem_comp") {
            let (ii, gi) = (c.col("id"), c.col("group"));
            for r in 0..c.nrows() {
                groups.insert(s(get(c, r, ii)).to_ascii_uppercase(), s(get(c, r, gi)));
            }
        }
    }
    let mut out = Vec::new();
    for b in &doc.blocks {
        let Some(name) = b.name.strip_prefix("comp_") else { continue };
        if name == "list" {
            continue;
        }
        let mut comp = Comp {
            id: name.to_string(),
            group: groups.get(&name.to_ascii_uppercase()).cloned().unwrap_or_default(),
            source: path.to_path_buf(),
            ..Default::default()
        };
        if let Some(c) = b.category("_chem_comp_atom") {
            let (ai, ti, ei) = (c.col("atom_id"), c.col("type_symbol"), c.col("type_energy"));
            for r in 0..c.nrows() {
                let te = get(c, r, ei);
                comp.atoms.push(CompAtom {
                    id: s(get(c, r, ai)),
                    type_symbol: s(get(c, r, ti)).to_ascii_uppercase(),
                    type_energy: if cif::is_null(te) { None } else { Some(te.to_string()) },
                });
            }
        }
        if let Some(c) = b.category("_chem_comp_bond") {
            let cols: Vec<Option<usize>> = ["atom_id_1", "atom_id_2", "type", "value_dist", "value_dist_esd", "value_dist_neutron"]
                .iter()
                .map(|t| c.col(t))
                .collect();
            let nuc = c.col("value_dist_nucleus");
            for r in 0..c.nrows() {
                let neutron = opt_f(Some(get(c, r, cols[5]))).or_else(|| opt_f(Some(get(c, r, nuc))));
                comp.bonds.push(CompBond {
                    a1: s(get(c, r, cols[0])),
                    a2: s(get(c, r, cols[1])),
                    type_: s(get(c, r, cols[2])),
                    value_dist: opt_f(Some(get(c, r, cols[3]))),
                    esd: opt_f(Some(get(c, r, cols[4]))),
                    value_dist_neutron: neutron,
                });
            }
        }
        if let Some(c) = b.category("_chem_comp_angle") {
            let cols: Vec<Option<usize>> =
                ["atom_id_1", "atom_id_2", "atom_id_3", "value_angle", "value_angle_esd"].iter().map(|t| c.col(t)).collect();
            for r in 0..c.nrows() {
                comp.angles.push(CompAngle {
                    a1: s(get(c, r, cols[0])),
                    a2: s(get(c, r, cols[1])),
                    a3: s(get(c, r, cols[2])),
                    value: opt_f(Some(get(c, r, cols[3]))),
                    esd: opt_f(Some(get(c, r, cols[4]))),
                });
            }
        }
        if let Some(c) = b.category("_chem_comp_tor") {
            let cols: Vec<Option<usize>> = [
                "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4", "value_angle", "value_angle_esd", "period", "alt_value_angle",
            ]
            .iter()
            .map(|t| c.col(t))
            .collect();
            for r in 0..c.nrows() {
                comp.tors.push(CompTor {
                    id: s(get(c, r, cols[0])),
                    a: [s(get(c, r, cols[1])), s(get(c, r, cols[2])), s(get(c, r, cols[3])), s(get(c, r, cols[4]))],
                    value: opt_f(Some(get(c, r, cols[5]))),
                    esd: opt_f(Some(get(c, r, cols[6]))),
                    period: get(c, r, cols[7]).trim().parse::<i32>().unwrap_or(0),
                    alt_values: alt_list(get(c, r, cols[8])),
                });
            }
        }
        if let Some(c) = b.category("_chem_comp_chir") {
            let cols: Vec<Option<usize>> =
                ["id", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "volume_sign"].iter().map(|t| c.col(t)).collect();
            for r in 0..c.nrows() {
                comp.chirs.push(CompChir {
                    id: s(get(c, r, cols[0])),
                    centre: s(get(c, r, cols[1])),
                    a: [s(get(c, r, cols[2])), s(get(c, r, cols[3])), s(get(c, r, cols[4]))],
                    volume_sign: s(get(c, r, cols[5])),
                });
            }
        }
        if let Some(c) = b.category("_chem_comp_plane_atom") {
            let cols: Vec<Option<usize>> = ["plane_id", "atom_id", "dist_esd"].iter().map(|t| c.col(t)).collect();
            for r in 0..c.nrows() {
                comp.planes.push(CompPlaneAtom {
                    plane_id: s(get(c, r, cols[0])),
                    atom: s(get(c, r, cols[1])),
                    esd: opt_f(Some(get(c, r, cols[2]))),
                });
            }
        }
        out.push(comp);
    }
    out
}

fn parse_link_block(b: &cif::Block, l: &mut ChemLink) {
    if let Some(c) = b.category("_chem_link_bond") {
        let cols: Vec<Option<usize>> =
            ["atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2", "value_dist", "value_dist_esd"].iter().map(|t| c.col(t)).collect();
        l.bonds.clear();
        for r in 0..c.nrows() {
            l.bonds.push(LinkBond {
                c1: comp_digit(get(c, r, cols[0])),
                a1: s(get(c, r, cols[1])),
                c2: comp_digit(get(c, r, cols[2])),
                a2: s(get(c, r, cols[3])),
                value_dist: opt_f(Some(get(c, r, cols[4]))),
                esd: opt_f(Some(get(c, r, cols[5]))),
            });
        }
    }
    if let Some(c) = b.category("_chem_link_angle") {
        let cols: Vec<Option<usize>> = [
            "atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2", "atom_3_comp_id", "atom_id_3", "value_angle", "value_angle_esd",
        ]
        .iter()
        .map(|t| c.col(t))
        .collect();
        l.angles.clear();
        for r in 0..c.nrows() {
            l.angles.push(LinkAngle {
                c: [comp_digit(get(c, r, cols[0])), comp_digit(get(c, r, cols[2])), comp_digit(get(c, r, cols[4]))],
                a: [s(get(c, r, cols[1])), s(get(c, r, cols[3])), s(get(c, r, cols[5]))],
                value: opt_f(Some(get(c, r, cols[6]))),
                esd: opt_f(Some(get(c, r, cols[7]))),
            });
        }
    }
    if let Some(c) = b.category("_chem_link_tor") {
        let cols: Vec<Option<usize>> = [
            "id", "atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2", "atom_3_comp_id", "atom_id_3", "atom_4_comp_id",
            "atom_id_4", "value_angle", "value_angle_esd", "period",
        ]
        .iter()
        .map(|t| c.col(t))
        .collect();
        l.tors.clear();
        for r in 0..c.nrows() {
            l.tors.push(LinkTor {
                id: s(get(c, r, cols[0])),
                c: [
                    comp_digit(get(c, r, cols[1])),
                    comp_digit(get(c, r, cols[3])),
                    comp_digit(get(c, r, cols[5])),
                    comp_digit(get(c, r, cols[7])),
                ],
                a: [s(get(c, r, cols[2])), s(get(c, r, cols[4])), s(get(c, r, cols[6])), s(get(c, r, cols[8]))],
                value: opt_f(Some(get(c, r, cols[9]))),
                esd: opt_f(Some(get(c, r, cols[10]))),
                period: get(c, r, cols[11]).trim().parse::<i32>().unwrap_or(0),
            });
        }
    }
    if let Some(c) = b.category("_chem_link_chir") {
        let cols: Vec<Option<usize>> = [
            "atom_centre_comp_id", "atom_id_centre", "atom_1_comp_id", "atom_id_1", "atom_2_comp_id", "atom_id_2", "atom_3_comp_id",
            "atom_id_3", "volume_sign",
        ]
        .iter()
        .map(|t| c.col(t))
        .collect();
        l.chirs.clear();
        for r in 0..c.nrows() {
            l.chirs.push(LinkChir {
                c: [
                    comp_digit(get(c, r, cols[0])),
                    comp_digit(get(c, r, cols[2])),
                    comp_digit(get(c, r, cols[4])),
                    comp_digit(get(c, r, cols[6])),
                ],
                a: [s(get(c, r, cols[1])), s(get(c, r, cols[3])), s(get(c, r, cols[5])), s(get(c, r, cols[7]))],
                volume_sign: s(get(c, r, cols[8])),
            });
        }
    }
    if let Some(c) = b.category("_chem_link_plane") {
        let cols: Vec<Option<usize>> = ["plane_id", "atom_comp_id", "atom_id", "dist_esd"].iter().map(|t| c.col(t)).collect();
        l.planes.clear();
        for r in 0..c.nrows() {
            l.planes.push(LinkPlaneAtom {
                plane_id: s(get(c, r, cols[0])),
                c: comp_digit(get(c, r, cols[1])),
                atom: s(get(c, r, cols[2])),
                esd: opt_f(Some(get(c, r, cols[3]))),
            });
        }
    }
}

fn parse_mod_block(b: &cif::Block, m: &mut ChemMod) {
    if let Some(c) = b.category("_chem_mod_atom") {
        let cols: Vec<Option<usize>> =
            ["function", "atom_id", "new_atom_id", "new_type_symbol", "new_type_energy"].iter().map(|t| c.col(t)).collect();
        for r in 0..c.nrows() {
            m.atoms.push(ModAtom {
                function: s(get(c, r, cols[0])),
                atom_id: s(get(c, r, cols[1])),
                new_atom_id: s(get(c, r, cols[2])),
                new_type_symbol: s(get(c, r, cols[3])).to_ascii_uppercase(),
                new_type_energy: s(get(c, r, cols[4])),
            });
        }
    }
    if let Some(c) = b.category("_chem_mod_bond") {
        let cols: Vec<Option<usize>> =
            ["function", "atom_id_1", "atom_id_2", "new_type", "new_value_dist", "new_value_dist_esd", "new_value_dist_neutron"]
                .iter()
                .map(|t| c.col(t))
                .collect();
        for r in 0..c.nrows() {
            m.bonds.push(ModBond {
                function: s(get(c, r, cols[0])),
                a1: s(get(c, r, cols[1])),
                a2: s(get(c, r, cols[2])),
                new_type: s(get(c, r, cols[3])),
                new_value_dist: opt_f(Some(get(c, r, cols[4]))),
                new_esd: opt_f(Some(get(c, r, cols[5]))),
                new_value_dist_neutron: opt_f(Some(get(c, r, cols[6]))),
            });
        }
    }
    if let Some(c) = b.category("_chem_mod_angle") {
        let cols: Vec<Option<usize>> =
            ["function", "atom_id_1", "atom_id_2", "atom_id_3", "new_value_angle", "new_value_angle_esd"].iter().map(|t| c.col(t)).collect();
        for r in 0..c.nrows() {
            m.angles.push(ModAngle {
                function: s(get(c, r, cols[0])),
                a: [s(get(c, r, cols[1])), s(get(c, r, cols[2])), s(get(c, r, cols[3]))],
                new_value: opt_f(Some(get(c, r, cols[4]))),
                new_esd: opt_f(Some(get(c, r, cols[5]))),
            });
        }
    }
    if let Some(c) = b.category("_chem_mod_tor") {
        let cols: Vec<Option<usize>> = [
            "function", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4", "new_value_angle", "new_value_angle_esd", "new_period",
            "new_alt_value_angle",
        ]
        .iter()
        .map(|t| c.col(t))
        .collect();
        for r in 0..c.nrows() {
            m.tors.push(ModTor {
                function: s(get(c, r, cols[0])),
                id: s(get(c, r, cols[1])),
                a: [s(get(c, r, cols[2])), s(get(c, r, cols[3])), s(get(c, r, cols[4])), s(get(c, r, cols[5]))],
                new_value: opt_f(Some(get(c, r, cols[6]))),
                new_esd: opt_f(Some(get(c, r, cols[7]))),
                new_period: get(c, r, cols[8]).trim().parse::<i32>().ok(),
                new_alt: alt_list(get(c, r, cols[9])),
            });
        }
    }
    if let Some(c) = b.category("_chem_mod_chir") {
        let cols: Vec<Option<usize>> =
            ["function", "id", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "new_volume_sign"].iter().map(|t| c.col(t)).collect();
        for r in 0..c.nrows() {
            m.chirs.push(ModChir {
                function: s(get(c, r, cols[0])),
                id: s(get(c, r, cols[1])),
                centre: s(get(c, r, cols[2])),
                a: [s(get(c, r, cols[3])), s(get(c, r, cols[4])), s(get(c, r, cols[5]))],
                new_volume_sign: s(get(c, r, cols[6])),
            });
        }
    }
    if let Some(c) = b.category("_chem_mod_plane_atom") {
        let cols: Vec<Option<usize>> = ["function", "plane_id", "atom_id", "new_dist_esd"].iter().map(|t| c.col(t)).collect();
        for r in 0..c.nrows() {
            m.planes.push(ModPlane {
                function: s(get(c, r, cols[0])),
                plane_id: s(get(c, r, cols[1])),
                atom: s(get(c, r, cols[2])),
                new_esd: opt_f(Some(get(c, r, cols[3]))),
            });
        }
    }
}

fn read_dictionary_sites(source: &Path, resname: &str) -> Option<FxHashMap<String, crate::geom::Vec3>> {
    let text = std::fs::read_to_string(source).ok()?;
    let doc = cif::parse(&text);
    for b in &doc.blocks {
        let Some(cat) = b.category("_chem_comp_atom") else { continue };
        let (Some(ai), Some(xi), Some(yi), Some(zi)) = (cat.col("atom_id"), cat.col("x"), cat.col("y"), cat.col("z")) else { continue };
        let ci = cat.col("comp_id");
        if ci.is_none() && b.name != format!("comp_{}", resname) {
            continue;
        }
        let mut sites = FxHashMap::default();
        for r in 0..cat.nrows() {
            if let Some(ci) = ci {
                if cat.get(r, ci).trim() != resname {
                    continue;
                }
            }
            let (x, y, z) = (cif::parse_f64(cat.get(r, xi)), cif::parse_f64(cat.get(r, yi)), cif::parse_f64(cat.get(r, zi)));
            if let (Some(x), Some(y), Some(z)) = (x, y, z) {
                sites.insert(cat.get(r, ai).trim_matches('"').to_string(), crate::geom::v3(x, y, z));
            }
        }
        if !sites.is_empty() {
            return Some(sites);
        }
    }
    None
}
