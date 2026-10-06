//! mmCIF reading and writing with the conventions of iotbx
//! (`iotbx.pdb.mmcif.pdb_hierarchy_builder` and `hierarchy.as_cif_block`).

use crate::cell::{space_group_info, UnitCell};
use crate::cif;
use crate::cifsource::{CifCell, CifSink, CifSource, CifTable, CifText};
use crate::geom::v3;
use crate::model::*;
use crate::resclass::{self, ResClass};
use std::fmt::Write as _;

/// `format_pdb_atom_name`: the 4-column PDB spelling of an mmCIF atom name.
fn format_pdb_atom_name(name: &str, element: &str) -> String {
    if name.len() >= 4 {
        return name.to_string();
    }
    let mut n = name.trim().to_string();
    let first = n.chars().next().map(|c| c.to_ascii_uppercase());
    if element.len() == 1 && first == element.chars().next().map(|c| c.to_ascii_uppercase()) {
        n = format!(" {}", n);
    }
    while n.len() < 4 {
        n.push(' ');
    }
    n
}

fn is_aa_or_rna_dna(resname: &str) -> bool {
    matches!(
        resclass::get_class(resname.trim()),
        ResClass::CommonAminoAcid | ResClass::ModifiedAminoAcid | ResClass::CommonRnaDna | ResClass::ModifiedRnaDna
    )
}

/// Read the first data block with an `_atom_site` loop.
pub fn read_mmcif(text: &str) -> Result<Structure, String> {
    let doc = cif::parse(text);
    let block = doc
        .blocks
        .iter()
        .find(|b| b.category("_atom_site").is_some())
        .ok_or_else(|| "no _atom_site loop found in the mmCIF file".to_string())?;
    structure_from_cif(block)
}

/// Build the model from a parsed mmCIF data block (`_atom_site`,
/// `_atom_site_anisotrop`, `_cell` and the space group items).
pub fn structure_from_cif<S: CifSource + ?Sized>(block: &S) -> Result<Structure, String> {
    let cat = block.table("atom_site").ok_or_else(|| "no _atom_site loop found in the mmCIF data".to_string())?;
    let col = |names: &[&str]| names.iter().find_map(|n| cat.column(n));
    let req = |names: &[&str]| col(names).ok_or_else(|| format!("_atom_site.{} is missing", names[0]));
    let c_type = req(&["type_symbol"])?;
    let c_name = req(&["label_atom_id", "auth_atom_id"])?;
    let c_alt = req(&["label_alt_id"])?;
    let c_auth_asym = col(&["auth_asym_id", "label_asym_id"]).ok_or("_atom_site.auth_asym_id is missing")?;
    let c_label_asym = col(&["label_asym_id", "auth_asym_id"]).unwrap_or(c_auth_asym);
    let c_comp = req(&["auth_comp_id", "label_comp_id"])?;
    let c_seq = req(&["auth_seq_id", "label_seq_id"])?;
    let c_ins = col(&["pdbx_PDB_ins_code"]);
    let c_model = col(&["pdbx_PDB_model_num"]);
    let c_id = col(&["id"]);
    let c_group = col(&["group_PDB"]);
    let (c_x, c_y, c_z) = (req(&["Cartn_x"])?, req(&["Cartn_y"])?, req(&["Cartn_z"])?);
    let c_occ = req(&["occupancy"])?;
    let c_b = req(&["B_iso_or_equiv"])?;
    let c_charge = col(&["pdbx_formal_charge"]);
    let c_segid = col(&["auth_segid"]);
    let c_break = col(&["auth_break"]);
    let get = |r: usize, c: usize| cat.cell(r, c);
    let num = |r: usize, c: usize| -> Result<f64, String> {
        cat.number(r, c).ok_or_else(|| format!("bad number '{}' in _atom_site row {}", cat.cell(r, c), r + 1))
    };
    // anisotropic displacements by atom id
    let mut aniso: rustc_hash::FxHashMap<String, [f64; 6]> = rustc_hash::FxHashMap::default();
    if let Some(an) = block.table("atom_site_anisotrop") {
        if let Some(ci) = an.column("id") {
            let u: Vec<Option<usize>> =
                ["U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]"].iter().map(|n| an.column(n)).collect();
            let b: Vec<Option<usize>> =
                ["B[1][1]", "B[2][2]", "B[3][3]", "B[1][2]", "B[1][3]", "B[2][3]"].iter().map(|n| an.column(n)).collect();
            let (cols, scale) = if u.iter().all(|x| x.is_some()) {
                (u, 1.0)
            } else if b.iter().all(|x| x.is_some()) {
                (b, 1.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI))
            } else {
                (vec![], 1.0)
            };
            if !cols.is_empty() {
                for r in 0..an.row_count() {
                    let mut v = [0.0; 6];
                    let mut ok = true;
                    for k in 0..6 {
                        match cols[k].and_then(|c| an.number(r, c)) {
                            Some(x) => v[k] = x * scale,
                            None => ok = false,
                        }
                    }
                    if ok {
                        aniso.insert(an.cell(r, ci).into_owned(), v);
                    }
                }
            }
        }
    }

    let mut models: Vec<Model> = Vec::new();
    let (mut cur_model, mut cur_label_asym, mut cur_auth_asym, mut cur_res, mut cur_ins): (
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
    ) = (None, None, None, None, None);
    let mut prev_comp: Option<String> = None;
    let mut first_in_chain = true;
    // (altloc, resname) -> atom group index, for the current residue group
    let mut ag_index: Vec<(String, String, usize)> = Vec::new();
    for r in 0..cat.row_count() {
        let model_id = c_model.map(|c| get(r, c).into_owned()).unwrap_or_else(|| "1".into());
        let new_model = cur_model.as_deref() != Some(model_id.as_str());
        if new_model {
            models.push(Model { id: model_id.clone(), chains: Vec::new() });
        }
        let label_asym = get(r, c_label_asym).into_owned();
        let auth_asym = get(r, c_auth_asym).into_owned();
        if matches!(auth_asym.as_str(), "." | "?" | " " | "") {
            return Err("mmCIF file contains a record with an empty auth_asym_id".into());
        }
        let label_changed = cur_label_asym.as_deref() != Some(label_asym.as_str());
        let auth_changed = cur_auth_asym.as_deref() != Some(auth_asym.as_str());
        let model = models.last_mut().unwrap();
        let mut new_chain = false;
        if auth_changed || new_model || (label_changed && r > 0 && prev_comp.as_deref().map_or(false, is_aa_or_rna_dna)) {
            model.chains.push(Chain { id: auth_asym.clone(), residue_groups: Vec::new() });
            first_in_chain = true;
            new_chain = true;
        }
        let res = get(r, c_seq).into_owned();
        let mut ins = c_ins.map(|c| get(r, c).into_owned()).unwrap_or_else(|| " ".into());
        if matches!(ins.as_str(), "?" | ".") {
            ins = " ".into();
        }
        let chain = model.chains.last_mut().unwrap();
        if new_chain
            || cur_res.as_deref() != Some(res.as_str())
            || cur_ins.as_deref() != Some(ins.as_str())
            || auth_changed
            || new_model
        {
            let resseq = match res.trim().parse::<i64>() {
                Ok(v) => hy36_encode(4, v).unwrap_or_else(|| res.clone()),
                Err(_) => {
                    if res.len() != 4 {
                        return Err(format!("bad residue number '{}' in _atom_site row {}", res, r + 1));
                    }
                    res.clone()
                }
            };
            let mut link = true;
            if let Some(cb) = c_break {
                if !first_in_chain && get(r, cb) == "1" {
                    link = false;
                }
            }
            chain.residue_groups.push(ResidueGroup {
                resseq,
                icode: ins.clone(),
                link_to_previous: link,
                atom_groups: Vec::new(),
            });
            first_in_chain = false;
            ag_index.clear();
        }
        let rg = chain.residue_groups.last_mut().unwrap();
        let mut alt = get(r, c_alt).into_owned();
        if alt == "." || alt == "?" {
            alt.clear();
        }
        let comp = get(r, c_comp).into_owned();
        let gi = match ag_index.iter().find(|(a, c, _)| *a == alt && *c == comp) {
            Some(&(_, _, i)) => i,
            None => {
                let ag = AtomGroup { altloc: alt.clone(), resname: format!("{:>3}", comp), atoms: Vec::new() };
                if alt.is_empty() {
                    rg.atom_groups.insert(0, ag);
                    for e in ag_index.iter_mut() {
                        e.2 += 1;
                    }
                    ag_index.push((alt.clone(), comp.clone(), 0));
                    0
                } else {
                    rg.atom_groups.push(ag);
                    ag_index.push((alt.clone(), comp.clone(), rg.atom_groups.len() - 1));
                    rg.atom_groups.len() - 1
                }
            }
        };
        let type_symbol = get(r, c_type);
        let mut atom =
            Atom::new(&format_pdb_atom_name(&get(r, c_name), &type_symbol), &type_symbol, v3(num(r, c_x)?, num(r, c_y)?, num(r, c_z)?));
        atom.occ = num(r, c_occ)?;
        atom.b = num(r, c_b)?;
        atom.src = r as u32;
        if let Some(ci) = c_id {
            let id = get(r, ci);
            atom.serial = id.trim().parse::<i64>().ok().and_then(|v| hy36_encode(5, v)).unwrap_or_else(|| id.to_string());
            if let Some(u) = aniso.get(&*id) {
                atom.uij = Some(*u);
            }
        }
        atom.segid = match c_segid {
            Some(cs) => {
                let s = get(r, cs);
                let s = if cif::is_null(&s) { "" } else { &*s };
                format!("{:<4}", &s[..s.len().min(4)])
            }
            None => "    ".into(),
        };
        atom.hetero = c_group.map_or(false, |c| get(r, c) == "HETATM");
        if let Some(cc) = c_charge {
            let ch = get(r, cc);
            if !matches!(&*ch, "?" | ".") {
                let sign = if ch.ends_with('-') || ch.starts_with('-') { "-" } else { "+" };
                if let Ok(v) = ch.trim_matches(|c| c == ' ' || c == '-' || c == '+').parse::<i32>() {
                    atom.charge = if v == 0 { "  ".into() } else { format!("{}{}", v, sign) };
                }
            }
        }
        rg.atom_groups[gi].atoms.push(atom);
        cur_model = Some(model_id);
        cur_label_asym = Some(label_asym);
        cur_auth_asym = Some(auth_asym);
        cur_res = Some(res);
        cur_ins = Some(ins);
        prev_comp = Some(comp);
    }
    if models.len() == 1 {
        models[0].id = String::new();
    }
    // crystal symmetry
    let item = |category: &str, tag: &str| -> Option<String> {
        let t = block.table(category)?;
        let c = t.column(tag)?;
        if t.row_count() == 0 { None } else { Some(t.cell(0, c).into_owned()) }
    };
    let cellv = |tag: &str| -> Option<f64> {
        let t = block.table("cell")?;
        let c = t.column(tag)?;
        if t.row_count() == 0 { None } else { t.number(0, c) }
    };
    let cell = [
        cellv("length_a"),
        cellv("length_b"),
        cellv("length_c"),
        cellv("angle_alpha"),
        cellv("angle_beta"),
        cellv("angle_gamma"),
    ];
    let mut crystal = None;
    if let [Some(a), Some(b), Some(c), Some(al), Some(be), Some(ga)] = cell {
        let sg = [("symmetry", "space_group_name_H-M"), ("space_group", "name_H-M_alt"), ("space_group", "name_H-M_full")]
            .iter()
            .find_map(|(c, t)| item(c, t).filter(|v| !cif::is_null(v)))
            .unwrap_or_else(|| "P 1".to_string());
        crystal = Some(CrystalSymmetry { cell: [a, b, c, al, be, ga], space_group: sg });
    }
    let had = crystal.is_some();
    let mut st = Structure { models, crystal, had_cell_record: had, links: Vec::new(), pdb_records: PdbRecords::default() };
    st.reset_i_seq();
    Ok(st)
}

/// `all_label_asym_ids()[n]`: A..Z, then AA, AB, ... (letters with repeats).
fn label_asym_id(mut n: usize) -> String {
    let mut len = 1;
    let mut count = 26usize;
    while n >= count {
        n -= count;
        len += 1;
        count *= 26;
    }
    let mut s = vec![b'A'; len];
    for k in (0..len).rev() {
        s[k] = b'A' + (n % 26) as u8;
        n /= 26;
    }
    String::from_utf8(s).unwrap()
}

/// `_residue_group_kinds`: 'p'olymer, 'w'ater or 'l'igand per residue group.
fn residue_group_kinds(chain: &Chain) -> Vec<char> {
    let atom = |rg: &ResidueGroup, name: &str| -> Option<Vec3> {
        rg.atom_groups.iter().find_map(|ag| ag.atoms.iter().find(|a| a.name.trim() == name).map(|a| a.xyz))
    };
    let linked = |a: &ResidueGroup, b: &ResidueGroup| -> bool {
        for (p, n) in [("C", "N"), ("O3'", "P")] {
            if let (Some(x), Some(y)) = (atom(a, p), atom(b, n)) {
                if x.dist(y) < 2.0 {
                    return true;
                }
            }
        }
        false
    };
    let rgs = &chain.residue_groups;
    (0..rgs.len())
        .map(|i| {
            let cls = resclass::get_class(rgs[i].atom_groups[0].resname.trim());
            if matches!(
                cls,
                ResClass::CommonAminoAcid | ResClass::ModifiedAminoAcid | ResClass::DAminoAcid | ResClass::CommonRnaDna | ResClass::ModifiedRnaDna
            ) {
                'p'
            } else if cls == ResClass::CommonWater {
                'w'
            } else if (i > 0 && linked(&rgs[i - 1], &rgs[i])) || (i + 1 < rgs.len() && linked(&rgs[i], &rgs[i + 1])) {
                'p'
            } else {
                'l'
            }
        })
        .collect()
}

use crate::geom::Vec3;

/// Write the model as mmCIF the way `model_as_mmcif` lays it out.
pub fn write_mmcif(st: &Structure) -> String {
    let mut w = CifText::with_capacity(st.atoms_size() * 100 + 4096);
    write_cif(st, &mut w);
    w.out
}

/// Every atom with its model, chain and residue group indices and its serial
/// number in the output (serials restart at 1 in each model).
fn output_atoms(st: &Structure) -> impl Iterator<Item = (usize, usize, usize, &ResidueGroup, &AtomGroup, &Atom, i64)> {
    st.models.iter().enumerate().flat_map(|(mi, m)| {
        let mut serial = 0i64;
        m.chains.iter().enumerate().flat_map(move |(ci, c)| {
            c.residue_groups
                .iter()
                .enumerate()
                .flat_map(move |(ri, rg)| rg.atom_groups.iter().flat_map(move |ag| ag.atoms.iter().map(move |a| (mi, ci, ri, rg, ag, a))))
        })
        .map(move |(mi, ci, ri, rg, ag, a)| {
            serial += 1;
            (mi, ci, ri, rg, ag, a, serial)
        })
    })
}

/// Send the model to `sink` as one data block, with the same items, loops and
/// values that `write_mmcif` writes.
pub fn write_cif<K: CifSink + ?Sized>(st: &Structure, sink: &mut K) {
    sink.begin_block("default");
    let mut symops: &[&str] = &[];
    if let Some(cs) = &st.crystal {
        let sg = space_group_info(&cs.space_group, &cs.cell);
        let uc = crate::cell::cryst1_rotations(&cs.space_group, &cs.cell)
            .and_then(|r| UnitCell::new(cs.cell).and_then(|u| u.averaged(r)))
            .or_else(|| UnitCell::new(cs.cell));
        if let Some(uc) = &uc {
            let p = uc.params;
            for (k, v) in ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"].iter().zip(p.iter()) {
                sink.item(&format!("_cell.{}", k), CifCell::Text(&format!("{:.3}", v)));
            }
            sink.item("_cell.volume", CifCell::Text(&format!("{:.3}", uc.volume)));
        }
        if let Some(sg) = &sg {
            let number = sg.number.to_string();
            sink.item("_space_group.crystal_system", CifCell::Text(sg.crystal_system));
            sink.item("_space_group.IT_number", CifCell::Text(&number));
            sink.item("_space_group.name_H-M_alt", CifCell::Text(sg.hm));
            sink.item("_space_group.name_Hall", CifCell::Text(sg.hall));
            sink.item("_symmetry.space_group_name_H-M", CifCell::Text(sg.hm));
            sink.item("_symmetry.space_group_name_Hall", CifCell::Text(sg.hall));
            sink.item("_symmetry.Int_Tables_number", CifCell::Text(&number));
            symops = sg.ops;
        }
    }
    if !symops.is_empty() {
        sink.begin_loop(&["_space_group_symop.id", "_space_group_symop.operation_xyz"], symops.len());
        for (k, op) in symops.iter().enumerate() {
            sink.row(&[CifCell::Text(&(k + 1).to_string()), CifCell::Text(op)]);
        }
        sink.end_loop();
    }
    // label_asym_id per residue group and label_seq_id per atom group
    let mut lai: Vec<Vec<Vec<String>>> = Vec::new();
    let mut lsi: Vec<Vec<Vec<String>>> = Vec::new();
    let mut n_lai = 0usize;
    let mut prev_key = String::new();
    for m in &st.models {
        let mut ml = Vec::new();
        let mut ms = Vec::new();
        for c in &m.chains {
            let kinds = residue_group_kinds(c);
            let mut prev: Option<char> = None;
            let mut cl = Vec::new();
            let mut cs = Vec::new();
            let mut seq = 0usize;
            for (rg, &kind) in c.residue_groups.iter().zip(kinds.iter()) {
                if kind == 'l' {
                    if prev.is_some() {
                        n_lai += 1;
                    }
                } else if prev.is_some() && prev != Some(kind) {
                    n_lai += 1;
                }
                cl.push(label_asym_id(n_lai));
                prev = Some(kind);
                let key = format!("{}{}{}", c.id, rg.resseq, rg.icode);
                if key != prev_key {
                    seq += 1;
                    prev_key = key;
                }
                cs.push(if kind == 'p' { seq.to_string() } else { ".".to_string() });
            }
            n_lai += 1;
            ml.push(cl);
            ms.push(cs);
        }
        n_lai += 1;
        lai.push(ml);
        lsi.push(ms);
    }
    let mut struct_asym: Vec<String> = Vec::new();
    let mut comps: Vec<String> = Vec::new();
    for (mi, m) in st.models.iter().enumerate() {
        for (ci, c) in m.chains.iter().enumerate() {
            for (ri, rg) in c.residue_groups.iter().enumerate() {
                let l = &lai[mi][ci][ri];
                if !rg.atom_groups.iter().all(|g| g.atoms.is_empty()) && !struct_asym.contains(l) {
                    struct_asym.push(l.clone());
                }
                for ag in &rg.atom_groups {
                    let id = ag.resname.trim().to_string();
                    if !ag.atoms.is_empty() && !comps.contains(&id) {
                        comps.push(id);
                    }
                }
            }
        }
    }
    comps.sort();
    sink.begin_loop(&["_struct_asym.id"], struct_asym.len());
    for s in &struct_asym {
        sink.row(&[CifCell::Text(s)]);
    }
    sink.end_loop();
    sink.begin_loop(&["_chem_comp.id"], comps.len());
    for c in &comps {
        sink.row(&[CifCell::Text(c)]);
    }
    sink.end_loop();

    let model_nums: Vec<String> =
        st.models.iter().map(|m| if m.id.trim().is_empty() { "1".to_string() } else { m.id.trim().to_string() }).collect();
    let auth_asyms: Vec<Vec<String>> = st
        .models
        .iter()
        .enumerate()
        .map(|(mi, m)| {
            m.chains
                .iter()
                .enumerate()
                .map(|(ci, c)| if c.id.trim().is_empty() { lai[mi][ci].first().cloned().unwrap_or_else(|| "A".into()) } else { c.id.clone() })
                .collect()
        })
        .collect();
    let seq_id = |rg: &ResidueGroup| -> String {
        let rs = rg.resseq.trim();
        if rs.len() == 4 { hy36_decode(rs).map(|v| v.to_string()).unwrap_or_else(|| rs.to_string()) } else { rs.to_string() }
    };
    fn icode(rg: &ResidueGroup) -> CifCell<'_> {
        if rg.icode.trim().is_empty() { CifCell::Unknown } else { CifCell::Text(&rg.icode) }
    }
    fn alt(ag: &AtomGroup) -> CifCell<'_> {
        if ag.altloc.is_empty() { CifCell::NotApplicable } else { CifCell::Text(&ag.altloc) }
    }
    fn cell_or(s: &str) -> CifCell<'_> {
        match s {
            "?" => CifCell::Unknown,
            "." => CifCell::NotApplicable,
            _ => CifCell::Text(s),
        }
    }
    // reused buffers for formatted numbers
    let mut f: [String; 6] = Default::default();
    fn set(f: &mut [String; 6], vals: &[(f64, usize)]) {
        for (k, &(v, prec)) in vals.iter().enumerate() {
            f[k].clear();
            let _ = write!(f[k], "{:.*}", prec, v);
        }
    }

    let header = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.type_symbol",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ];
    sink.begin_loop(&header, st.atoms_size());
    let mut n_aniso = 0usize;
    let mut serial_text = String::new();
    let mut seq_text = String::new();
    let mut seq_of = None;
    for (mi, ci, ri, rg, ag, a, serial) in output_atoms(st) {
        if seq_of != Some((mi, ci, ri)) {
            seq_text = seq_id(rg);
            seq_of = Some((mi, ci, ri));
        }
        serial_text.clear();
        let _ = write!(serial_text, "{}", serial);
        let charge = match charge_tidy(&a.charge) {
            Some(s) if s.trim().len() == 2 => {
                let b = s.trim().as_bytes();
                if b[1] == b'-' { format!("-{}", b[0] as char) } else { (b[0] as char).to_string() }
            }
            _ => String::new(),
        };
        set(&mut f, &[(a.xyz.x, 5), (a.xyz.y, 5), (a.xyz.z, 5), (a.occ, 3), (a.b, 5)]);
        let name = a.name.trim();
        sink.row(&[
            CifCell::Text(if a.hetero { "HETATM" } else { "ATOM" }),
            CifCell::Text(&serial_text),
            CifCell::Text(name),
            alt(ag),
            CifCell::Text(ag.resname.trim()),
            CifCell::Text(&auth_asyms[mi][ci]),
            CifCell::Text(&seq_text),
            icode(rg),
            CifCell::Text(&f[0]),
            CifCell::Text(&f[1]),
            CifCell::Text(&f[2]),
            CifCell::Text(&f[3]),
            CifCell::Text(&f[4]),
            CifCell::Text(a.elem()),
            if charge.is_empty() { CifCell::Unknown } else { CifCell::Text(&charge) },
            CifCell::Text(&lai[mi][ci][ri]),
            CifCell::Unknown,
            cell_or(&lsi[mi][ci][ri]),
            CifCell::Text(name),
            CifCell::Text(&model_nums[mi]),
        ]);
        if a.uij.is_some() {
            n_aniso += 1;
        }
    }
    sink.end_loop();
    if n_aniso > 0 {
        sink.begin_loop(
            &[
                "_atom_site_anisotrop.id",
                "_atom_site_anisotrop.pdbx_auth_atom_id",
                "_atom_site_anisotrop.pdbx_label_alt_id",
                "_atom_site_anisotrop.pdbx_auth_comp_id",
                "_atom_site_anisotrop.pdbx_auth_asym_id",
                "_atom_site_anisotrop.pdbx_auth_seq_id",
                "_atom_site_anisotrop.pdbx_PDB_ins_code",
                "_atom_site_anisotrop.U[1][1]",
                "_atom_site_anisotrop.U[2][2]",
                "_atom_site_anisotrop.U[3][3]",
                "_atom_site_anisotrop.U[1][2]",
                "_atom_site_anisotrop.U[1][3]",
                "_atom_site_anisotrop.U[2][3]",
            ],
            n_aniso,
        );
        for (mi, ci, _, rg, ag, a, serial) in output_atoms(st) {
            let Some(u) = a.uij else { continue };
            serial_text.clear();
            let _ = write!(serial_text, "{}", serial);
            set(&mut f, &[(u[0], 5), (u[1], 5), (u[2], 5), (u[3], 5), (u[4], 5), (u[5], 5)]);
            let seq = seq_id(rg);
            let name = a.name.trim();
            sink.row(&[
                CifCell::Text(&serial_text),
                CifCell::Text(name),
                alt(ag),
                CifCell::Text(ag.resname.trim()),
                CifCell::Text(&auth_asyms[mi][ci]),
                CifCell::Text(&seq),
                icode(rg),
                CifCell::Text(&f[0]),
                CifCell::Text(&f[1]),
                CifCell::Text(&f[2]),
                CifCell::Text(&f[3]),
                CifCell::Text(&f[4]),
                CifCell::Text(&f[5]),
            ]);
        }
        sink.end_loop();
    }
}

/// `_atom_site` items that describe the residue rather than the atom: a new
/// hydrogen takes them from an input atom of its residue.
const RESIDUE_ITEMS: &[&str] = &[
    "group_pdb",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_seq_id",
    "pdbx_pdb_ins_code",
    "pdbx_pdb_model_num",
    "pdbx_label_index",
    "pdbx_pdb_residue_no",
    "pdbx_pdb_residue_name",
    "pdbx_pdb_strand_id",
    "pdbx_sifts_xref_db_acc",
    "pdbx_sifts_xref_db_name",
    "pdbx_sifts_xref_db_num",
    "pdbx_sifts_xref_db_res",
];

#[derive(Clone, Copy, PartialEq)]
enum AtomItem {
    Id,
    TypeSymbol,
    Name,
    AltId,
    Coord(usize),
    Occupancy,
    BIso,
    Charge,
    Residue,
    Other,
}

/// A value to send: from the source table, or made here.
enum Val {
    Source(usize, usize),
    /// A byte range of the row's scratch text.
    Made(usize, usize),
    Missing(CifCell<'static>),
}

/// Append `text` to the row's scratch text and name its range.
fn made(scratch: &mut String, text: &str) -> Val {
    let s = scratch.len();
    scratch.push_str(text);
    Val::Made(s, scratch.len())
}

/// `format!("{:.prec$}", v)` into the row's scratch text.
fn made_fixed(scratch: &mut String, v: f64, prec: u32) -> Val {
    let s = scratch.len();
    crate::fastfmt::push_fixed(scratch, v, prec, 0);
    Val::Made(s, scratch.len())
}

fn made_int(scratch: &mut String, v: usize) -> Val {
    let s = scratch.len();
    let mut buf = [0u8; 20];
    let mut k = buf.len();
    let mut x = v;
    loop {
        k -= 1;
        buf[k] = b'0' + (x % 10) as u8;
        x /= 10;
        if x == 0 {
            break;
        }
    }
    scratch.push_str(std::str::from_utf8(&buf[k..]).unwrap());
    Val::Made(s, scratch.len())
}

fn send_row<T: CifTable + ?Sized, K: CifSink + ?Sized>(t: &T, vals: &[Val], scratch: &str, sink: &mut K) {
    let cells: smallvec::SmallVec<[(Option<CifCell<'static>>, std::borrow::Cow<'_, str>); 32]> = vals
        .iter()
        .map(|v| match v {
            Val::Source(r, c) => (t.missing(*r, *c), t.cell(*r, *c)),
            Val::Made(a, b) => (None, std::borrow::Cow::Borrowed(&scratch[*a..*b])),
            Val::Missing(m) => (Some(*m), std::borrow::Cow::Borrowed("")),
        })
        .collect();
    let row: smallvec::SmallVec<[CifCell<'_>; 32]> = cells.iter().map(|(m, s)| m.unwrap_or(CifCell::Text(s))).collect();
    sink.row(&row);
}

/// Send a source category unchanged.
fn send_table<T: CifTable + ?Sized, K: CifSink + ?Sized>(t: &T, sink: &mut K) {
    let tags: Vec<String> = t.tags().iter().map(|s| s.to_string()).collect();
    let ncol = tags.len();
    if t.is_loop() {
        let refs: Vec<&str> = tags.iter().map(|s| s.as_str()).collect();
        sink.begin_loop(&refs, t.row_count());
        let mut vals: Vec<Val> = Vec::with_capacity(ncol);
        for r in 0..t.row_count() {
            vals.clear();
            vals.extend((0..ncol).map(|c| Val::Source(r, c)));
            send_row(t, &vals, "", sink);
        }
        sink.end_loop();
    } else if t.row_count() > 0 {
        for (c, tag) in tags.iter().enumerate() {
            let cell = t.cell(0, c);
            sink.item(tag, t.missing(0, c).unwrap_or(CifCell::Text(&cell)));
        }
    }
}

/// Write the model into its source data block: every category of the source in
/// its order, with `_atom_site` rebuilt and `_atom_site_anisotrop` renumbered.
///
/// Input atoms keep their source values, including the label identifiers that
/// `_struct_conn`, the sequence schemes and other categories refer to; only the
/// atom ids are renumbered (1 to N over all models), and coordinates,
/// occupancies and B values are rewritten where Reduce3 changed them. A new
/// hydrogen takes the residue items ([`RESIDUE_ITEMS`]) from an input atom of
/// its residue and is unknown (`?`) in every other item it has no value for.
/// `_atom_type` gains the elements it lacks. The other categories are sent
/// unchanged (or copied by the sink itself, see [`CifSink::copy_category`]).
pub fn write_cif_preserving<S: CifSource + ?Sized, K: CifSink + ?Sized>(
    st: &Structure,
    source: &S,
    code: &str,
    sink: &mut K,
) -> Result<(), String> {
    let atoms = source.table("atom_site").ok_or_else(|| "the source has no _atom_site loop".to_string())?;
    let tags: Vec<String> = atoms.tags().iter().map(|t| t.to_string()).collect();
    let item_of = |t: &str| t.split_once('.').map(|x| x.1).unwrap_or("").to_ascii_lowercase();
    let roles: Vec<AtomItem> = tags
        .iter()
        .map(|t| match item_of(t).as_str() {
            "id" => AtomItem::Id,
            "type_symbol" => AtomItem::TypeSymbol,
            "label_atom_id" | "auth_atom_id" => AtomItem::Name,
            "label_alt_id" => AtomItem::AltId,
            "cartn_x" => AtomItem::Coord(0),
            "cartn_y" => AtomItem::Coord(1),
            "cartn_z" => AtomItem::Coord(2),
            "occupancy" => AtomItem::Occupancy,
            "b_iso_or_equiv" => AtomItem::BIso,
            "pdbx_formal_charge" => AtomItem::Charge,
            i if RESIDUE_ITEMS.contains(&i) => AtomItem::Residue,
            _ => AtomItem::Other,
        })
        .collect();
    let c_id = roles.iter().position(|&r| r == AtomItem::Id);
    // output atoms in hierarchy order, each with an input row to take values from
    let mut out: Vec<(&AtomGroup, &Atom, Option<usize>)> = Vec::new();
    for m in &st.models {
        for c in &m.chains {
            for rg in &c.residue_groups {
                let rg_row = rg.atom_groups.iter().flat_map(|g| &g.atoms).find(|a| a.src != Atom::NEW).map(|a| a.src as usize);
                for ag in &rg.atom_groups {
                    let ag_row = ag.atoms.iter().find(|a| a.src != Atom::NEW).map(|a| a.src as usize).or(rg_row);
                    for a in &ag.atoms {
                        out.push((ag, a, if a.src != Atom::NEW { Some(a.src as usize) } else { ag_row }));
                    }
                }
            }
        }
    }
    sink.begin_block(code);
    for cat in source.categories() {
        match cat.to_ascii_lowercase().as_str() {
            "atom_site" => {
                let refs: Vec<&str> = tags.iter().map(|s| s.as_str()).collect();
                sink.begin_loop(&refs, out.len());
                let mut scratch = String::new();
                let mut vals: Vec<Val> = Vec::with_capacity(roles.len());
                for (k, &(ag, a, row)) in out.iter().enumerate() {
                    let input = a.src != Atom::NEW;
                    scratch.clear();
                    vals.clear();
                    for (c, &role) in roles.iter().enumerate() {
                        let v = match role {
                            AtomItem::Id => made_int(&mut scratch, k + 1),
                            AtomItem::Coord(i) => {
                                let v = [a.xyz.x, a.xyz.y, a.xyz.z][i];
                                match row {
                                    Some(r) if input && atoms.number(r, c) == Some(v) => Val::Source(r, c),
                                    _ => made_fixed(&mut scratch, v, 3),
                                }
                            }
                            AtomItem::Occupancy | AtomItem::BIso => {
                                let v = if role == AtomItem::Occupancy { a.occ } else { a.b };
                                match row {
                                    Some(r) if input && atoms.number(r, c) == Some(v) => Val::Source(r, c),
                                    _ => made_fixed(&mut scratch, v, 2),
                                }
                            }
                            AtomItem::Name => match row {
                                Some(r) if input && atoms.cell(r, c).trim() == a.name.trim() => Val::Source(r, c),
                                _ => made(&mut scratch, a.name.trim()),
                            },
                            _ if input => Val::Source(row.unwrap_or(0), c),
                            AtomItem::TypeSymbol => made(&mut scratch, a.elem()),
                            AtomItem::AltId if ag.altloc.is_empty() => Val::Missing(CifCell::NotApplicable),
                            AtomItem::AltId => made(&mut scratch, &ag.altloc),
                            AtomItem::Residue => match row {
                                Some(r) => Val::Source(r, c),
                                None => Val::Missing(CifCell::Unknown),
                            },
                            AtomItem::Charge | AtomItem::Other => Val::Missing(CifCell::Unknown),
                        };
                        vals.push(v);
                    }
                    send_row(&atoms, &vals, &scratch, sink);
                }
                sink.end_loop();
            }
            "atom_site_anisotrop" => {
                let Some(an) = source.table(&cat) else { continue };
                let (Some(c_id), Some(an_id)) = (c_id, an.column("id")) else { continue };
                let by_id: rustc_hash::FxHashMap<String, usize> =
                    (0..an.row_count()).map(|r| (an.cell(r, an_id).into_owned(), r)).collect();
                let rows: Vec<(usize, usize)> = out
                    .iter()
                    .enumerate()
                    .filter(|(_, (_, a, _))| a.src != Atom::NEW)
                    .filter_map(|(k, (_, a, _))| by_id.get(&*atoms.cell(a.src as usize, c_id)).map(|&r| (k, r)))
                    .collect();
                if rows.is_empty() {
                    continue;
                }
                let an_tags: Vec<String> = an.tags().iter().map(|s| s.to_string()).collect();
                let refs: Vec<&str> = an_tags.iter().map(|s| s.as_str()).collect();
                sink.begin_loop(&refs, rows.len());
                let mut scratch = String::new();
                for (k, r) in rows {
                    scratch.clear();
                    let vals: Vec<Val> =
                        (0..an_tags.len()).map(|c| if c == an_id { made_int(&mut scratch, k + 1) } else { Val::Source(r, c) }).collect();
                    send_row(&an, &vals, &scratch, sink);
                }
                sink.end_loop();
            }
            "atom_type" => {
                let Some(t) = source.table(&cat) else { continue };
                let present: Vec<String> = match t.column("symbol") {
                    Some(c) => (0..t.row_count()).map(|r| t.cell(r, c).trim().to_ascii_uppercase()).collect(),
                    None => Vec::new(),
                };
                let mut missing: Vec<String> = Vec::new();
                for (_, a, _) in &out {
                    let e = a.elem().to_ascii_uppercase();
                    if !present.contains(&e) && !missing.contains(&e) {
                        missing.push(e);
                    }
                }
                let Some(c_sym) = t.column("symbol").filter(|_| !missing.is_empty()) else {
                    if !sink.copy_category(&cat) {
                        send_table(&t, sink);
                    }
                    continue;
                };
                let t_tags: Vec<String> = t.tags().iter().map(|s| s.to_string()).collect();
                let refs: Vec<&str> = t_tags.iter().map(|s| s.as_str()).collect();
                sink.begin_loop(&refs, t.row_count() + missing.len());
                for r in 0..t.row_count() {
                    let vals: Vec<Val> = (0..t_tags.len()).map(|c| Val::Source(r, c)).collect();
                    send_row(&t, &vals, "", sink);
                }
                let mut scratch = String::new();
                for e in missing {
                    scratch.clear();
                    let vals: Vec<Val> = (0..t_tags.len())
                        .map(|c| if c == c_sym { made(&mut scratch, &e) } else { Val::Missing(CifCell::Unknown) })
                        .collect();
                    send_row(&t, &vals, &scratch, sink);
                }
                sink.end_loop();
            }
            _ => {
                if sink.copy_category(&cat) {
                    continue;
                }
                if let Some(t) = source.table(&cat) {
                    send_table(&t, sink);
                }
            }
        }
    }
    Ok(())
}
