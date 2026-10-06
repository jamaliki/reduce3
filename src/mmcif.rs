//! mmCIF reading and writing with the conventions of iotbx
//! (`iotbx.pdb.mmcif.pdb_hierarchy_builder` and `hierarchy.as_cif_block`).

use crate::cell::{space_group_info, UnitCell};
use crate::cif;
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
    let cat = block.category("_atom_site").unwrap();
    let col = |names: &[&str]| names.iter().find_map(|n| cat.col(n));
    let req = |names: &[&str]| col(names).ok_or_else(|| format!("_atom_site.{} is missing", names[0]));
    let c_type = req(&["type_symbol"])?;
    let c_name = req(&["label_atom_id", "auth_atom_id"])?;
    let c_alt = req(&["label_alt_id"])?;
    let c_auth_asym = col(&["auth_asym_id", "label_asym_id"]).ok_or("_atom_site.auth_asym_id is missing")?;
    let c_label_asym = col(&["label_asym_id", "auth_asym_id"]).unwrap();
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
    let get = |r: usize, c: usize| -> &str { cat.get(r, c) };
    let num = |r: usize, c: usize| -> Result<f64, String> {
        cif::parse_f64(get(r, c)).ok_or_else(|| format!("bad number '{}' in _atom_site row {}", get(r, c), r + 1))
    };
    // anisotropic displacements by atom id
    let mut aniso: rustc_hash::FxHashMap<String, [f64; 6]> = rustc_hash::FxHashMap::default();
    if let Some(an) = block.category("_atom_site_anisotrop") {
        if let Some(ci) = an.col("id") {
            let u: Vec<Option<usize>> =
                ["U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]"].iter().map(|n| an.col(n)).collect();
            let b: Vec<Option<usize>> =
                ["B[1][1]", "B[2][2]", "B[3][3]", "B[1][2]", "B[1][3]", "B[2][3]"].iter().map(|n| an.col(n)).collect();
            let (cols, scale) = if u.iter().all(|x| x.is_some()) {
                (u, 1.0)
            } else if b.iter().all(|x| x.is_some()) {
                (b, 1.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI))
            } else {
                (vec![], 1.0)
            };
            if !cols.is_empty() {
                for r in 0..an.nrows() {
                    let mut v = [0.0; 6];
                    let mut ok = true;
                    for k in 0..6 {
                        match cif::parse_f64(an.get(r, cols[k].unwrap())) {
                            Some(x) => v[k] = x * scale,
                            None => ok = false,
                        }
                    }
                    if ok {
                        aniso.insert(an.get(r, ci).to_string(), v);
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
    for r in 0..cat.nrows() {
        let model_id = c_model.map(|c| get(r, c).to_string()).unwrap_or_else(|| "1".into());
        let new_model = cur_model.as_deref() != Some(model_id.as_str());
        if new_model {
            models.push(Model { id: model_id.clone(), chains: Vec::new() });
        }
        let label_asym = get(r, c_label_asym).to_string();
        let auth_asym = get(r, c_auth_asym).to_string();
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
        let res = get(r, c_seq).to_string();
        let mut ins = c_ins.map(|c| get(r, c).to_string()).unwrap_or_else(|| " ".into());
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
        let mut alt = get(r, c_alt).to_string();
        if alt == "." || alt == "?" {
            alt.clear();
        }
        let comp = get(r, c_comp).to_string();
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
        let mut atom = Atom::new(&format_pdb_atom_name(get(r, c_name), type_symbol), type_symbol, v3(num(r, c_x)?, num(r, c_y)?, num(r, c_z)?));
        atom.occ = num(r, c_occ)?;
        atom.b = num(r, c_b)?;
        if let Some(ci) = c_id {
            let id = get(r, ci);
            atom.serial = id.trim().parse::<i64>().ok().and_then(|v| hy36_encode(5, v)).unwrap_or_else(|| id.to_string());
            if let Some(u) = aniso.get(id) {
                atom.uij = Some(*u);
            }
        }
        atom.segid = match c_segid {
            Some(cs) => {
                let s = get(r, cs);
                let s = if cif::is_null(s) { "" } else { s };
                format!("{:<4}", &s[..s.len().min(4)])
            }
            None => "    ".into(),
        };
        atom.hetero = c_group.map_or(false, |c| get(r, c) == "HETATM");
        if let Some(cc) = c_charge {
            let ch = get(r, cc);
            if !matches!(ch, "?" | ".") {
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
    let cellv = |tag: &str| -> Option<f64> { block.item(tag).and_then(cif::parse_f64) };
    let cell = [
        cellv("_cell.length_a"),
        cellv("_cell.length_b"),
        cellv("_cell.length_c"),
        cellv("_cell.angle_alpha"),
        cellv("_cell.angle_beta"),
        cellv("_cell.angle_gamma"),
    ];
    let mut crystal = None;
    if cell.iter().all(|x| x.is_some()) {
        let sg = ["_symmetry.space_group_name_H-M", "_space_group.name_H-M_alt", "_space_group.name_H-M_full"]
            .iter()
            .find_map(|t| block.item(t).filter(|v| !cif::is_null(v)))
            .unwrap_or("P 1")
            .to_string();
        crystal = Some(CrystalSymmetry { cell: cell.map(|x| x.unwrap()), space_group: sg });
    }
    let had = crystal.is_some();
    let mut st = Structure { models, crystal, had_cell_record: had, links: Vec::new() };
    st.reset_i_seq();
    Ok(st)
}

/// A CIF value, quoted when it would otherwise not read back as one token.
fn cif_value(s: &str) -> String {
    if s.is_empty() {
        return "''".into();
    }
    let special = s.starts_with(['_', '#', '$', '\'', '"', '[', ']', ';'])
        || s.contains([' ', '\t'])
        || matches!(s.to_ascii_lowercase().as_str(), "loop_" | "stop_" | "global_")
        || s.to_ascii_lowercase().starts_with("data_")
        || s.to_ascii_lowercase().starts_with("save_");
    if !special {
        s.to_string()
    } else if !s.contains('\'') {
        format!("'{}'", s)
    } else {
        format!("\"{}\"", s)
    }
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

fn pad_key(out: &mut String, key: &str, value: &str) {
    let _ = writeln!(out, "{:<34}{}", key, value);
}

/// Write the model as mmCIF the way `model_as_mmcif` lays it out.
pub fn write_mmcif(st: &Structure) -> String {
    let mut out = String::with_capacity(st.atoms_size() * 100 + 4096);
    out.push_str("data_default\n");
    let mut symops: &[&str] = &[];
    if let Some(cs) = &st.crystal {
        let sg = space_group_info(&cs.space_group, &cs.cell);
        let uc = crate::cell::cryst1_rotations(&cs.space_group, &cs.cell)
            .and_then(|r| UnitCell::new(cs.cell).and_then(|u| u.averaged(r)))
            .or_else(|| UnitCell::new(cs.cell));
        if let Some(uc) = &uc {
            let p = uc.params;
            for (k, v) in ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"].iter().zip(p.iter()) {
                pad_key(&mut out, &format!("_cell.{}", k), &format!("{:.3}", v));
            }
            pad_key(&mut out, "_cell.volume", &format!("{:.3}", uc.volume));
        }
        if let Some(sg) = &sg {
            pad_key(&mut out, "_space_group.crystal_system", sg.crystal_system);
            pad_key(&mut out, "_space_group.IT_number", &sg.number.to_string());
            pad_key(&mut out, "_space_group.name_H-M_alt", &cif_value(sg.hm));
            pad_key(&mut out, "_space_group.name_Hall", &cif_value(sg.hall));
            pad_key(&mut out, "_symmetry.space_group_name_H-M", &cif_value(sg.hm));
            pad_key(&mut out, "_symmetry.space_group_name_Hall", &cif_value(sg.hall));
            pad_key(&mut out, "_symmetry.Int_Tables_number", &sg.number.to_string());
            symops = sg.ops;
        }
    }
    if !symops.is_empty() {
        out.push_str("loop_\n  _space_group_symop.id\n  _space_group_symop.operation_xyz\n");
        for (k, op) in symops.iter().enumerate() {
            let _ = writeln!(out, "   {} {}", k + 1, op);
        }
        out.push('\n');
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
    out.push_str("loop_\n  _struct_asym.id\n");
    for s in &struct_asym {
        let _ = writeln!(out, "   {}", cif_value(s));
    }
    out.push_str("\nloop_\n  _chem_comp.id\n");
    for c in &comps {
        let _ = writeln!(out, "   {}", cif_value(c));
    }
    out.push('\n');
    let header = [
        "group_PDB", "id", "label_atom_id", "label_alt_id", "label_comp_id", "auth_asym_id", "auth_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol",
        "pdbx_formal_charge", "label_asym_id", "label_entity_id", "label_seq_id", "auth_atom_id", "pdbx_PDB_model_num",
    ];
    out.push_str("loop_\n");
    for h in header {
        let _ = writeln!(out, "  _atom_site.{}", h);
    }
    let mut aniso = String::new();
    for (mi, m) in st.models.iter().enumerate() {
        let model_num = if m.id.trim().is_empty() { "1".to_string() } else { m.id.trim().to_string() };
        let mut serial = 1i64;
        for (ci, c) in m.chains.iter().enumerate() {
            let auth_asym = if c.id.trim().is_empty() { lai[mi][ci].first().cloned().unwrap_or_else(|| "A".into()) } else { c.id.clone() };
            for (ri, rg) in c.residue_groups.iter().enumerate() {
                let rs = rg.resseq.trim();
                let seq_id = if rs.len() == 4 { hy36_decode(rs).map(|v| v.to_string()).unwrap_or_else(|| rs.to_string()) } else { rs.to_string() };
                let icode = if rg.icode.trim().is_empty() { "?".to_string() } else { rg.icode.clone() };
                for ag in &rg.atom_groups {
                    let comp = ag.resname.trim();
                    let alt = if ag.altloc.is_empty() { "." } else { ag.altloc.as_str() };
                    for a in &ag.atoms {
                        let charge = match charge_tidy(&a.charge) {
                            Some(s) if s.trim().len() == 2 => {
                                let b = s.trim().as_bytes();
                                if b[1] == b'-' { format!("-{}", b[0] as char) } else { (b[0] as char).to_string() }
                            }
                            _ => "?".to_string(),
                        };
                        let name = a.name.trim();
                        let _ = writeln!(
                            out,
                            "   {} {} {} {} {} {} {} {} {:.5} {:.5} {:.5} {:.3} {:.5} {} {} {} ? {} {} {}",
                            if a.hetero { "HETATM" } else { "ATOM" },
                            serial,
                            cif_value(name),
                            cif_value(alt),
                            cif_value(comp),
                            cif_value(&auth_asym),
                            seq_id,
                            cif_value(&icode),
                            a.xyz.x,
                            a.xyz.y,
                            a.xyz.z,
                            a.occ,
                            a.b,
                            a.elem(),
                            charge,
                            lai[mi][ci][ri],
                            lsi[mi][ci][ri],
                            cif_value(name),
                            model_num
                        );
                        if let Some(u) = a.uij {
                            let _ = writeln!(
                                aniso,
                                "   {} {} {} {} {} {} {} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5}",
                                serial,
                                cif_value(name),
                                cif_value(alt),
                                cif_value(comp),
                                cif_value(&auth_asym),
                                seq_id,
                                cif_value(&icode),
                                u[0],
                                u[1],
                                u[2],
                                u[3],
                                u[4],
                                u[5]
                            );
                        }
                        serial += 1;
                    }
                }
            }
        }
    }
    out.push('\n');
    if !aniso.is_empty() {
        out.push_str("loop_\n");
        for h in [
            "id", "pdbx_auth_atom_id", "pdbx_label_alt_id", "pdbx_auth_comp_id", "pdbx_auth_asym_id", "pdbx_auth_seq_id",
            "pdbx_PDB_ins_code", "U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]",
        ] {
            let _ = writeln!(out, "  _atom_site_anisotrop.{}", h);
        }
        out.push_str(&aniso);
        out.push('\n');
    }
    out
}
