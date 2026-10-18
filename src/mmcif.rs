//! mmCIF reading and writing with the conventions of iotbx
//! (`iotbx.pdb.mmcif.pdb_hierarchy_builder` and `hierarchy.as_cif_block`).

use crate::cell::{space_group_info, UnitCell};
use crate::cif;
use crate::cifsource::{CifCell, CifSink, CifSource, CifTable, CifText, ColumnValue};
use crate::geom::v3;
use crate::model::*;
use crate::resclass::{self, ResClass};
use std::borrow::Cow;
use std::fmt::Write as _;

/// `format_pdb_atom_name`: the 4-column PDB spelling of an mmCIF atom name.
fn format_pdb_atom_name(name: &str, element: &str) -> ShortStr {
    if name.len() >= 4 {
        return name.into();
    }
    let b = name.as_bytes();
    // the usual name: ASCII without blanks, which `trim` leaves as it is
    if b.iter().all(u8::is_ascii_graphic) {
        let lead = element.len() == 1
            && b.first().map(u8::to_ascii_uppercase) == element.as_bytes().first().map(u8::to_ascii_uppercase);
        let mut field = [b' '; 4];
        let at = usize::from(lead);
        field[at..at + b.len()].copy_from_slice(b);
        return ShortStr::from_ascii_field(field);
    }
    let t = name.trim();
    let first = t.chars().next().map(|c| c.to_ascii_uppercase());
    let mut n = ShortStr::new();
    if element.len() == 1 && first == element.chars().next().map(|c| c.to_ascii_uppercase()) {
        n.push(' ');
    }
    n.push_str(t);
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
    read_structure(block, None)
}

/// [`structure_from_cif`], with the numbers it read for
/// [`write_cif_preserving_with`].
pub fn structure_and_numbers_from_cif<S: CifSource + ?Sized>(block: &S) -> Result<(Structure, SourceNumbers), String> {
    let mut numbers = SourceNumbers(Default::default());
    let st = read_structure(block, Some(&mut numbers))?;
    Ok((st, numbers))
}

fn read_structure<S: CifSource + ?Sized>(block: &S, mut numbers: Option<&mut SourceNumbers>) -> Result<Structure, String> {
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

    let nrows = cat.row_count();
    // The first pass reads the hierarchy items and gives each row an atom
    // group; the second reads the atoms into groups sized by the first.
    let mut models: Vec<Model> = Vec::new();
    // the atom groups of each residue group in their final order, by model,
    // chain and residue group
    let mut res_groups: Vec<Vec<Vec<Vec<u32>>>> = Vec::new();
    // altloc, resname and number of rows of each atom group
    let mut groups: Vec<(String, String, usize)> = Vec::new();
    let mut group_of_row: Vec<u32> = Vec::with_capacity(nrows);
    // Rows may return to a chain or residue earlier in the same model block.
    let mut chain_index = rustc_hash::FxHashMap::default();
    let mut residue_index: rustc_hash::FxHashMap<(usize, usize, ShortStr, ShortStr), usize> = rustc_hash::FxHashMap::default();
    // the hierarchy items of the previous row and its atom group: a row with
    // the same items belongs to the same group
    let mut last: Option<[Cow<'_, str>; 7]> = None;
    let mut last_group = 0u32;
    let mut cur_chain = 0;
    // (model, chain, residue number, insertion code) of the last new residue
    // label and its residue group
    let mut cur_res: Option<(usize, usize, ShortStr, ShortStr, usize)> = None;
    let mut resname = String::new();
    // the first row the first pass rejects, and why
    let mut failure: Option<(usize, String)> = None;
    for r in 0..nrows {
        let items = [
            c_model.map(|c| get(r, c)).unwrap_or(Cow::Borrowed("1")),
            get(r, c_label_asym),
            get(r, c_auth_asym),
            get(r, c_seq),
            c_ins.map(|c| get(r, c)).unwrap_or(Cow::Borrowed("\0")),
            get(r, c_alt),
            get(r, c_comp),
        ];
        if last.as_ref().is_some_and(|l| l.iter().zip(&items).all(|(a, b)| same_text(a, b))) {
            groups[last_group as usize].2 += 1;
            group_of_row.push(last_group);
            continue;
        }
        let [model_id, label_asym, auth_asym, res, ins_cell, alt_cell, comp] = &items;
        let prev = last.as_ref();
        let new_model = prev.map(|p| &*p[0]) != Some(&**model_id);
        if new_model {
            models.push(Model { id: model_id.to_string(), chains: Vec::new() });
            res_groups.push(Vec::new());
        }
        let mi = models.len() - 1;
        if matches!(&**auth_asym, "." | "?" | " " | "") {
            failure = Some((r, "mmCIF file contains a record with an empty auth_asym_id".into()));
            break;
        }
        let label_changed = prev.map(|p| &*p[1]) != Some(&**label_asym);
        let auth_changed = prev.map(|p| &*p[2]) != Some(&**auth_asym);
        let prev_comp = prev.map(|p| &*p[6]);
        let ci = if !(new_model || label_changed || auth_changed) { cur_chain } else {
            let (model, chain_groups) = (&mut models[mi], &mut res_groups[mi]);
            *chain_index.entry((mi, auth_asym.to_string(), label_asym.to_string())).or_insert_with(|| {
                if auth_changed || new_model || (label_changed && r > 0 && prev_comp.map_or(false, is_aa_or_rna_dna)) {
                    model.chains.push(Chain { id: auth_asym.to_string(), residue_groups: Vec::new() });
                    chain_groups.push(Vec::new());
                    model.chains.len() - 1
                } else { cur_chain }
            })
        };
        let ins = match (c_ins, &**ins_cell) {
            (None, _) | (Some(_), "?" | ".") => " ",
            (Some(_), s) => s,
        };
        let ri = match &cur_res {
            Some((m, c, rs, is, ri)) if *m == mi && *c == ci && *rs == **res && *is == *ins => *ri,
            _ => {
                let key = (mi, ci, ShortStr::from(&**res), ShortStr::from(ins));
                let ri = if let Some(&index) = residue_index.get(&key) { index } else {
                    let resseq = match res.trim().parse::<i64>() {
                        Ok(v) => hy36_encode(4, v).unwrap_or_else(|| res.to_string()),
                        Err(_) => {
                            if res.len() != 4 {
                                failure = Some((r, format!("bad residue number '{}' in _atom_site row {}", res, r + 1)));
                                break;
                            }
                            res.to_string()
                        }
                    };
                    let chain = &mut models[mi].chains[ci];
                    let mut link = true;
                    if let Some(cb) = c_break {
                        if !chain.residue_groups.is_empty() && get(r, cb) == "1" {
                            link = false;
                        }
                    }
                    chain.residue_groups.push(ResidueGroup {
                        resseq,
                        icode: ins.to_string(),
                        link_to_previous: link,
                        atom_groups: Vec::new(),
                    });
                    res_groups[mi][ci].push(Vec::new());
                    let index = chain.residue_groups.len() - 1;
                    residue_index.insert(key.clone(), index);
                    index
                };
                cur_res = Some((key.0, key.1, key.2, key.3, ri));
                ri
            }
        };
        let alt = if matches!(&**alt_cell, "." | "?") { "" } else { &**alt_cell };
        // `format!("{:>3}", comp)`, which pads by characters
        resname.clear();
        resname.extend(std::iter::repeat_n(' ', 3usize.saturating_sub(comp.chars().count())));
        resname.push_str(comp);
        let order = &mut res_groups[mi][ci][ri];
        let g = match order.iter().find(|&&g| groups[g as usize].0 == alt && groups[g as usize].1 == resname) {
            Some(&g) => g,
            None => {
                let g = groups.len() as u32;
                groups.push((alt.to_string(), resname.clone(), 0));
                if alt.is_empty() {
                    order.insert(0, g);
                } else {
                    order.push(g);
                }
                g
            }
        };
        groups[g as usize].2 += 1;
        group_of_row.push(g);
        last = Some(items);
        last_group = g;
        cur_chain = ci;
    }
    // each group's first atom in structure order, which sets `i_seq`
    let mut first_atom = vec![0usize; groups.len()];
    let mut next = 0;
    for &g in res_groups.iter().flatten().flatten().flatten() {
        first_atom[g as usize] = next;
        next += groups[g as usize].2;
    }
    let mut atoms: Vec<Vec<Atom>> = groups.iter().map(|g| Vec::with_capacity(g.2)).collect();
    if let Some(SourceNumbers(items)) = numbers.as_deref_mut() {
        items.iter_mut().for_each(|values| values.reserve(nrows));
    }
    // Atom::new formats the element symbol, which depends on nothing else;
    // keep one per spelling, with one-byte spellings also found by their byte
    let mut elements: Vec<(String, ShortStr)> = Vec::new();
    let mut element_of_byte: Vec<Option<ShortStr>> = vec![None; 128];
    let (blank_segid, blank_charge) = (ShortStr::from("    "), ShortStr::from("  "));
    // the last segid and charge cells, with the fields they give
    let mut last_segid: Option<(Cow<'_, str>, ShortStr)> = None;
    let mut last_charge: Option<(Cow<'_, str>, ShortStr)> = None;
    let rows = failure.as_ref().map_or(nrows, |f| f.0);
    for (r, &g) in group_of_row.iter().enumerate().take(rows) {
        let type_symbol = get(r, c_type);
        let byte = match type_symbol.as_bytes() {
            &[c] => element_of_byte.get_mut(c as usize),
            _ => None,
        };
        let element = match byte {
            Some(Some(e)) => e.clone(),
            byte => {
                let e = match elements.iter().find(|(t, _)| same_text(t, &type_symbol)) {
                    Some((_, e)) => e.clone(),
                    None => {
                        let e = Atom::new("", &type_symbol, v3(0.0, 0.0, 0.0)).element;
                        elements.push((type_symbol.to_string(), e.clone()));
                        e
                    }
                };
                if let Some(slot) = byte {
                    *slot = Some(e.clone());
                }
                e
            }
        };
        let name = format_pdb_atom_name(&get(r, c_name), &type_symbol);
        let xyz = v3(num(r, c_x)?, num(r, c_y)?, num(r, c_z)?);
        let occ = num(r, c_occ)?;
        let b = num(r, c_b)?;
        if let Some(SourceNumbers(items)) = numbers.as_deref_mut() {
            for (values, v) in items.iter_mut().zip([xyz.x, xyz.y, xyz.z, occ, b]) {
                // rows are read in order, so row r goes to index r
                values.resize(r, f64::NAN);
                values.push(v);
            }
        }
        let (mut serial, mut uij) = (ShortStr::new(), None);
        if let Some(ci) = c_id {
            let id = get(r, ci);
            serial = serial_number(&id);
            if !aniso.is_empty() {
                uij = aniso.get(&*id).copied();
            }
        }
        let segid = match c_segid {
            Some(cs) => {
                let s = get(r, cs);
                match &last_segid {
                    Some((t, seg)) if same_text(t, &s) => seg.clone(),
                    _ => {
                        let seg = segid_field(&s);
                        last_segid = Some((s, seg.clone()));
                        seg
                    }
                }
            }
            None => blank_segid.clone(),
        };
        let hetero = c_group.map_or(false, |c| get(r, c) == "HETATM");
        let charge = match c_charge {
            Some(cc) => {
                let ch = get(r, cc);
                match &last_charge {
                    Some((t, c)) if same_text(t, &ch) => c.clone(),
                    _ => {
                        let c = charge_field(&ch);
                        last_charge = Some((ch, c.clone()));
                        c
                    }
                }
            }
            None => blank_charge.clone(),
        };
        let group = &mut atoms[g as usize];
        let i_seq = first_atom[g as usize] + group.len();
        group.push(Atom { name, element, charge, serial, xyz, occ, b, segid, hetero, uij, i_seq, src: r as u32 });
    }
    if let Some((_, message)) = failure {
        return Err(message);
    }
    for (model, model_groups) in models.iter_mut().zip(&res_groups) {
        for (chain, chain_groups) in model.chains.iter_mut().zip(model_groups) {
            for (rg, order) in chain.residue_groups.iter_mut().zip(chain_groups) {
                rg.atom_groups = order.iter().map(|&g| {
                    let (altloc, resname, _) = std::mem::take(&mut groups[g as usize]);
                    AtomGroup { altloc, resname, atoms: std::mem::take(&mut atoms[g as usize]) }
                }).collect();
            }
        }
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
    let st = Structure { models, crystal, had_cell_record: had, links: covalent_links(block), pdb_records: PdbRecords::default() };
    Ok(st)
}

/// `a == b` by a loop over the bytes: the values compared are a few bytes
/// long, and a call to `memcmp` takes longer than the loop.
#[inline]
fn same_text(a: &str, b: &str) -> bool {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() != b.len() {
        return false;
    }
    let mut k = 0;
    while k < a.len() {
        if a[k] != b[k] {
            return false;
        }
        k += 1;
    }
    true
}

/// The serial number field of `_atom_site.id`: `hy36_encode(5, id)`, or the
/// text itself when it is not a number.
fn serial_number(id: &str) -> ShortStr {
    let b = id.as_bytes();
    // a plain decimal of at most five digits is the number right-justified
    if (1..=5).contains(&b.len()) && b.iter().all(u8::is_ascii_digit) && (b.len() == 1 || b[0] != b'0') {
        let mut field = [b' '; 5];
        field[5 - b.len()..].copy_from_slice(b);
        return ShortStr::from_ascii_field(field);
    }
    id.trim().parse::<i64>().ok().and_then(|v| hy36_short(5, v)).unwrap_or_else(|| id.into())
}

/// The four-column segment field of `_atom_site.auth_segid`.
fn segid_field(s: &str) -> ShortStr {
    let s = if cif::is_null(s) { "" } else { s };
    // `format!("{:<4}", ..)`, which pads by characters
    let s = &s[..s.len().min(4)];
    let mut seg = ShortStr::from(s);
    for _ in s.chars().count()..4 {
        seg.push(' ');
    }
    seg
}

/// The two-column charge field of `_atom_site.pdbx_formal_charge`.
fn charge_field(ch: &str) -> ShortStr {
    if !matches!(ch, "?" | ".") {
        let sign = if ch.ends_with('-') || ch.starts_with('-') { "-" } else { "+" };
        if let Ok(v) = ch.trim_matches(|c| c == ' ' || c == '-' || c == '+').parse::<i32>() {
            return if v == 0 {
                "  ".into()
            } else {
                let mut c = ShortStr::from(v.to_string());
                c.push_str(sign);
                c
            };
        }
    }
    "  ".into()
}

fn covalent_links<S: CifSource + ?Sized>(source: &S) -> Vec<LinkRecord> {
    let Some(table) = source.table("struct_conn") else { return Vec::new() };
    let value = |row, name: &str| {
        table.column(name).map(|col| table.cell(row, col).into_owned())
            .filter(|s| !cif::is_null(s)).unwrap_or_default()
    };
    let partner = |row, side| {
        let field = |stem: &str| value(row, &format!("ptnr{side}_{stem}"));
        let fallback = |auth, label| {
            let result = field(auth);
            if result.is_empty() { field(label) } else { result }
        };
        AtomLabel {
            name: field("label_atom_id"),
            altloc: value(row, &format!("pdbx_ptnr{side}_label_alt_id")),
            resname: fallback("auth_comp_id", "label_comp_id"),
            chain: fallback("auth_asym_id", "label_asym_id"),
            resseq: fallback("auth_seq_id", "label_seq_id"),
            icode: value(row, &format!("pdbx_ptnr{side}_PDB_ins_code")),
        }
    };
    (0..table.row_count()).filter_map(|row| {
        let kind = value(row, "conn_type_id");
        if !matches!(kind.as_str(), "covale" | "disulf")
            || value(row, "ptnr1_symmetry") != "1_555"
            || value(row, "ptnr2_symmetry") != "1_555" {
            return None;
        }
        let (atom1, atom2) = (partner(row, 1), partner(row, 2));
        if [&atom1, &atom2].iter().any(|a| a.name.is_empty() || a.resname.is_empty() || a.chain.is_empty() || a.resseq.is_empty()) {
            return None;
        }
        Some(LinkRecord {
            atom1, atom2, distance: table.column("pdbx_dist_value").and_then(|col| table.number(row, col)), kind,
            order: Some(value(row, "pdbx_value_order")).filter(|v| !v.is_empty()),
        })
    }).collect()
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
#[derive(Clone, Copy)]
enum Val {
    Source(usize, usize),
    /// A byte range of the row's scratch text.
    Made(usize, usize),
    Unknown,
    NotApplicable,
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
    push_fixed_short(scratch, v, prec);
    Val::Made(s, scratch.len())
}

/// `format!("{:.prec$}", x)`, as [`crate::fastfmt::push_fixed`] writes it,
/// in 64-bit arithmetic when `prec` is at most 3 and `|x|` is below 2^52.
fn push_fixed_short(out: &mut String, x: f64, prec: u32) {
    const SCALE: [u64; 4] = [1, 10, 100, 1000];
    let bits = x.to_bits();
    let exp = ((bits >> 52) & 0x7ff) as u32;
    // below exponent 1075 the binary exponent e is negative, so m * 10^prec
    // (under 2^63) needs no left shift; this also excludes NaN and infinity
    let Some(&scale) = SCALE.get(prec as usize).filter(|_| exp < 1075) else {
        crate::fastfmt::push_fixed(out, x, prec, 0);
        return;
    };
    let frac = bits & ((1u64 << 52) - 1);
    let (m, sh) = if exp == 0 { (frac, 1074) } else { (frac | (1u64 << 52), 1075 - exp) };
    let p = m * scale;
    // n = round_half_even(p / 2^sh); for sh >= 64, p / 2^sh is below 1/2
    let n = if sh >= 64 {
        0
    } else {
        let q = p >> sh;
        let r = p & ((1u64 << sh) - 1);
        let half = 1u64 << (sh - 1);
        if r > half || (r == half && q & 1 == 1) { q + 1 } else { q }
    };
    let (mut int, mut fr) = (n / scale, n % scale);
    let mut buf = [0u8; 32];
    let mut k = buf.len();
    for _ in 0..prec {
        k -= 1;
        buf[k] = b'0' + (fr % 10) as u8;
        fr /= 10;
    }
    if prec > 0 {
        k -= 1;
        buf[k] = b'.';
    }
    loop {
        k -= 1;
        buf[k] = b'0' + (int % 10) as u8;
        int /= 10;
        if int == 0 {
            break;
        }
    }
    if bits >> 63 == 1 {
        k -= 1;
        buf[k] = b'-';
    }
    // the bytes are ASCII, so the check passes
    if let Ok(text) = std::str::from_utf8(&buf[k..]) {
        out.push_str(text);
    }
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
            Val::Source(r, c) => t.cell_with_missing(*r, *c),
            Val::Made(a, b) => (None, std::borrow::Cow::Borrowed(&scratch[*a..*b])),
            Val::Unknown => (Some(CifCell::Unknown), std::borrow::Cow::Borrowed("")),
            Val::NotApplicable => (Some(CifCell::NotApplicable), std::borrow::Cow::Borrowed("")),
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

/// The coordinates, occupancy and B value of each row of a source
/// `_atom_site` loop as [`structure_and_numbers_from_cif`] read them.
/// [`write_cif_preserving_with`] compares the output with these numbers
/// instead of parsing the source text again.
pub struct SourceNumbers([Vec<f64>; 5]);

/// The items of [`SourceNumbers`], as `structure_from_cif` finds their columns.
const NUMBER_ITEMS: [&str; 5] = ["Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"];

/// `s.trim()`, with a byte loop for ASCII ends.
fn trim(s: &str) -> &str {
    // the ASCII characters `char::is_whitespace` accepts
    let blank = |b: &u8| matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
    let b = s.as_bytes();
    let start = b.iter().position(|c| !blank(c)).unwrap_or(b.len());
    let end = b.iter().rposition(|c| !blank(c)).map_or(start, |i| i + 1);
    match s.get(start..end) {
        Some(t) if t.as_bytes().first().is_none_or(u8::is_ascii) && t.as_bytes().last().is_none_or(u8::is_ascii) => t,
        _ => s.trim(),
    }
}

/// Chooses each value of the output `_atom_site` loop.
struct AtomSiteValues<'s, T: ?Sized> {
    atoms: &'s T,
    roles: Vec<AtomItem>,
    /// Output atoms in hierarchy order.
    out: Vec<(&'s AtomGroup, &'s Atom)>,
    /// For each output atom, the input row to take values from (`u32::MAX`
    /// for none) and whether it is an input atom.
    rows: Vec<(u32, bool)>,
    /// x, y, z, occupancy and B of each output atom.
    numbers: [Vec<f64>; 5],
    /// The name of each output atom, trimmed.
    names: Vec<&'s str>,
    source_numbers: Option<&'s SourceNumbers>,
    /// The column `structure_from_cif` read for each item of `source_numbers`.
    number_columns: [Option<usize>; 5],
}

impl<T: CifTable + ?Sized> AtomSiteValues<'_, T> {
    /// Whether `cell`, trimmed, equals the trimmed name `name`.
    fn same_name(cell: &str, name: &str) -> bool {
        cell == name || trim(cell) == name
    }

    /// Whether the number in row `r`, column `c` of the source equals `v`;
    /// `item` indexes [`NUMBER_ITEMS`].
    fn unchanged(&self, r: usize, c: usize, item: usize, v: f64) -> bool {
        let read = self.source_numbers.filter(|_| self.number_columns[item] == Some(c));
        match read.and_then(|n| n.0[item].get(r)) {
            Some(&x) if !x.is_nan() => x == v,
            _ => self.atoms.number(r, c) == Some(v),
        }
    }

    /// The value of output atom `k` in column `c`; made text goes to `scratch`.
    fn value(&self, k: usize, c: usize, scratch: &mut String) -> Val {
        self.value_as(self.roles[c], k, c, scratch)
    }

    /// The values of column `c` in output order; made text goes to `scratch`.
    fn column(&self, c: usize, scratch: &mut String, column: &mut Vec<ColumnValue>) {
        let n = self.rows.len();
        // one loop for each role, in which `value_as` reduces to that role's case
        macro_rules! fill {
            ($role:expr) => {
                column.extend((0..n).map(|k| match self.value_as($role, k, c, scratch) {
                    Val::Source(r, _) => ColumnValue::Source(r),
                    Val::Made(a, b) => ColumnValue::Made(a, b),
                    Val::Unknown => ColumnValue::Unknown,
                    Val::NotApplicable => ColumnValue::NotApplicable,
                }))
            };
        }
        match self.roles[c] {
            AtomItem::Id => fill!(AtomItem::Id),
            AtomItem::TypeSymbol => fill!(AtomItem::TypeSymbol),
            AtomItem::Name => fill!(AtomItem::Name),
            AtomItem::AltId => fill!(AtomItem::AltId),
            AtomItem::Coord(i) => fill!(AtomItem::Coord(i)),
            AtomItem::Occupancy => fill!(AtomItem::Occupancy),
            AtomItem::BIso => fill!(AtomItem::BIso),
            AtomItem::Charge => fill!(AtomItem::Charge),
            AtomItem::Residue => fill!(AtomItem::Residue),
            AtomItem::Other => fill!(AtomItem::Other),
        }
    }

    /// [`AtomSiteValues::value`] for a column with role `role`.
    #[inline(always)]
    fn value_as(&self, role: AtomItem, k: usize, c: usize, scratch: &mut String) -> Val {
        let (row, input) = self.rows[k];
        let row = (row != u32::MAX).then_some(row as usize);
        match role {
            AtomItem::Id => made_int(scratch, k + 1),
            AtomItem::Coord(i) => {
                let v = self.numbers[i][k];
                match row {
                    Some(r) if input && self.unchanged(r, c, i, v) => Val::Source(r, c),
                    _ => made_fixed(scratch, v, 3),
                }
            }
            AtomItem::Occupancy | AtomItem::BIso => {
                let item = if role == AtomItem::Occupancy { 3 } else { 4 };
                let v = self.numbers[item][k];
                match row {
                    Some(r) if input && self.unchanged(r, c, item, v) => Val::Source(r, c),
                    _ => made_fixed(scratch, v, 2),
                }
            }
            AtomItem::Name => match row {
                Some(r) if input && Self::same_name(&self.atoms.cell(r, c), self.names[k]) => Val::Source(r, c),
                _ => made(scratch, self.names[k]),
            },
            _ if input => Val::Source(row.unwrap_or(0), c),
            AtomItem::TypeSymbol => made(scratch, self.out[k].1.elem()),
            AtomItem::AltId if self.out[k].0.altloc.is_empty() => Val::NotApplicable,
            AtomItem::AltId => made(scratch, &self.out[k].0.altloc),
            AtomItem::Residue => match row {
                Some(r) => Val::Source(r, c),
                None => Val::Unknown,
            },
            AtomItem::Charge | AtomItem::Other => Val::Unknown,
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
    write_cif_preserving_with(st, source, None, code, sink)
}

/// [`write_cif_preserving`], with the numbers
/// [`structure_and_numbers_from_cif`] read from `source`. The output is the
/// same; the writer only skips parsing the source coordinates again.
pub fn write_cif_preserving_with<S: CifSource + ?Sized, K: CifSink + ?Sized>(
    st: &Structure,
    source: &S,
    numbers_read: Option<&SourceNumbers>,
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
    // one pass over the atoms collects everything the columns need
    let atom_count = st.atoms_size();
    let mut out: Vec<(&AtomGroup, &Atom)> = Vec::with_capacity(atom_count);
    let mut rows: Vec<(u32, bool)> = Vec::with_capacity(atom_count);
    let mut numbers: [Vec<f64>; 5] = std::array::from_fn(|_| Vec::with_capacity(atom_count));
    let mut names: Vec<&str> = Vec::with_capacity(atom_count);
    // the element spellings, in order of first use, with their keys
    let mut elements: Vec<&str> = Vec::new();
    let mut element_keys: Vec<u64> = Vec::new();
    // the length and bytes of a spelling of up to 7 bytes, else 0
    let key = |s: &str| match s.len() {
        n @ 1..=7 => s.bytes().enumerate().fold(n as u64, |k, (i, b)| k | u64::from(b) << (8 * (i + 1))),
        _ => 0,
    };
    for m in &st.models {
        for c in &m.chains {
            for rg in &c.residue_groups {
                let rg_row = rg.atom_groups.iter().flat_map(|g| &g.atoms).find(|a| a.src != Atom::NEW).map(|a| a.src);
                for ag in &rg.atom_groups {
                    let ag_row = ag.atoms.iter().find(|a| a.src != Atom::NEW).map(|a| a.src).or(rg_row);
                    for a in &ag.atoms {
                        out.push((ag, a));
                        let input = a.src != Atom::NEW;
                        rows.push((if input { a.src } else { ag_row.unwrap_or(Atom::NEW) }, input));
                        for (values, v) in numbers.iter_mut().zip([a.xyz.x, a.xyz.y, a.xyz.z, a.occ, a.b]) {
                            values.push(v);
                        }
                        names.push(trim(&a.name));
                        let k = key(&a.element);
                        let seen = if k == 0 { elements.contains(&a.element.as_str()) } else { element_keys.contains(&k) };
                        if !seen {
                            elements.push(a.element.as_str());
                            element_keys.push(k);
                        }
                    }
                }
            }
        }
    }
    let values = AtomSiteValues {
        atoms: &atoms,
        roles,
        names,
        out,
        rows,
        numbers,
        source_numbers: numbers_read,
        number_columns: NUMBER_ITEMS.map(|item| atoms.column(item)),
    };
    let out = &values.out;
    sink.begin_block(code);
    for cat in source.categories() {
        match cat.to_ascii_lowercase().as_str() {
            "atom_site" => {
                let refs: Vec<&str> = tags.iter().map(|s| s.as_str()).collect();
                let mut scratch = String::new();
                if sink.begin_source_loop("atom_site", &refs, out.len()) {
                    let mut column: Vec<ColumnValue> = Vec::with_capacity(out.len());
                    for c in 0..refs.len() {
                        scratch.clear();
                        column.clear();
                        values.column(c, &mut scratch, &mut column);
                        sink.source_column(c, &column, &scratch);
                    }
                } else {
                    sink.begin_loop(&refs, out.len());
                    let mut vals: Vec<Val> = Vec::with_capacity(refs.len());
                    for k in 0..out.len() {
                        scratch.clear();
                        vals.clear();
                        vals.extend((0..refs.len()).map(|c| values.value(k, c, &mut scratch)));
                        send_row(&atoms, &vals, &scratch, sink);
                    }
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
                    .filter(|(_, (_, a))| a.src != Atom::NEW)
                    .filter_map(|(k, (_, a))| by_id.get(&*atoms.cell(a.src as usize, c_id)).map(|&r| (k, r)))
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
                for element in &elements {
                    // `Atom::elem`
                    let e = element.trim().to_ascii_uppercase();
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
                        .map(|c| if c == c_sym { made(&mut scratch, &e) } else { Val::Unknown })
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

#[cfg(test)]
mod fixed_tests {
    use super::{push_fixed_short, trim};

    #[test]
    fn trim_matches_str_trim() {
        let chars = [' ', '\t', '\n', '\r', '\u{b}', '\u{c}', '\u{85}', '\u{a0}', '\u{2003}', 'C', 'a', '1', '\u{e9}', '\''];
        let mut s: u64 = 0x9e3779b97f4a7c15;
        for _ in 0..100_000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let text: String = (0..(s % 7)).map(|i| chars[((s >> (8 + 4 * i)) % chars.len() as u64) as usize]).collect();
            assert_eq!(trim(&text), text.trim(), "{:?}", text);
        }
    }

    #[test]
    fn short_fixed_matches_format() {
        let mut s: u64 = 0x2545f4914f6cdd1d;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut values = vec![0.0, -0.0, 0.5, -0.5, 1.5, 2.5, 0.0005, 0.0015, -0.0025, 0.125, 4503599627370495.5, 4503599627370496.0, f64::MIN_POSITIVE, 5e-324, 1e300, f64::NAN, f64::INFINITY, -f64::INFINITY];
        for _ in 0..50_000 {
            let r = next();
            values.push(f64::from_bits(r));
            // coordinates and B values as a file gives them
            values.push((r % 4_000_000) as f64 / 1000.0 - 2000.0);
            // halfway cases at each precision
            values.push((r % 2_000_001) as f64 / 2000.0 - 500.0);
            values.push(((r >> 8) % 1_000_000) as f64 * 0.37 - 1e5);
        }
        for x in values {
            for prec in 0..5 {
                let mut a = String::new();
                push_fixed_short(&mut a, x, prec);
                assert_eq!(a, format!("{:.p$}", x, p = prec as usize), "x = {:e} ({:#x}), prec {}", x, x.to_bits(), prec);
            }
        }
    }
}
