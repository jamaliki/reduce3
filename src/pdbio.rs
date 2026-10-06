//! PDB-format reading and writing with iotbx.pdb hierarchy-construction rules.

use crate::geom::v3;
use crate::model::*;
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Raw atom record with its labels, prior to hierarchy construction.
#[derive(Clone, Debug)]
pub(crate) struct RawAtom {
    pub atom: Atom,
    pub altloc: String,
    pub resname: String,
    pub chain: String,
    pub resseq: String,
    pub icode: String,
}

fn col(line: &str, a: usize, b: usize) -> &str {
    // 1-based inclusive columns, tolerant of short lines
    let bytes = line.as_bytes();
    let s = (a - 1).min(bytes.len());
    let e = b.min(bytes.len());
    if s >= e { "" } else { &line[s..e] }
}

fn pad_right(s: &str, n: usize) -> String {
    format!("{:<n$}", s, n = n)
}

/// Build an iotbx hierarchy from raw atoms grouped by model and chain ranges.
/// `chains` lists, per model, the atom index ranges forming chains;
/// `breaks` holds indices of atoms that follow a BREAK record.
pub(crate) fn build_hierarchy(
    raw: Vec<RawAtom>,
    model_ids: Vec<String>,
    chains: Vec<Vec<(usize, usize)>>,
    breaks: &[usize],
) -> Vec<Model> {
    let mut raw: Vec<Option<RawAtom>> = raw.into_iter().map(Some).collect();
    let mut models = Vec::with_capacity(model_ids.len());
    for (mi, mid) in model_ids.into_iter().enumerate() {
        let mut model = Model { id: mid, chains: Vec::new() };
        for &(cb, ce) in &chains[mi] {
            if cb >= ce {
                continue;
            }
            let chain_id = raw[cb].as_ref().unwrap().chain.clone();
            let mut chain = Chain { id: chain_id, residue_groups: Vec::new() };
            let mut rg_start = cb;
            let mut prev: Option<(String, String)> = None; // (resid, resname)
            let mut open_run_has_blank = false;
            let mut link_to_previous = false;
            let mut i = cb;
            while i < ce {
                let ra = raw[i].as_ref().unwrap();
                let resid = format!("{}{}", ra.resseq, ra.icode);
                let resname = ra.resname.clone();
                let cur_blank = ra.altloc.is_empty() || ra.altloc == " ";
                let first_after_break = breaks.contains(&i);
                if let Some((presid, presname)) = &prev {
                    let mut boundary = *presid != resid;
                    if !boundary && *presname != resname {
                        if open_run_has_blank || cur_blank {
                            boundary = true;
                        } else {
                            for j in (i + 1)..ce {
                                let f = raw[j].as_ref().unwrap();
                                if f.resname != resname {
                                    break;
                                }
                                if format!("{}{}", f.resseq, f.icode) != resid {
                                    break;
                                }
                                if f.altloc.is_empty() || f.altloc == " " {
                                    boundary = true;
                                    break;
                                }
                            }
                        }
                    }
                    if boundary {
                        let rg = make_residue_group(&mut raw, rg_start, i, link_to_previous);
                        chain.residue_groups.push(rg);
                        rg_start = i;
                        link_to_previous = !first_after_break;
                        open_run_has_blank = false;
                    }
                }
                prev = Some((resid, resname));
                if cur_blank {
                    open_run_has_blank = true;
                }
                i += 1;
            }
            if prev.is_some() {
                let rg = make_residue_group(&mut raw, rg_start, ce, link_to_previous);
                chain.residue_groups.push(rg);
            }
            merge_disconnected_pure_altloc(&mut chain);
            model.chains.push(chain);
        }
        models.push(model);
    }
    models
}

fn make_residue_group(raw: &mut [Option<RawAtom>], b: usize, e: usize, link: bool) -> ResidueGroup {
    let first = raw[b].as_ref().unwrap();
    let resseq = first.resseq.clone();
    let icode = first.icode.clone();
    // Group atoms by confid (altloc + resname), ordered by first appearance.
    let mut groups: Vec<(String, String, Vec<Atom>)> = Vec::new();
    for slot in raw[b..e].iter_mut() {
        let ra = slot.take().unwrap();
        let alt = if ra.altloc.is_empty() { " ".to_string() } else { ra.altloc.clone() };
        match groups.iter_mut().find(|g| g.0 == alt && g.1 == ra.resname) {
            Some(g) => g.2.push(ra.atom),
            None => groups.push((alt, ra.resname.clone(), vec![ra.atom])),
        }
    }
    let mut rg = ResidueGroup {
        resseq,
        icode,
        link_to_previous: link,
        atom_groups: groups
            .into_iter()
            .map(|(alt, resname, atoms)| AtomGroup { altloc: alt, resname, atoms })
            .collect(),
    };
    edit_blank_altloc(&mut rg);
    rg
}

/// iotbx `residue_group::edit_blank_altloc()`: blank-altloc atom groups go to
/// the front with altloc "", and blank atoms whose names also occur in an
/// alternate of the same residue name move to a separate " " atom group.
fn edit_blank_altloc(rg: &mut ResidueGroup) {
    let (mut blank, mut alt): (Vec<AtomGroup>, Vec<AtomGroup>) =
        rg.atom_groups.drain(..).partition(|ag| ag.altloc == " " || ag.altloc.is_empty());
    if blank.is_empty() {
        rg.atom_groups = alt;
        return;
    }
    for ag in &mut blank {
        ag.altloc = String::new();
    }
    // names in blank groups per resname that also occur in alternates
    let mut blank_but_alt: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for a in &alt {
        if let Some(bg) = blank.iter().find(|g| g.resname == a.resname) {
            for at in &a.atoms {
                if bg.atoms.iter().any(|x| x.name == at.name) {
                    blank_but_alt.entry(a.resname.clone()).or_default().push(at.name.clone());
                }
            }
        }
    }
    let mut extra: Vec<AtomGroup> = Vec::new();
    if !blank_but_alt.is_empty() {
        for bg in &mut blank {
            if let Some(names) = blank_but_alt.get(&bg.resname) {
                let (moved, kept): (Vec<Atom>, Vec<Atom>) =
                    bg.atoms.drain(..).partition(|a| names.contains(&a.name));
                bg.atoms = kept;
                if !moved.is_empty() {
                    extra.push(AtomGroup { altloc: " ".into(), resname: bg.resname.clone(), atoms: moved });
                }
            }
        }
        blank.retain(|g| !g.atoms.is_empty());
    }
    blank.extend(extra);
    blank.append(&mut alt);
    rg.atom_groups = blank;
}

/// iotbx `chain.merge_disconnected_residue_groups_with_pure_altloc()`: residue
/// groups sharing resseq+icode whose atom groups all carry distinct, non-blank
/// altlocs are merged into the first of them.
fn merge_disconnected_pure_altloc(chain: &mut Chain) {
    let n = chain.residue_groups.len();
    let mut by_id: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, rg) in chain.residue_groups.iter().enumerate() {
        by_id.entry((rg.resseq.clone(), rg.icode.clone())).or_default().push(i);
    }
    let mut remove = vec![false; n];
    for (_, idx) in by_id {
        if idx.len() < 2 {
            continue;
        }
        let mut seen: Vec<&str> = vec!["", " "];
        let mut ok = true;
        'outer: for &i in &idx {
            for ag in &chain.residue_groups[i].atom_groups {
                if seen.contains(&ag.altloc.as_str()) {
                    ok = false;
                    break 'outer;
                }
                seen.push(&ag.altloc);
            }
        }
        if !ok {
            continue;
        }
        let first = idx[0];
        for &i in &idx[1..] {
            let ags = std::mem::take(&mut chain.residue_groups[i].atom_groups);
            chain.residue_groups[first].atom_groups.extend(ags);
            remove[i] = true;
        }
    }
    if remove.iter().any(|&r| r) {
        let mut k = 0;
        chain.residue_groups.retain(|_| {
            let keep = !remove[k];
            k += 1;
            keep
        });
    }
}

fn parse_f(s: &str) -> f64 {
    s.trim().parse::<f64>().unwrap_or(0.0)
}

/// Normalize a chain id taken from PDB columns 21-22.
fn chain_from_cols(c2: &str) -> String {
    let s = format!("{:>2}", c2);
    s.trim_start().to_string()
}

pub fn read_pdb(text: &str) -> Structure {
    let mut raw: Vec<RawAtom> = Vec::new();
    let mut model_ids: Vec<String> = Vec::new();
    let mut chains: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut breaks: Vec<usize> = Vec::new();
    let mut st = Structure::default();

    // chain tracking per model
    let mut cur_chain_start: Option<usize> = None;
    let mut prev_chain_segid: Option<(String, String)> = None;
    let mut model_chains: Vec<(usize, usize)> = Vec::new();
    let mut in_model = false;
    let mut segids_seen: Vec<String> = Vec::new();
    let mut chain_splits_ignoring_segid: Vec<(usize, usize)> = Vec::new();
    let mut ignore_segid_start: Option<usize> = None;

    let close_chain = |start: &mut Option<usize>, end: usize, list: &mut Vec<(usize, usize)>| {
        if let Some(s) = start.take() {
            if end > s {
                list.push((s, end));
            }
        }
    };

    for line in text.lines() {
        let rec = col(line, 1, 6);
        match rec {
            "ATOM  " | "HETATM" | "ATOM" | "HETATM " => {
                let name = pad_right(col(line, 13, 16), 4);
                let altloc = col(line, 17, 17).to_string();
                let altloc = if altloc == " " { String::new() } else { altloc };
                let resname = col(line, 18, 20).to_string();
                let chain = chain_from_cols(col(line, 21, 22));
                let resseq = format!("{:>4}", col(line, 23, 26));
                let icode = {
                    let c = col(line, 27, 27);
                    if c.is_empty() { " ".to_string() } else { c.to_string() }
                };
                let x = parse_f(col(line, 31, 38));
                let y = parse_f(col(line, 39, 46));
                let z = parse_f(col(line, 47, 54));
                let occ_s = col(line, 55, 60);
                let occ = if occ_s.trim().is_empty() { 1.0 } else { parse_f(occ_s) };
                let b = parse_f(col(line, 61, 66));
                let segid = col(line, 73, 76).to_string();
                let mut element = col(line, 77, 78).to_ascii_uppercase();
                let charge = format!("{:<2}", col(line, 79, 80));
                if element.trim().is_empty() {
                    element = infer_element(&name);
                }
                let atom = Atom {
                    name,
                    element: format!("{:>2}", element.trim()),
                    charge,
                    serial: col(line, 7, 11).to_string(),
                    xyz: v3(x, y, z),
                    occ,
                    b,
                    segid: segid.clone(),
                    hetero: rec.starts_with("HETATM"),
                    uij: None,
                    i_seq: 0,
                    src: raw.len() as u32,
                };
                let idx = raw.len();
                // chain boundary detection
                let key = (chain.clone(), segid.clone());
                match &prev_chain_segid {
                    None => {
                        cur_chain_start = Some(idx);
                        ignore_segid_start = Some(idx);
                    }
                    Some((pc, ps)) => {
                        if *pc != chain {
                            close_chain(&mut cur_chain_start, idx, &mut model_chains);
                            close_chain(&mut ignore_segid_start, idx, &mut chain_splits_ignoring_segid);
                            cur_chain_start = Some(idx);
                            ignore_segid_start = Some(idx);
                        } else if chain.len() <= 1 && *ps != segid {
                            close_chain(&mut cur_chain_start, idx, &mut model_chains);
                            cur_chain_start = Some(idx);
                        }
                    }
                }
                if segids_seen.last() != Some(&segid) {
                    segids_seen.push(segid.clone());
                }
                prev_chain_segid = Some(key);
                in_model = true;
                st.pdb_records.serials.push(col(line, 7, 11).to_string());
                raw.push(RawAtom { atom, altloc, resname, chain, resseq, icode });
            }
            "ANISOU" => {
                if let Some(last) = raw.last_mut() {
                    let vals: Vec<f64> = (0..6)
                        .map(|k| col(line, 29 + 7 * k, 35 + 7 * k).trim().parse::<f64>().unwrap_or(0.0) * 1e-4)
                        .collect();
                    last.atom.uij = Some([vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]]);
                }
            }
            "MODEL " | "MODEL" => {
                let id = col(line, 11, 14).trim().to_string();
                model_ids.push(if id.is_empty() { col(line, 7, 80).trim().to_string() } else { id });
            }
            "ENDMDL" => {
                close_chain(&mut cur_chain_start, raw.len(), &mut model_chains);
                close_chain(&mut ignore_segid_start, raw.len(), &mut chain_splits_ignoring_segid);
                prev_chain_segid = None;
                finish_model_chains(&mut chains, &mut model_chains, &mut chain_splits_ignoring_segid, &mut segids_seen);
                in_model = false;
            }
            "TER   " | "TER" => {
                close_chain(&mut cur_chain_start, raw.len(), &mut model_chains);
                close_chain(&mut ignore_segid_start, raw.len(), &mut chain_splits_ignoring_segid);
                prev_chain_segid = None;
            }
            "BREAK " | "BREAK" => breaks.push(raw.len()),
            "CRYST1" => {
                st.had_cell_record = true;
                let a = parse_f(col(line, 7, 15));
                let b = parse_f(col(line, 16, 24));
                let c = parse_f(col(line, 25, 33));
                let al = parse_f(col(line, 34, 40));
                let be = parse_f(col(line, 41, 47));
                let ga = parse_f(col(line, 48, 54));
                let sg = col(line, 56, 66).trim().to_string();
                st.crystal = Some(CrystalSymmetry { cell: [a, b, c, al, be, ga], space_group: sg });
            }
            "LINK  " | "LINKR " => {
                let l1 = AtomLabel {
                    name: col(line, 13, 16).trim().to_string(),
                    altloc: col(line, 17, 17).trim().to_string(),
                    resname: col(line, 18, 20).trim().to_string(),
                    chain: chain_from_cols(col(line, 21, 22)),
                    resseq: col(line, 23, 26).trim().to_string(),
                    icode: col(line, 27, 27).trim().to_string(),
                };
                let l2 = AtomLabel {
                    name: col(line, 43, 46).trim().to_string(),
                    altloc: col(line, 47, 47).trim().to_string(),
                    resname: col(line, 48, 50).trim().to_string(),
                    chain: chain_from_cols(col(line, 51, 52)),
                    resseq: col(line, 53, 56).trim().to_string(),
                    icode: col(line, 57, 57).trim().to_string(),
                };
                let d = col(line, 74, 78).trim().parse::<f64>().ok();
                st.links.push(LinkRecord { atom1: l1, atom2: l2, distance: d, kind: "LINK".into() });
                st.pdb_records.links.push(line.to_string());
            }
            "SSBOND" => {
                let l1 = AtomLabel {
                    name: "SG".into(),
                    altloc: String::new(),
                    resname: col(line, 12, 14).trim().to_string(),
                    chain: chain_from_cols(col(line, 15, 16)),
                    resseq: col(line, 18, 21).trim().to_string(),
                    icode: col(line, 22, 22).trim().to_string(),
                };
                let l2 = AtomLabel {
                    name: "SG".into(),
                    altloc: String::new(),
                    resname: col(line, 26, 28).trim().to_string(),
                    chain: chain_from_cols(col(line, 29, 30)),
                    resseq: col(line, 32, 35).trim().to_string(),
                    icode: col(line, 36, 36).trim().to_string(),
                };
                let d = col(line, 74, 78).trim().parse::<f64>().ok();
                st.links.push(LinkRecord { atom1: l1, atom2: l2, distance: d, kind: "SSBOND".into() });
                st.pdb_records.links.push(line.to_string());
            }
            "REMARK" | "END   " | "END" | "MASTER" | "ATOM 1" => {}
            _ => {
                if rec.starts_with("CONECT") {
                    st.pdb_records.conect.push(line.to_string());
                }
                if !line.trim().is_empty() {
                    // any other record type ends the current chain (iotbx transition)
                    if prev_chain_segid.is_some() && !in_model_atom_record(rec) {
                        close_chain(&mut cur_chain_start, raw.len(), &mut model_chains);
                        close_chain(&mut ignore_segid_start, raw.len(), &mut chain_splits_ignoring_segid);
                        prev_chain_segid = None;
                    }
                }
            }
        }
    }
    let _ = in_model;
    close_chain(&mut cur_chain_start, raw.len(), &mut model_chains);
    close_chain(&mut ignore_segid_start, raw.len(), &mut chain_splits_ignoring_segid);
    if !model_chains.is_empty() || chains.len() < model_ids.len() || model_ids.is_empty() {
        if !model_chains.is_empty() || model_ids.is_empty() {
            finish_model_chains(&mut chains, &mut model_chains, &mut chain_splits_ignoring_segid, &mut segids_seen);
        }
    }
    if model_ids.is_empty() {
        model_ids.push(String::new());
    }
    while chains.len() < model_ids.len() {
        chains.push(Vec::new());
    }
    chains.truncate(model_ids.len());
    st.models = build_hierarchy(raw, model_ids, chains, &breaks);
    st.reset_i_seq();
    st
}

fn in_model_atom_record(rec: &str) -> bool {
    matches!(rec, "SIGATM" | "SIGUIJ")
}

fn finish_model_chains(
    chains: &mut Vec<Vec<(usize, usize)>>,
    model_chains: &mut Vec<(usize, usize)>,
    ignoring_segid: &mut Vec<(usize, usize)>,
    segids: &mut Vec<String>,
) {
    // If segids are not unique across consecutive runs, chains are not split on segid.
    let mut set = std::collections::HashSet::new();
    let dup = segids.iter().any(|s| !set.insert(s.clone()));
    if dup {
        chains.push(std::mem::take(ignoring_segid));
        model_chains.clear();
    } else {
        chains.push(std::mem::take(model_chains));
        ignoring_segid.clear();
    }
    segids.clear();
}

/// Element from the atom name when the element columns are blank (iotbx
/// `determine_chemical_element_simple`, simplified).
fn infer_element(name: &str) -> String {
    let b = name.as_bytes();
    if b.len() < 2 {
        return String::new();
    }
    let c0 = b[0] as char;
    let c1 = b[1] as char;
    if c0 == ' ' || c0.is_ascii_digit() {
        c1.to_string()
    } else if c0 == 'H' && name.trim().len() == 4 {
        "H".to_string()
    } else {
        format!("{}{}", c0, c1)
    }
}

// ----------------------------------------------------------------------------
// Writing

/// Fractionalization matrix for a unit cell (PDB/cctbx orthogonalization convention).
/// cos of an angle in degrees, exact for multiples of 90 (as cctbx.uctbx does).
fn cos_deg(d: f64) -> f64 {
    let r = d.rem_euclid(360.0);
    if r == 90.0 || r == 270.0 { 0.0 } else if r == 0.0 { 1.0 } else if r == 180.0 { -1.0 } else { d.to_radians().cos() }
}
fn sin_deg(d: f64) -> f64 {
    let r = d.rem_euclid(360.0);
    if r == 0.0 || r == 180.0 { 0.0 } else if r == 90.0 { 1.0 } else if r == 270.0 { -1.0 } else { d.to_radians().sin() }
}

pub fn fractionalization_matrix(cell: &[f64; 6]) -> [[f64; 3]; 3] {
    let [a, b, c, al, be, ga] = *cell;
    let (ca, cb, cg) = (cos_deg(al), cos_deg(be), cos_deg(ga));
    let sg = sin_deg(ga);
    let v = (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt();
    [
        [1.0 / a, -cg / (a * sg), (ca * cg - cb) / (a * v * sg)],
        [0.0, 1.0 / (b * sg), (cb * cg - ca) / (b * v * sg)],
        [0.0, 0.0, sg / (c * v)],
    ]
}

pub fn format_atom_line(out: &mut String, a: &Atom, ag: &AtomGroup, rg: &ResidueGroup, chain_id: &str) {
    format_atom_line_serial(out, a, &a.serial, ag, rg, chain_id)
}

/// An ATOM/HETATM (+ANISOU) record with the given serial field.
pub fn format_atom_line_serial(out: &mut String, a: &Atom, serial: &str, ag: &AtomGroup, rg: &ResidueGroup, chain_id: &str) {
    let start = out.len();
    out.push_str(if a.hetero { "HETATM" } else { "ATOM  " });
    let _ = write!(out, "{:>5} ", serial);
    let _ = write!(out, "{:<4}", a.name);
    let alt = if ag.altloc.is_empty() { " " } else { &ag.altloc };
    let _ = write!(out, "{:<1}", alt);
    let _ = write!(out, "{:>3}", ag.resname);
    let _ = write!(out, "{:>2}", chain_id);
    let _ = write!(out, "{:>4}", rg.resseq);
    let _ = write!(out, "{:<1}", if rg.icode.is_empty() { " " } else { &rg.icode });
    out.push_str("   ");
    for k in 0..3 {
        crate::fastfmt::push_fixed(out, a.xyz[k].clamp(-1e7, 1e8), 3, 8);
    }
    crate::fastfmt::push_fixed(out, a.occ.clamp(-1e5, 1e6), 2, 6);
    crate::fastfmt::push_fixed(out, a.b.clamp(-1e5, 1e6), 2, 6);
    let _ = write!(out, "      {:<4}{:>2}{:<2}", a.segid, a.element.trim(), a.charge.trim_end());
    // right-trim blanks
    while out.len() > start && out.ends_with(' ') {
        out.pop();
    }
    out.push('\n');
    if let Some(u) = a.uij {
        let ls = out.len();
        out.push_str("ANISOU");
        let _ = write!(out, "{:>5} ", serial);
        let _ = write!(out, "{:<4}{:<1}{:>3}{:>2}{:>4}{:<1}", a.name, alt, ag.resname, chain_id, rg.resseq,
            if rg.icode.is_empty() { " " } else { &rg.icode });
        out.push(' ');
        for k in 0..6 {
            // "%7.0f": round half to even, and a tiny negative prints "-0"
            crate::fastfmt::push_fixed(out, (u[k] * 10000.0).clamp(-1.0e7, 1.0e8), 0, 7);
        }
        let _ = write!(out, "  {:<4}{:>2}{:<2}", a.segid, a.element.trim(), a.charge.trim_end());
        while out.len() > ls && out.ends_with(' ') {
            out.pop();
        }
        out.push('\n');
    }
}

pub fn write_pdb(st: &Structure, write_cryst: bool) -> String {
    let mut out = String::with_capacity(st.atoms_size() * 82 + 1024);
    if write_cryst {
        if let Some(cs) = &st.crystal {
            let c = cs.cell;
            let _ = writeln!(
                out,
                "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {}",
                c[0], c[1], c[2], c[3], c[4], c[5], cs.space_group
            );
            let f = fractionalization_matrix(&c);
            for (i, row) in f.iter().enumerate() {
                let _ = writeln!(out, "SCALE{}    {:10.6}{:10.6}{:10.6}     {:10.5}", i + 1,
                    row[0] + 0.0, row[1] + 0.0, row[2] + 0.0, 0.0);
            }
        }
    }
    let nm = st.models.len();
    for m in &st.models {
        if nm != 1 {
            let _ = writeln!(out, "MODEL     {:>4}", m.id);
        }
        // serial numbers restart in every model, as iotbx writes them
        let mut serial = 1i64;
        for ch in &m.chains {
            for rg in &ch.residue_groups {
                for ag in &rg.atom_groups {
                    for a in &ag.atoms {
                        let sn = crate::model::hy36_encode(5, serial).unwrap_or_else(|| "*****".into());
                        format_atom_line_serial(&mut out, a, &sn, ag, rg, &ch.id);
                        serial += 1;
                    }
                }
            }
            if ch.is_polymer_chain() {
                out.push_str("TER\n");
            }
        }
        if nm != 1 {
            out.push_str("ENDMDL\n");
        }
    }
    out.push_str("END\n");
    out
}

/// `write_pdb` plus the input's SSBOND and LINK records (before CRYST1) and its
/// CONECT records, renumbered to the output serials of the first model;
/// partners that are no longer there are left out.
pub fn write_pdb_preserving(st: &Structure) -> String {
    let body = write_pdb(st, true);
    let recs = &st.pdb_records;
    if recs.links.is_empty() && recs.conect.is_empty() {
        return body;
    }
    let mut out = String::with_capacity(body.len() + 81 * (recs.links.len() + recs.conect.len()));
    for l in &recs.links {
        out.push_str(l.trim_end());
        out.push('\n');
    }
    let mut conect = String::new();
    if !recs.conect.is_empty() {
        let mut src_of: std::collections::HashMap<&str, u32> = std::collections::HashMap::new();
        for (k, s) in recs.serials.iter().enumerate() {
            src_of.entry(s.trim()).or_insert(k as u32);
        }
        let mut serial_of: std::collections::HashMap<u32, String> = std::collections::HashMap::new();
        if let Some(m) = st.models.first() {
            let mut serial = 1i64;
            for a in m.chains.iter().flat_map(|c| &c.residue_groups).flat_map(|r| &r.atom_groups).flat_map(|g| &g.atoms) {
                if a.src != Atom::NEW {
                    serial_of.insert(a.src, crate::model::hy36_encode(5, serial).unwrap_or_else(|| "*****".into()));
                }
                serial += 1;
            }
        }
        let new_serial = |field: &str| src_of.get(field.trim()).and_then(|k| serial_of.get(k)).cloned();
        for line in &recs.conect {
            let fields: Vec<&str> = (0..5).filter_map(|k| line.get(6 + 5 * k..(11 + 5 * k).min(line.len()))).collect();
            let Some(base) = fields.first().and_then(|f| new_serial(f)) else { continue };
            let partners: Vec<String> = fields[1..].iter().filter(|f| !f.trim().is_empty()).filter_map(|f| new_serial(f)).collect();
            if partners.is_empty() {
                continue;
            }
            let _ = write!(conect, "CONECT{:>5}", base);
            for p in partners {
                let _ = write!(conect, "{:>5}", p);
            }
            conect.push('\n');
        }
    }
    match body.strip_suffix("END\n") {
        Some(head) => {
            out.push_str(head);
            out.push_str(&conect);
            out.push_str("END\n");
        }
        None => {
            out.push_str(&body);
            out.push_str(&conect);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conect_records_follow_the_new_serials() {
        let txt = "SSBOND   1 CYS A    1    CYS A    1                          1555   2555  2.03  \n\
ATOM     10  N   CYS A   1       1.000   1.000   1.000  1.00 10.00           N  \n\
ATOM     11  SG  CYS A   1       2.000   1.000   1.000  1.00 10.00           S  \n\
ATOM     12  H   CYS A   1       0.500   1.000   1.000  1.00 10.00           H  \n\
HETATM   20 ZN    ZN A 101       4.000   1.000   1.000  1.00 10.00          ZN  \n\
CONECT   11   20\n\
CONECT   12   10\n\
END\n";
        let mut st = read_pdb(txt);
        // the H goes, which shifts nothing before it but drops its CONECT
        st.retain_atoms(|a| a.name.trim() != "H");
        let out = write_pdb_preserving(&st);
        assert!(out.starts_with("SSBOND   1 CYS A    1"));
        let conect: Vec<&str> = out.lines().filter(|l| l.starts_with("CONECT")).collect();
        assert_eq!(conect, ["CONECT    2    3"]);
        assert!(out.ends_with("CONECT    2    3\nEND\n"));
    }

    #[test]
    fn roundtrip_simple() {
        let txt = "CRYST1   40.960   18.650   22.520  90.00  90.77  90.00 P 1 21 1      2\n\
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N  \n\
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C  \n\
ATOM      3  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N  \n\
END\n";
        let st = read_pdb(txt);
        assert_eq!(st.models.len(), 1);
        assert_eq!(st.models[0].chains[0].residue_groups.len(), 2);
        let out = write_pdb(&st, true);
        assert!(out.contains("ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N\n"));
        assert!(out.contains("SCALE1      0.024414  0.000000  0.000328        0.00000"));
    }

    #[test]
    fn altloc_groups() {
        let txt = "ATOM      1  N   HIS A   1      26.965  32.911   7.593  1.00  7.19           N\n\
ATOM      2  N  AHIS A   2      26.965  32.911   7.593  1.00  7.19           N\n\
ATOM      3  N  BHIS A   2      26.965  32.911   7.593  1.00  7.19           N\n\
ATOM      4  CA  HIS A   2      26.965  32.911   7.593  1.00  7.19           C\n";
        let st = read_pdb(txt);
        let rg = &st.models[0].chains[0].residue_groups[1];
        assert_eq!(rg.atom_groups.len(), 3);
        assert_eq!(rg.atom_groups[0].altloc, "");
        assert_eq!(rg.atom_groups[1].altloc, "A");
        assert_eq!(st.models[0].chains[0].conformer_altlocs(), vec!["A".to_string(), "B".to_string()]);
    }
}
