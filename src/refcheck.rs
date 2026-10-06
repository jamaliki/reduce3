//! Developer tool: run the optimizer on the exact post-placement state captured
//! from the original Reduce2 (JSON dumps from ref/harness/dump_ref.py) and
//! compare decisions, scores, coordinates and deletions.

use crate::geom::v3;
use crate::movers::RidingRef;
use crate::optimizer::{self, ConformerInput, OptParams};
use crate::probe::AtomInfo;
use crate::resclass::{self, ResClass};
use crate::world::{AtomLabels, World};
use rustc_hash::FxHashMap;
use serde_json::Value;

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}
fn s(v: &Value) -> String {
    v.as_str().unwrap_or("").to_string()
}

fn info_from(v: &Value) -> AtomInfo {
    let alt = s(&v["altLoc"]);
    AtomInfo {
        vdw_radius: f(&v["vdwRadius"]),
        is_acceptor: v["isAcceptor"].as_bool().unwrap_or(false),
        is_donor: v["isDonor"].as_bool().unwrap_or(false),
        is_dummy_hydrogen: v["isDummyHydrogen"].as_bool().unwrap_or(false),
        is_ion: v["isIon"].as_bool().unwrap_or(false),
        charge: v["charge"].as_i64().unwrap_or(0) as i8,
        alt: alt.bytes().next().unwrap_or(0),
    }
}

pub fn world_from_dump(d: &Value) -> World {
    let atoms = d["after_h_placement"].as_array().unwrap();
    let n = atoms.len();
    let mut labels = Vec::with_capacity(n);
    let mut pos = Vec::with_capacity(n);
    let mut occ = Vec::with_capacity(n);
    let mut b = Vec::with_capacity(n);
    let mut ag = 0u32;
    let mut rg = 0u32;
    let mut prev_ag_key = String::new();
    let mut prev_rg_key = String::new();
    let mut model_index = 0u32;
    let mut prev_model = None;
    for a in atoms {
        let model_id = s(&a["model_id"]);
        if prev_model.as_ref() != Some(&model_id) {
            if prev_model.is_some() {
                model_index += 1;
            }
            prev_model = Some(model_id.clone());
        }
        let rg_key = format!("{}|{}|{}|{}", model_id, s(&a["chain_id"]), s(&a["resseq"]), s(&a["icode"]));
        let ag_key = format!("{}|{}|{}", rg_key, s(&a["altloc"]), s(&a["resname"]));
        if rg_key != prev_rg_key {
            rg += 1;
            prev_rg_key = rg_key;
        }
        if ag_key != prev_ag_key {
            ag += 1;
            prev_ag_key = ag_key;
        }
        let resname_raw = s(&a["resname"]);
        let class = resclass::get_class(&resname_raw);
        labels.push(AtomLabels {
            name: s(&a["name"]).trim().to_ascii_uppercase(),
            raw_name: s(&a["name"]),
            element: s(&a["element"]).trim().to_ascii_uppercase(),
            resname: resname_raw.trim().to_ascii_uppercase().into(),
            resname_raw: resname_raw.clone().into(),
            chain: s(&a["chain_id"]).into(),
            resseq: a["resseq_int"].as_i64().unwrap_or(0) as i32,
            resseq_raw: s(&a["resseq"]).into(),
            icode: s(&a["icode"]).trim().into(),
            altloc: s(&a["altloc"]).into(),
            model_index,
            model_id: model_id.into(),
            ag,
            rg,
            hetero: a["hetero"].as_bool().unwrap_or(false),
            is_water: class == ResClass::CommonWater,
            class: Some(class),
        });
        let x = a["xyz"].as_array().unwrap();
        pos.push(v3(f(&x[0]), f(&x[1]), f(&x[2])));
        occ.push(f(&a["occ"]));
        b.push(f(&a["b"]));
    }
    let mut info = vec![AtomInfo::default(); n];
    let mut bonded = vec![Vec::new(); n];
    for ai in d["atom_info"].as_array().unwrap() {
        let i = ai["i_seq"].as_u64().unwrap() as usize;
        info[i] = info_from(&ai["initial"]);
        bonded[i] = ai["bonded"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect();
    }
    let mut riding = vec![None; n];
    for r in d["riding"].as_array().unwrap() {
        let ih = r["ih"].as_u64().unwrap() as usize;
        riding[ih] = Some(RidingRef { n: r["n"].as_i64().unwrap() as i32, a2: r["a2"].as_i64().unwrap() });
    }
    World { labels, pos, occ, b, info, bonded, riding, n_real: n }
}

pub fn run(path: &str, compat: bool) {
    let text = std::fs::read_to_string(path).expect("read dump");
    let d: Value = serde_json::from_str(&text).expect("parse dump");
    let mut w = world_from_dump(&d);
    let flips = d["options"]["add_flip_movers"].as_bool().unwrap_or(false);
    let rotatable: Vec<u32> = d["riding"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|r| r["htype"].as_str() == Some("alg1b"))
        .map(|r| r["ih"].as_u64().unwrap() as u32)
        .collect();
    let mut p = OptParams::default();
    p.add_flip_movers = flips;
    p.compat = compat;
    // group runs by model
    let runs = d["optimizer_runs"].as_array().unwrap();
    let mut models: Vec<ConformerInput> = Vec::new();
    for r in runs {
        let alt = s(&r["alt"]);
        let mut atoms: Vec<u32> = r["conformer_atom_i_seqs"].as_array().unwrap().iter().map(|x| x.as_u64().unwrap() as u32).collect();
        // The dump lists the atoms after low-quality waters were dropped; put
        // those waters back (GetAtomsForConformer order is i_seq order here).
        if let Some(&first) = atoms.first() {
            let mi = w.labels[first as usize].model_index;
            let present: std::collections::HashSet<u32> = atoms.iter().copied().collect();
            for a in 0..w.n_real as u32 {
                let l = &w.labels[a as usize];
                if l.model_index != mi || present.contains(&a) || !l.is_water || l.element != "O" {
                    continue;
                }
                let bad = !(w.occ[a as usize] >= 0.66 && w.b[a as usize] < 40.0);
                if bad && (l.altloc.is_empty() || *l.altloc == *alt || alt.is_empty()) {
                    atoms.push(a);
                }
            }
            atoms.sort_unstable();
        }
        let mid = s(&r["model_id"]);
        let mi = atoms.first().map(|&a| w.labels[a as usize].model_index as usize).unwrap_or(0);
        if models.last().map(|m| m.model_index) != Some(mi) {
            let model_atoms: Vec<u32> =
                (0..w.n_real as u32).filter(|&a| w.labels[a as usize].model_index as usize == mi).collect();
            let _ = mid;
            models.push(ConformerInput {
                model_index: mi,
                report_model: r["report_model_index"].as_i64().unwrap_or(-1),
                alts: vec![],
                model_atoms,
            });
        }
        models.last_mut().unwrap().alts.push((alt, atoms));
    }
    let nmod = w.labels.iter().map(|l| l.model_index).max().unwrap_or(0) as usize + 1;
    let t0 = std::time::Instant::now();
    let mut debug = optimizer::DebugMovers::default();
    let out = optimizer::optimize_debug(&mut w, &models, &rotatable, &p, nmod, Some(&mut debug));
    // ---- compare Mover coarse positions for the last run
    {
        let dm = d["movers"].as_array().unwrap();
        let last_run = dm.iter().map(|m| m["run"].as_i64().unwrap()).max().unwrap_or(0);
        let theirs: Vec<&Value> = dm.iter().filter(|m| m["run"].as_i64().unwrap() == last_run).collect();
        if theirs.len() != debug.movers.len() {
            println!("MOVER COUNT DIFF {} vs {}", debug.movers.len(), theirs.len());
        }
        for (k, (mine, th)) in debug.movers.iter().zip(theirs.iter()).enumerate() {
            let cps = th["coarse_positions"].as_array().unwrap();
            if cps.len() != mine.coarse_pos.len() {
                println!("MOVER {} ({}) coarse count {} vs {}", k, s(&th["info_at_placement"]), mine.coarse_pos.len(), cps.len());
                if let Some(r) = th["rotator"].as_object() {
                    println!("   theirs angles {:?}", r["coarse_angles"]);
                    println!("   mine angles   {:?}", mine.rot.as_ref().map(|r| r.coarse_angles.clone()));
                }
                continue;
            }
            let mut dev = 0.0f64;
            for (c, st) in cps.iter().enumerate() {
                for (j, x) in st.as_array().unwrap().iter().enumerate() {
                    let x = x.as_array().unwrap();
                    let pt = v3(f(&x[0]), f(&x[1]), f(&x[2]));
                    if j < mine.coarse_pos[c].len() {
                        dev = dev.max(pt.dist(mine.coarse_pos[c][j]));
                    }
                }
            }
            if dev > 1e-6 {
                println!("MOVER {} ({}) coarse positions differ by {:.5}", k, s(&th["info_at_placement"]), dev);
                if let Some(r) = th["rotator"].as_object() {
                    println!("   theirs angles {:?}", r["coarse_angles"]);
                    println!("   mine angles   {:?}", mine.rot.as_ref().map(|r| r.coarse_angles.clone()));
                }
            }
            let fin = &th["final"];
            let mc = debug.final_cfg[k];
            if fin["coarse_index"].as_i64().unwrap() != mc.0 as i64 || fin["fine_index"].as_i64().unwrap() != mc.1 {
                println!(
                    "MOVER {} ({}) final state mine ({}, {}) theirs ({}, {}) score mine {:.4} theirs {:.4}",
                    k,
                    s(&th["info_at_placement"]),
                    mc.0,
                    mc.1,
                    fin["coarse_index"],
                    fin["fine_index"],
                    debug.final_score[k],
                    f(&fin["score"])
                );
            }
        }
    }
    let dt = t0.elapsed().as_secs_f64();

    if std::env::var("REDUCE3_SHOW_CLIQUES").is_ok() {
        let ours: Vec<&str> = out.info.lines().filter(|l| l.contains("Clique optimized")).collect();
        let th = s(&d["info_text"]);
        let theirs_c: Vec<&str> = th.lines().filter(|l| l.contains("Clique optimized")).collect();
        for (i, (a, b)) in ours.iter().zip(theirs_c.iter()).enumerate() {
            println!("CLIQUE {}: mine {} | theirs {}", i, a.trim(), b.trim());
        }
    }
    // ---- compare report lines
    let theirs = s(&d["info_text"]);
    let pick = |t: &str| -> Vec<String> {
        let mut v = Vec::new();
        let mut inside = false;
        for line in t.lines() {
            if line.contains("BEGIN REPORT") {
                inside = true;
            }
            if inside {
                v.push(line.to_string());
            }
            if line.contains("END REPORT") {
                inside = false;
            }
        }
        v
    };
    let a = pick(&out.info);
    let b = pick(&theirs);
    let mut diffs = 0;
    for i in 0..a.len().max(b.len()) {
        let x = a.get(i).map(|s| s.as_str()).unwrap_or("<missing>");
        let y = b.get(i).map(|s| s.as_str()).unwrap_or("<missing>");
        if x != y {
            if diffs < 25 {
                println!("REPORT DIFF line {}:\n  mine:   {}\n  theirs: {}", i, x, y);
            }
            diffs += 1;
        }
    }
    // ---- compare other info lines (placement messages)
    let filt = |t: &str| -> Vec<String> {
        t.lines()
            .filter(|l| l.contains("Added Mover") || l.contains("Could not") || l.contains("Did not") || l.contains("Not adding")
                || l.contains("Set MoverHisFlip") || l.contains("Inserted") || l.contains("Marked") || l.contains("Found ")
                || l.contains("phantom") || l.contains("Ignored") || l.contains("stagger"))
            .filter(|l| !l.contains("Time to"))
            .map(|s| s.to_string())
            .collect()
    };
    let pa = filt(&out.info);
    let pb = filt(&theirs);
    let mut pdiffs = 0;
    for i in 0..pa.len().max(pb.len()) {
        let x = pa.get(i).map(|s| s.as_str()).unwrap_or("<missing>");
        let y = pb.get(i).map(|s| s.as_str()).unwrap_or("<missing>");
        if x != y {
            if pdiffs < 10 {
                println!("PLACEMENT DIFF {}:\n  mine:   {}\n  theirs: {}", i, x, y);
            }
            pdiffs += 1;
        }
    }
    // ---- coordinates after optimization
    let mut maxdev = 0.0f64;
    let mut worst = 0usize;
    let mut nbad = 0;
    for e in d["after_optimization_before_deletion"].as_array().unwrap() {
        let i = e["i_seq"].as_u64().unwrap() as usize;
        let x = e["xyz"].as_array().unwrap();
        let pt = v3(f(&x[0]), f(&x[1]), f(&x[2]));
        let dev = pt.dist(w.pos[i]);
        if dev > 0.001 {
            nbad += 1;
        }
        if dev > maxdev {
            maxdev = dev;
            worst = i;
        }
    }
    // ---- deletions
    let theirs_del: Vec<u32> = d["hydrogens_to_delete"].as_array().unwrap().iter().map(|e| e["i_seq"].as_u64().unwrap() as u32).collect();
    let mut td = theirs_del.clone();
    td.sort_unstable();
    let del_ok = td == out.hydrogens_to_delete;
    println!(
        "{}: report lines {} vs {} ({} differ); placement msgs {} differ; coords: {} atoms > 0.001 A (max {:.4} at {} {} {}); deletions {} ({} vs {}); time {:.3}s",
        path,
        a.len(),
        b.len(),
        diffs,
        pdiffs,
        nbad,
        maxdev,
        worst,
        w.labels.get(worst).map(|l| l.name.as_str()).unwrap_or(""),
        w.labels.get(worst).map(|l| l.resseq).unwrap_or(0),
        if del_ok { "match" } else { "DIFFER" },
        out.hydrogens_to_delete.len(),
        td.len(),
        dt
    );
}

/// Compare hydrogen placement with the dump's `after_h_placement` atoms.
pub fn hcheck(pdb: &str, dump: &str, chem_data: &str, compat: bool) {
    let text = std::fs::read_to_string(pdb).expect("read pdb");
    let mut st = crate::pdbio::read_pdb(&text);
    st.retain_atoms(|a| a.elem() != "X");
    let t0 = std::time::Instant::now();
    let ml = crate::monlib::MonLib::load(std::path::Path::new(chem_data)).expect("monlib");
    let t1 = std::time::Instant::now();
    let p = crate::hplace::HPlaceParams {
        neutron: false,
        n_terminal_charge: crate::hplace::NTermCharge::ResidueOne,
        exclude_water: true,
        keep_existing_h: false,
        adp_scale: 1.0,
        compat,
        cell: crate::cell::processing_cell(&st),
    };
    let placed = crate::hplace::place_hydrogens(&mut st, &ml, &p);
    let t2 = std::time::Instant::now();
    let dtext = std::fs::read_to_string(dump).expect("read dump");
    let d: Value = serde_json::from_str(&dtext).expect("parse dump");
    let theirs = d["after_h_placement"].as_array().unwrap();
    let paths = st.atom_paths();
    let mine: Vec<(String, String, String, i32, String, crate::geom::Vec3)> = paths
        .iter()
        .map(|&pp| {
            let a = st.atom(pp);
            (
                format!("{}#{}", st.models[pp.model as usize].id.trim(), st.chain(pp).id),
                st.atom_group(pp).resname.clone(),
                st.atom_group(pp).altloc.clone(),
                st.residue_group(pp).resseq_as_int(),
                a.name.clone(),
                a.xyz,
            )
        })
        .collect();
    println!(
        "{}: atoms mine {} theirs {}; H mine {} theirs {}; load monlib {:.3}s place {:.3}s; missing residues {:?}",
        pdb,
        mine.len(),
        theirs.len(),
        placed.n_h_final,
        d["h_placement"]["n_H_final"],
        (t1 - t0).as_secs_f64(),
        (t2 - t1).as_secs_f64(),
        placed.no_h_placed
    );
    // match by label
    let key = |c: &str, rn: &str, alt: &str, rs: i32, nm: &str| format!("{}|{}|{}|{}|{}", c, rn.trim(), alt, rs, nm.trim());
    let mut mine_map: FxHashMap<String, (usize, crate::geom::Vec3)> = FxHashMap::default();
    for (k, m) in mine.iter().enumerate() {
        mine_map.insert(key(&m.0, &m.1, &m.2, m.3, &m.4), (k, m.5));
    }
    let mut missing = 0;
    let mut maxdev = 0.0f64;
    let mut ndev = 0;
    let mut order_mismatch = 0;
    let mut shown = 0;
    let mut seen: rustc_hash::FxHashSet<String> = rustc_hash::FxHashSet::default();
    for (k, a) in theirs.iter().enumerate() {
        let kk = key(&format!("{}#{}", s(&a["model_id"]).trim(), s(&a["chain_id"])), &s(&a["resname"]), &s(&a["altloc"]), a["resseq_int"].as_i64().unwrap() as i32, &s(&a["name"]));
        seen.insert(kk.clone());
        let x = a["xyz"].as_array().unwrap();
        let pt = v3(f(&x[0]), f(&x[1]), f(&x[2]));
        match mine_map.get(&kk) {
            None => {
                missing += 1;
                if shown < 15 {
                    println!("  missing in mine: {}", kk);
                    shown += 1;
                }
            }
            Some(&(mk, mp)) => {
                if mk != k {
                    order_mismatch += 1;
                }
                let dv = mp.dist(pt);
                if mp == pt {
                    BITEXACT.with(|c| c.set(c.get() + 1));
                } else if std::env::var_os("R3BITS").is_some() {
                    println!("  inexact {} {:?} vs {:?}", kk, mp, pt);
                }
                if dv > 0.002 {
                    ndev += 1;
                    if shown < 15 {
                        println!("  deviation {:.3} at {}", dv, kk);
                        shown += 1;
                    }
                }
                maxdev = maxdev.max(dv);
            }
        }
    }
    println!("  bit-exact coordinates: {} of {}", BITEXACT.with(|c| c.get()), mine.len());
    // riding types per H (by label)
    {
        let mut theirs_type: FxHashMap<String, String> = FxHashMap::default();
        for r in d["riding"].as_array().unwrap() {
            let ih = r["ih"].as_u64().unwrap() as usize;
            let a = &theirs[ih];
            let kk = key(&format!("{}#{}", s(&a["model_id"]).trim(), s(&a["chain_id"])), &s(&a["resname"]), &s(&a["altloc"]), a["resseq_int"].as_i64().unwrap() as i32, &s(&a["name"]));
            let lab = |i: i64| -> String {
                if i < 0 {
                    return "-".into();
                }
                let a = &theirs[i as usize];
                format!("{}{}", s(&a["name"]).trim(), a["resseq_int"].as_i64().unwrap())
            };
            theirs_type.insert(
                kk,
                format!(
                    "{} n={} a0={} a1={} a2={}",
                    s(&r["htype"]),
                    r["n"].as_i64().unwrap(),
                    lab(r["a0"].as_i64().unwrap()),
                    lab(r["a1"].as_i64().unwrap()),
                    lab(r["a2"].as_i64().unwrap())
                ),
            );
        }
        let mut ndiff = 0;
        for (k, m) in mine.iter().enumerate() {
            let kk = key(&m.0, &m.1, &m.2, m.3, &m.4);
            let mlab = |i: i64| -> String {
                if i < 0 {
                    return "-".into();
                }
                let m = &mine[i as usize];
                format!("{}{}", m.4.trim(), m.3)
            };
            if std::env::var_os("R3COEF").is_some() {
                if let Some(c) = placed.riding[k] {
                    let th = d["riding"].as_array().unwrap().iter().find(|r| {
                        let a = &theirs[r["ih"].as_u64().unwrap() as usize];
                        key(&format!("{}#{}", s(&a["model_id"]).trim(), s(&a["chain_id"])), &s(&a["resname"]), &s(&a["altloc"]), a["resseq_int"].as_i64().unwrap() as i32, &s(&a["name"])) == kk
                    });
                    if let Some(r) = th {
                        let (ta, tb, thh, td) = (f(&r["a"]), f(&r["b"]), f(&r["h"]), f(&r["disth"]));
                        let same_set = c.a == ta && c.b == tb && c.h == thh && c.disth == td;
                        let swapped = c.a == tb && c.b == ta;
                        if !same_set && !swapped {
                            println!("  coef differs {} {}: mine a={:?} b={:?} h={:?} d={:?} theirs a={:?} b={:?} h={:?} d={:?}", kk, c.htype.name(), c.a, c.b, c.h, c.disth, ta, tb, thh, td);
                        }
                    }
                }
            }
            let mt = placed.riding[k].map(|c| {
                format!("{} n={} a0={} a1={} a2={}", c.htype.name(), c.n, mlab(c.a0 as i64), mlab(c.a1 as i64), mlab(c.a2))
            });
            let tt = theirs_type.get(&kk).cloned();
            // a1/a2/a3 order of the symmetric types follows Python set order and
            // does not change positions: compare those as sets
            let norm = |t: &Option<String>| -> Option<String> {
                t.as_ref().map(|t| {
                    let f: Vec<&str> = t.split(' ').collect();
                    if matches!(f[0], "3neigbs" | "2neigbs" | "flat_2neigbs" | "2tetra") {
                        format!("{} {} {}", f[0], f[1], f[2])
                    } else {
                        t.clone()
                    }
                })
            };
            if norm(&mt) != norm(&tt) && (mt.is_some() || tt.is_some()) {
                if ndiff < 10 {
                    println!("  riding type differs at {}: mine {:?} theirs {:?}", kk, mt, tt);
                }
                ndiff += 1;
            }
        }
        if ndiff > 0 {
            println!("  riding type differences: {}", ndiff);
        }
    }
    let extra: Vec<&String> = mine_map.keys().filter(|k| !seen.contains(*k)).collect();
    for e in extra.iter().take(10) {
        println!("  extra in mine: {}", e);
    }
    println!(
        "  missing {} extra {} order-mismatch {} deviating(>0.002) {} max dev {:.4}",
        missing,
        extra.len(),
        order_mismatch,
        ndev,
        maxdev
    );
}

thread_local! {
    static BITEXACT: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Compare the optimizer inputs Reduce3 builds itself (bonded lists and
/// extra atom info after reinterpretation) with the dump.
pub fn wcheck(pdb: &str, dump: &str, chem_data: &str) {
    let text = std::fs::read_to_string(pdb).expect("read pdb");
    let mut st = crate::pdbio::read_pdb(&text);
    let ml = crate::monlib::MonLib::load(std::path::Path::new(chem_data)).expect("monlib");
    let d: Value = serde_json::from_str(&std::fs::read_to_string(dump).expect("read dump")).expect("parse");
    let hp = crate::hplace::HPlaceParams {
        neutron: false,
        n_terminal_charge: crate::hplace::NTermCharge::ResidueOne,
        exclude_water: true,
        keep_existing_h: false,
        adp_scale: 1.0,
        compat: true,
        cell: crate::cell::processing_cell(&st),
    };
    let _placed = crate::hplace::place_hydrogens(&mut st, &ml, &hp);
    st.sort_atoms_in_place();
    st.reset_i_seq();
    let flat = crate::interp::FlatAtoms::from_structure(&st);
    let it = crate::interp::interpret(&st, &flat, &ml, &crate::interp::InterpParams { neutron: false, link_distance_cutoff: 3.0, compat: true, auto_comps: Default::default() });
    let n = flat.pos.len();
    let bonded = crate::atominfo::bonded_lists(n, it.bonds.iter().map(|b| (b.i, b.j)));
    let ex = crate::atominfo::extra_atom_info(&st, &flat, &it.etype, &ml, &bonded, true);
    let theirs = world_from_dump(&d);
    println!("{}: atoms mine {} theirs {}", pdb, n, theirs.labels.len());
    let (mut nb_set, mut nb_order, mut ninfo) = (0, 0, 0);
    for k in 0..n.min(theirs.labels.len()) {
        let mut a = bonded[k].clone();
        let mut b = theirs.bonded[k].clone();
        if a != b {
            a.sort_unstable();
            b.sort_unstable();
            if a != b {
                if nb_set < 10 {
                    println!("  bonded set differs at {} {} {}: mine {:?} theirs {:?}", k, theirs.labels[k].resname, theirs.labels[k].name, bonded[k], theirs.bonded[k]);
                }
                nb_set += 1;
            } else {
                if nb_order < 10 {
                    println!("  bonded order differs at {} {} {}: mine {:?} theirs {:?}", k, theirs.labels[k].resname, theirs.labels[k].name, bonded[k], theirs.bonded[k]);
                }
                nb_order += 1;
            }
        }
        let (x, y) = (&ex.info[k], &theirs.info[k]);
        if x.vdw_radius != y.vdw_radius || x.is_acceptor != y.is_acceptor || x.is_donor != y.is_donor || x.is_ion != y.is_ion || x.charge != y.charge || x.alt != y.alt {
            if ninfo < 10 {
                println!("  info differs at {} {} {} {}: mine {:?} theirs {:?}", k, theirs.labels[k].resname, theirs.labels[k].resseq, theirs.labels[k].name, x, y);
            }
            ninfo += 1;
        }
    }
    let tw = s(&d["initial_extra_atom_info_warnings"]);
    println!("  bonded set diffs {} order diffs {} info diffs {} warnings equal {}", nb_set, nb_order, ninfo, tw == ex.warnings);
    if tw != ex.warnings {
        for (i, (a, b)) in ex.warnings.lines().zip(tw.lines()).enumerate() {
            if a != b {
                println!("  first warning diff at line {}: mine {:?} theirs {:?}", i, a, b);
                break;
            }
        }
        println!("  warning lines mine {} theirs {}", ex.warnings.lines().count(), tw.lines().count());
    }
}

/// Compare the CCD-derived dictionaries with `reference/harness/dump_ccd_restraints.py`.
pub fn ccdcheck(chem_data: &str, jsonl: &str) {
    let ml = crate::monlib::MonLib::load(std::path::Path::new(chem_data)).expect("chem_data");
    let text = std::fs::read_to_string(jsonl).expect("jsonl");
    let (mut n, mut same) = (0, 0);
    let mut problems: std::collections::BTreeMap<&str, Vec<String>> = Default::default();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let d: serde_json::Value = serde_json::from_str(line).expect("json");
        let id = d["id"].as_str().unwrap().to_string();
        n += 1;
        let mine = ml.ccd_comp(&id, true);
        let ok = d["ok"].as_bool().unwrap();
        let mut bad: Vec<(&str, String)> = Vec::new();
        match (&mine, ok) {
            (None, false) => {}
            (Some(_), false) => bad.push(("accepted, Reduce2 rejects", id.clone())),
            (None, true) => bad.push(("rejected, Reduce2 accepts", id.clone())),
            (Some(c), true) => {
                let s = |v: &serde_json::Value| v.as_str().unwrap_or("").to_string();
                let f = |v: &serde_json::Value| v.as_f64().or_else(|| v.as_str().and_then(|s| s.parse::<f64>().ok())).unwrap();
                let same = |x: f64, y: f64| x == y || (x.is_nan() && y.is_nan());
                let atoms: Vec<(String, String)> = d["atoms"].as_array().unwrap().iter().map(|a| (s(&a[0]), s(&a[1]))).collect();
                let mine_atoms: Vec<(String, String)> = c.atoms.iter().map(|a| (a.id.clone(), a.type_symbol.clone())).collect();
                if atoms != mine_atoms {
                    bad.push(("atoms", id.clone()));
                }
                let bonds = d["bonds"].as_array().unwrap();
                if bonds.len() != c.bonds.len() {
                    bad.push(("bond count", id.clone()));
                } else {
                    for (b, m) in bonds.iter().zip(&c.bonds) {
                        if (s(&b[0]), s(&b[1]), s(&b[2])) != (m.a1.clone(), m.a2.clone(), m.type_.clone()) {
                            bad.push(("bond atoms", format!("{} {}-{}", id, m.a1, m.a2)));
                        } else if !same(f(&b[3]), m.value_dist.unwrap()) {
                            bad.push(("bond value", format!("{} {}-{} {} vs {}", id, m.a1, m.a2, f(&b[3]), m.value_dist.unwrap())));
                        }
                    }
                }
                let angles = d["angles"].as_array().unwrap();
                if angles.len() != c.angles.len() {
                    bad.push(("angle count", id.clone()));
                } else {
                    for (a, m) in angles.iter().zip(&c.angles) {
                        if (s(&a[0]), s(&a[1]), s(&a[2])) != (m.a1.clone(), m.a2.clone(), m.a3.clone()) {
                            bad.push(("angle order", id.clone()));
                            break;
                        } else if !same(f(&a[3]), m.value.unwrap()) {
                            bad.push(("angle value", format!("{} {}-{}-{} {} vs {}", id, m.a1, m.a2, m.a3, f(&a[3]), m.value.unwrap())));
                        }
                    }
                }
                let tors = d["tors"].as_array().unwrap();
                if tors.len() != c.tors.len() {
                    bad.push(("torsion count", id.clone()));
                } else {
                    for (t, m) in tors.iter().zip(&c.tors) {
                        let names = [s(&t[1]), s(&t[2]), s(&t[3]), s(&t[4])];
                        if names != m.a || s(&t[0]) != m.id {
                            bad.push(("torsion order", id.clone()));
                            break;
                        }
                        let (x, y) = (f(&t[5]), m.value.unwrap());
                        if !same(x, y) {
                            let kind = if x.abs() == 180.0 && y.abs() == 180.0 { "torsion +-180" } else { "torsion value" };
                            bad.push((kind, format!("{} {} {} vs {}", id, m.id, x, y)));
                        }
                    }
                }
            }
        }
        if bad.is_empty() {
            same += 1;
        }
        for (k, v) in bad {
            problems.entry(k).or_default().push(v);
        }
    }
    println!("{} entries, {} identical", n, same);
    for (k, v) in problems {
        println!("  {}: {} (e.g. {})", k, v.len(), v.iter().take(4).cloned().collect::<Vec<_>>().join("; "));
    }
}
