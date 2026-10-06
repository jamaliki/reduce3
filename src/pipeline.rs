//! The Reduce3 program flow (`mmtbx/programs/reduce2.py` `Program.run`):
//! add hydrogens, reinterpret, optimize movable groups, delete hydrogens the
//! optimizer rejects, and produce the output model and description text.

use crate::atominfo;
use crate::hplace::{self, HPlaceParams, NTermCharge};
use crate::interp::{self, FlatAtoms, InterpParams};
use crate::model::Structure;
use crate::monlib::MonLib;
use crate::movers::RidingRef;
use crate::optimizer::{self, ConformerInput, OptParams};
use crate::riding::HType;
use crate::world::World;
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Approach {
    Add,
    Remove,
    Optimize,
}

#[derive(Clone, Debug)]
pub struct Params {
    pub approach: Approach,
    /// Reproduce Reduce2 exactly, including its bugs. `run` applies it to
    /// `opt.compat` as well.
    pub compat: bool,
    pub keep_existing_h: bool,
    pub n_terminal_charge: NTermCharge,
    pub exclude_water: bool,
    pub model_id: Option<usize>,
    pub ignore_missing_restraints: bool,
    pub stop_on_any_missing_hydrogen: bool,
    pub opt: OptParams,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            approach: Approach::Add,
            compat: false,
            keep_existing_h: false,
            n_terminal_charge: NTermCharge::ResidueOne,
            exclude_water: true,
            model_id: None,
            ignore_missing_restraints: false,
            stop_on_any_missing_hydrogen: false,
            opt: OptParams::default(),
        }
    }
}

pub struct Output {
    pub structure: Structure,
    /// The description file text (without the program/version header).
    pub description: String,
    /// Messages for the terminal.
    pub log: String,
}

fn res_name_and_id(st: &Structure, p: crate::model::AtomPath) -> String {
    let ag = st.atom_group(p);
    let rg = st.residue_group(p);
    let ch = st.chain(p);
    format!("chain {} {}{} {}{}", ch.id, ag.altloc, ag.resname.trim().to_ascii_uppercase(), rg.resseq_as_int(), rg.icode.trim())
}

/// `_RemoveModelsExceptIndex` followed by taking the first remaining model.
fn select_model(st: &mut Structure, model_id: usize, compat: bool) -> Result<(), String> {
    if model_id == 0 {
        return Err("Model ID must be >=1 if specified (None means all models)".into());
    }
    let n = st.models.len();
    let keep = if compat {
        // Reduce2 indexes models() with the 1-based id and keeps everything
        // when it is out of range, then uses the first remaining model.
        if model_id < n { model_id } else { 0 }
    } else {
        if model_id > n {
            return Err(format!("Model ID {} is out of range: the file has {} model(s)", model_id, n));
        }
        model_id - 1
    };
    let m = st.models.swap_remove(keep);
    st.models = vec![m];
    Ok(())
}

/// Run the whole program on a structure.
pub fn run(mut st: Structure, ml: &MonLib, p: &Params) -> Result<Output, String> {
    // `compat` governs the optimizer too
    let p = &{
        let mut p = p.clone();
        p.opt.compat = p.compat;
        p
    };
    let mut desc = String::new();
    let mut log = String::new();
    // element X atoms are dropped first
    st.retain_atoms(|a| a.elem() != "X");
    if let Some(mid) = p.model_id {
        select_model(&mut st, mid, p.compat)?;
    }
    st.reset_i_seq();
    let cell = if p.compat { crate::cell::processing_cell(&st) } else { None };

    match p.approach {
        Approach::Add | Approach::Optimize => {
            let t_add = Instant::now();
            let mut prior: Option<Prior> = None;
            if p.approach == Approach::Add {
                let hp = HPlaceParams {
                    neutron: p.opt.use_neutron_distances,
                    n_terminal_charge: p.n_terminal_charge,
                    exclude_water: p.exclude_water,
                    keep_existing_h: p.keep_existing_h,
                    adp_scale: 1.0,
                    compat: p.compat,
                    cell: cell.clone(),
                };
                let placed = hplace::place_hydrogens(&mut st, ml, &hp);
                log += &placed.log;
                if !p.ignore_missing_restraints && !placed.no_h_placed.is_empty() {
                    let mut bad: Vec<String> = placed.no_h_placed.clone();
                    bad.dedup();
                    return Err(format!("Restraints were not found for the following residues: {}", bad.join(" ")));
                }
                if p.stop_on_any_missing_hydrogen && !placed.site_labels_no_para.is_empty() {
                    return Err(format!(
                        "Insufficient restraints were found for the following atoms:{}",
                        placed.site_labels_no_para.join(",")
                    ));
                }
                if !st.has_hydrogens() {
                    return Err("It was not possible to place any H atoms. Is this a single atom model?".into());
                }
                prior = Some(Prior { bonds: placed.bonds, etype: placed.etype, riding: placed.riding });
            }
            let add_time = t_add.elapsed().as_secs_f64();
            let t_opt = Instant::now();
            if prior.is_none() {
                // _ReinterpretModel sorts before processing; after placement the
                // model already has its restraints and is used as it is
                st.sort_atoms_in_place();
                st.reset_serial();
                st.reset_i_seq();
                if let Some(uc) = &cell {
                    hplace::round_off_like_cctbx(&mut st, uc);
                }
            }
            let (opt_info, deletes) = optimize_structure(&mut st, ml, p, prior)?;
            let opt_time = t_opt.elapsed().as_secs_f64();
            desc += &opt_info;
            if !deletes.is_empty() {
                desc += " Deleting hydrogens requested for deletion by optimization:\n";
                let paths = st.atom_paths();
                for &k in &deletes {
                    let pth = paths[k as usize];
                    desc += &format!("  Deleting {} {}\n", res_name_and_id(&st, pth), st.atom(pth).name.trim().to_ascii_uppercase());
                }
                let del: rustc_hash::FxHashSet<usize> = deletes.iter().map(|&k| k as usize).collect();
                let mut k = 0usize;
                st.retain_atoms(|_| {
                    let keep = !del.contains(&k);
                    k += 1;
                    keep
                });
            }
            if p.approach == Approach::Add {
                desc += &format!("Time to Add Hydrogen = {:.3} sec\n", add_time);
            }
            desc += &format!("Time to Optimize = {:.3} sec\n", opt_time);
        }
        Approach::Remove => {
            st.retain_atoms(|a| a.elem() != "H");
        }
    }
    st.sort_atoms_in_place();
    st.reset_serial();
    st.reset_i_seq();
    Ok(Output { structure: st, description: desc, log })
}

/// Restraint data that hydrogen placement leaves on the model. Reduce2 keeps
/// it for the optimizer instead of interpreting the model again.
pub struct Prior {
    pub bonds: Vec<(u32, u32, u16)>,
    pub etype: Vec<interp::EType>,
    pub riding: Vec<Option<crate::riding::RidingCoef>>,
}

/// Build the optimizer's view of the model, optimize, and write the new
/// coordinates back. Returns the optimizer's text and the flat indices of
/// hydrogens to delete. With `prior` (after placement) its bonds, energy
/// types and riding parameterization are used; without it the model is
/// interpreted here and a riding manager is set up, which idealizes the
/// hydrogens first.
pub fn optimize_structure(st: &mut Structure, ml: &MonLib, p: &Params, prior: Option<Prior>) -> Result<(String, Vec<u32>), String> {
    let mut info = String::new();
    let mut tm = Instant::now();
    let v = p.opt.verbosity;
    let timing = |what: &str, tm: &mut Instant| -> String {
        let s = if v >= 2 { format!("  Time to {}: {:.3}\n", what, tm.elapsed().as_secs_f64()) } else { String::new() };
        *tm = Instant::now();
        s
    };
    let flat = FlatAtoms::from_structure(st);
    let n = flat.pos.len();
    let (it, prior_riding) = match prior {
        Some(pr) => {
            let mut it = interp::Interp::empty(n);
            for &(i, j, o) in &pr.bonds {
                it.add_bond(i, j, 0.0, o);
            }
            it.etype = pr.etype;
            (it, Some(pr.riding))
        }
        None => {
            // _ReinterpretModel: process the model with its hydrogens
            let it = interp::interpret(st, &flat, ml, &InterpParams {
                    neutron: p.opt.use_neutron_distances,
                    link_distance_cutoff: 3.0,
                    compat: p.compat,
                    auto_comps: Default::default(),
                });
            (it, None)
        }
    };
    info += &timing("get coordinates", &mut tm);
    info += &timing("compute bond proxies", &mut tm);
    let bonded = atominfo::bonded_lists(n, it.bonds.iter().map(|b| (b.i, b.j)));
    info += &timing("compute bonded neighbor lists", &mut tm);
    let ex = atominfo::extra_atom_info(st, &flat, &it.etype, ml, &bonded, p.opt.probe.set_polar_hydrogen_radius);
    info += &ex.warnings;
    info += &timing("get extra atom info", &mut tm);

    let coefs = match prior_riding {
        Some(c) => c,
        None => {
            // setup_riding_h_manager(): parameterize and idealize
            let mut sites = flat.pos.clone();
            let occ: Vec<f64> = flat.path.iter().map(|&pp| st.atom(pp).occ).collect();
            let rg_of: Vec<u64> =
                flat.path.iter().map(|pp| ((pp.model as u64) << 40) | ((pp.chain as u64) << 20) | pp.rg as u64).collect();
            let resnames: Vec<String> = flat.path.iter().map(|&pp| st.atom_group(pp).resname.trim().to_string()).collect();
            let table = hplace::expected_heavy_table(ml, &Default::default(), &resnames);
            let expected_heavy = |a: u32| -> Option<usize> {
                table.get(&resnames[a as usize])?.as_ref()?.get(flat.name[a as usize].trim()).copied()
            };
            let ra = crate::riding::RidingAtoms {
                is_h: &flat.is_h,
                altloc: &flat.altloc,
                name: &flat.name,
                occ: &occ,
                rg_of: &rg_of,
                expected_heavy: &expected_heavy,
                dictionary_nh2_torsion: !p.compat,
            };
            let rr = crate::riding::riding(&it, &mut sites, &ra, false);
            for (k, pth) in flat.path.iter().enumerate() {
                st.atom_mut(*pth).xyz = sites[k];
            }
            rr.coef
        }
    };
    let mut rotatable: Vec<u32> = Vec::new();
    let mut riding: Vec<Option<RidingRef>> = vec![None; n];
    for (k, c) in coefs.iter().enumerate() {
        if let Some(c) = c {
            riding[k] = Some(RidingRef { n: c.n, a2: c.a2 });
            if c.htype == HType::Alg1b {
                rotatable.push(k as u32);
            }
        }
    }
    info += &timing("select rotatable hydrogens", &mut tm);

    crate::model::mem_checkpoint("world: before");
    let mut w = World::from_structure(st);
    crate::model::mem_checkpoint("world: built");
    w.info = ex.info;
    w.bonded = bonded;
    w.riding = riding;

    // models and alternates to run
    let n_models = st.models.len();
    let model_range: Vec<usize> = if p.compat { vec![n_models - 1] } else { (0..n_models).collect() };
    let paths = st.atom_paths();
    let mut models: Vec<ConformerInput> = Vec::new();
    for mi in model_range {
        let m = &st.models[mi];
        let mut alts: Vec<String> = vec![String::new()];
        for c in &m.chains {
            for a in c.conformer_altlocs() {
                if !alts.contains(&a) {
                    alts.push(a);
                }
            }
        }
        if alts.len() > 1 {
            alts.retain(|a| !a.is_empty() && a != " ");
        }
        alts.sort();
        alts.reverse();
        if let Some(a) = &p.opt.alt_id {
            alts = vec![a.clone()];
        }
        let model_atoms: Vec<u32> = (0..n as u32).filter(|&k| paths[k as usize].model as usize == mi).collect();
        // atoms of each chain (in atom order) and each chain's conformers
        let mut chain_atoms: Vec<Vec<u32>> = vec![Vec::new(); m.chains.len()];
        for &k in &model_atoms {
            chain_atoms[paths[k as usize].chain as usize].push(k);
        }
        let chain_confs: Vec<Vec<String>> = m.chains.iter().map(|c| c.conformer_altlocs()).collect();
        let mut runs = Vec::new();
        for alt in alts {
            // GetAtomsForConformer
            let mut atoms = Vec::new();
            for (ci, confs) in chain_confs.iter().enumerate() {
                let which = (1..confs.len()).find(|&i| confs[i] == alt).unwrap_or(0);
                let Some(conf) = confs.get(which) else { continue };
                for &k in &chain_atoms[ci] {
                    let al = &st.atom_group(paths[k as usize]).altloc;
                    if al.is_empty() || al == conf {
                        atoms.push(k);
                    }
                }
            }
            runs.push((alt, atoms));
        }
        info += &timing("compute alternates", &mut tm);
        models.push(ConformerInput {
            model_index: mi,
            report_model: if p.compat { -1 } else { mi as i64 },
            alts: runs,
            model_atoms,
        });
    }
    crate::model::mem_checkpoint("optimizer: start");
    let out = optimizer::optimize(&mut w, &models, &rotatable, &p.opt, n_models);
    crate::model::mem_checkpoint("optimizer: done");
    info += &out.info;
    // write the coordinates back
    for (k, pth) in paths.iter().enumerate() {
        st.atom_mut(*pth).xyz = w.pos[k];
    }
    Ok((info, out.hydrogens_to_delete))
}
