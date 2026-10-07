//! Reduce3: a fast Rust reimplementation of cctbx Reduce2 (hydrogen addition
//! and optimization of rotatable and flippable groups).

use reduce3::{hplace, mmcif, model, monlib, pdbio, pipeline};
#[cfg(feature = "refcheck")]
use reduce3::refcheck;
use pipeline::{Approach, Params};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

const USAGE: &str = "\
Reduce3: add hydrogens to a macromolecular model and optimize the rotatable
and flippable groups (a Rust reimplementation of cctbx Reduce2).

Usage: reduce3 [options] model.pdb [name=value ...]
       reduce3 [options] --out-dir DIR model1.cif model2.pdb ... [name=value ...]
       reduce3 [options] --out-dir DIR --batch list.txt [name=value ...]

Batch mode (several models, or --batch): the monomer library is loaded once and
models are processed in parallel, each on one thread. A model that fails is
reported and skipped; the exit status is nonzero if any failed. Inputs may be
gzip-compressed (model.cif.gz); outputs are written uncompressed.

Options:
  --compat              reproduce Reduce2 exactly, including its known bugs
  --chem-data DIR       chem_data directory (monomer library, CCD). Default:
                        $REDUCE3_CHEM_DATA, $CHEM_DATA, or the active conda env
  -o, --output FILE     output model (same as output.filename=FILE)
  --threads N           worker threads (default: all cores)
  --batch FILE          read model paths from FILE, one per line (- for stdin)
  --out-dir DIR         write outputs to DIR (batch mode; default: .)
  --jobs N              models processed at once in batch mode (default: all cores)
  --no-description      do not write the description (report) file
  -q, --quiet           no progress messages
  -h, --help            this help
  -V, --version         print the version

Reduce2 parameters (name=value), with the same defaults:
  approach=add|remove|optimize     add_flip_movers=False
  n_terminal_charge=residue_one|first_in_chain|no_charge
  keep_existing_H=False            exclude_water=True
  use_neutron_distances=False      preference_magnitude=1.0
  non_flip_preference=0.5          skip_bond_fix_up=False
  set_flip_states=None             model_id=None   alt_id=None
  bonded_neighbor_depth=4          verbosity=2
  stop_on_any_missing_hydrogen=False   ignore_missing_restraints=False
  output.filename=<input>H.pdb (FH with flips)
  output.description_file_name=<output>.txt   output.write_files=True
  probe.probe_radius=0.25 ... (all probe.* scoring parameters)
";

fn parse_bool(v: &str) -> Result<bool, String> {
    match v.to_ascii_lowercase().as_str() {
        "true" | "1" | "yes" => Ok(true),
        "false" | "0" | "no" => Ok(false),
        _ => Err(format!("not a boolean: {}", v)),
    }
}
fn parse_f(v: &str) -> Result<f64, String> {
    v.parse::<f64>().map_err(|_| format!("not a number: {}", v))
}
fn none_or(v: &str) -> Option<String> {
    if v.eq_ignore_ascii_case("none") || v.is_empty() { None } else { Some(v.to_string()) }
}

struct Cli {
    inputs: Vec<PathBuf>,
    batch_list: Option<String>,
    out_dir: Option<PathBuf>,
    jobs: Option<usize>,
    no_description: bool,
    output: Option<PathBuf>,
    description: Option<PathBuf>,
    write_files: bool,
    chem_data: Option<String>,
    threads: Option<usize>,
    quiet: bool,
    params: Params,
}

fn parse_args(args: &[String]) -> Result<Cli, String> {
    let mut cli = Cli {
        inputs: Vec::new(),
        batch_list: None,
        out_dir: None,
        jobs: None,
        no_description: false,
        output: None,
        description: None,
        write_files: true,
        chem_data: None,
        threads: None,
        quiet: false,
        params: Params::default(),
    };
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        let next = |i: &mut usize| -> Result<String, String> {
            *i += 1;
            args.get(*i).cloned().ok_or_else(|| format!("{} needs a value", args[*i - 1]))
        };
        match a.as_str() {
            "-h" | "--help" => return Err(String::new()),
            "-V" | "--version" => {
                println!("reduce3 {}", env!("CARGO_PKG_VERSION"));
                std::process::exit(0);
            }
            "--compat" => cli.params.compat = true,
            "--fixed" => cli.params.compat = false,
            "-q" | "--quiet" => cli.quiet = true,
            "--chem-data" => cli.chem_data = Some(next(&mut i)?),
            "-o" | "--output" => cli.output = Some(PathBuf::from(next(&mut i)?)),
            "--threads" => cli.threads = Some(next(&mut i)?.parse().map_err(|_| "bad --threads value".to_string())?),
            "--batch" => cli.batch_list = Some(next(&mut i)?),
            "--out-dir" => cli.out_dir = Some(PathBuf::from(next(&mut i)?)),
            "--jobs" => cli.jobs = Some(next(&mut i)?.parse().map_err(|_| "bad --jobs value".to_string())?),
            "--no-description" => cli.no_description = true,
            _ if a.contains('=') && !a.starts_with('-') => {
                let (k, v) = a.split_once('=').unwrap();
                let key = k.trim().trim_start_matches("reduce2.");
                let v = v.trim();
                let p = &mut cli.params;
                match key {
                    "approach" => {
                        p.approach = match v {
                            "add" => Approach::Add,
                            "remove" => Approach::Remove,
                            "optimize" => Approach::Optimize,
                            _ => return Err(format!("approach must be add, remove or optimize, not {}", v)),
                        }
                    }
                    "keep_existing_H" => p.keep_existing_h = parse_bool(v)?,
                    "n_terminal_charge" => {
                        p.n_terminal_charge = match v {
                            "residue_one" => hplace::NTermCharge::ResidueOne,
                            "first_in_chain" => hplace::NTermCharge::FirstInChain,
                            "no_charge" => hplace::NTermCharge::NoCharge,
                            _ => return Err(format!("bad n_terminal_charge: {}", v)),
                        }
                    }
                    "exclude_water" => p.exclude_water = parse_bool(v)?,
                    "use_neutron_distances" => p.opt.use_neutron_distances = parse_bool(v)?,
                    "preference_magnitude" => p.opt.preference_magnitude = parse_f(v)?,
                    "alt_id" => p.opt.alt_id = none_or(v),
                    "model_id" => {
                        p.model_id = match none_or(v) {
                            None => None,
                            Some(s) => Some(s.parse().map_err(|_| format!("bad model_id: {}", s))?),
                        }
                    }
                    "add_flip_movers" => p.opt.add_flip_movers = parse_bool(v)?,
                    "non_flip_preference" => p.opt.non_flip_preference = parse_f(v)?,
                    "skip_bond_fix_up" => p.opt.skip_bond_fixup = parse_bool(v)?,
                    "set_flip_states" => p.opt.flip_states = none_or(v).unwrap_or_default(),
                    "verbosity" => p.opt.verbosity = v.parse().map_err(|_| format!("bad verbosity: {}", v))?,
                    "bonded_neighbor_depth" => {
                        p.opt.bonded_neighbor_depth = v.parse().map_err(|_| format!("bad bonded_neighbor_depth: {}", v))?
                    }
                    "stop_on_any_missing_hydrogen" => p.stop_on_any_missing_hydrogen = parse_bool(v)?,
                    "ignore_missing_restraints" => p.ignore_missing_restraints = parse_bool(v)?,
                    "output.filename" | "filename" => cli.output = none_or(v).map(PathBuf::from),
                    "output.description_file_name" | "description_file_name" => {
                        cli.description = none_or(v).map(PathBuf::from)
                    }
                    "output.write_files" | "write_files" => cli.write_files = parse_bool(v)?,
                    "overwrite" | "output.overwrite" | "profile" | "output.print_atom_info" | "output.flipkin_directory"
                    | "output.clique_outline_file_name" | "comparison_file" => {}
                    _ if key.starts_with("probe.") => {
                        let pr = &mut p.opt.probe;
                        match &key[6..] {
                            "probe_radius" => pr.probe_radius = parse_f(v)?,
                            "density" => pr.density = parse_f(v)?,
                            "worse_clash_cutoff" => pr.worse_clash_cutoff = parse_f(v)?,
                            "clash_cutoff" => pr.clash_cutoff = parse_f(v)?,
                            "contact_cutoff" => pr.contact_cutoff = parse_f(v)?,
                            "uncharged_hydrogen_cutoff" => pr.uncharged_hydrogen_cutoff = parse_f(v)?,
                            "charged_hydrogen_cutoff" => pr.charged_hydrogen_cutoff = parse_f(v)?,
                            "bump_weight" => pr.bump_weight = parse_f(v)?,
                            "hydrogen_bond_weight" => pr.hydrogen_bond_weight = parse_f(v)?,
                            "gap_weight" => pr.gap_weight = parse_f(v)?,
                            "allow_weak_hydrogen_bonds" => pr.allow_weak_hydrogen_bonds = parse_bool(v)?,
                            "ignore_ion_interactions" => pr.ignore_ion_interactions = parse_bool(v)?,
                            "set_polar_hydrogen_radius" => pr.set_polar_hydrogen_radius = parse_bool(v)?,
                            other => return Err(format!("unknown parameter: probe.{}", other)),
                        }
                    }
                    _ => return Err(format!("unknown parameter: {}", key)),
                }
            }
            _ if a.starts_with('-') && a.len() > 1 => return Err(format!("unknown option: {}", a)),
            _ => cli.inputs.push(PathBuf::from(a)),
        }
        i += 1;
    }
    if cli.inputs.is_empty() && cli.batch_list.is_none() {
        return Err("no input model given".to_string());
    }
    cli.params.opt.compat = cli.params.compat;
    Ok(cli)
}

/// Seconds since the epoch as "YYYY-MM-DD HH:MM:SS" (UTC).
fn utc_timestamp() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0) as i64;
    let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    // civil_from_days (Howard Hinnant)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + if m <= 2 { 1 } else { 0 };
    format!("{:04}-{:02}-{:02} {:02}:{:02}:{:02} UTC", y, m, d, rem / 3600, (rem % 3600) / 60, rem % 60)
}

/// The input's name without a trailing `.gz` (outputs are written uncompressed).
fn model_name(input: &Path) -> PathBuf {
    let name = input.to_string_lossy();
    match name.len().checked_sub(3) {
        Some(k) if name.is_char_boundary(k) && name[k..].eq_ignore_ascii_case(".gz") => PathBuf::from(&name[..k]),
        _ => input.to_path_buf(),
    }
}

fn default_output(input: &Path, flips: bool) -> PathBuf {
    let input = model_name(input);
    let stem = input.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "model".into());
    let ext = input.extension().map(|s| format!(".{}", s.to_string_lossy())).unwrap_or_default();
    PathBuf::from(format!("{}{}{}", stem, if flips { "FH" } else { "H" }, ext))
}

/// A model file's text; gzip-compressed files (by their magic bytes) are
/// decompressed.
fn read_model_text(path: &Path) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {}: {}", path.display(), e))?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        #[cfg(feature = "gzip")]
        {
            let mut text = String::new();
            std::io::Read::read_to_string(&mut flate2::read::MultiGzDecoder::new(&bytes[..]), &mut text)
                .map_err(|e| format!("cannot decompress {}: {}", path.display(), e))?;
            return Ok(text);
        }
        #[cfg(not(feature = "gzip"))]
        return Err(format!("{} is gzip-compressed; this build reads uncompressed files only", path.display()));
    }
    String::from_utf8(bytes).map_err(|e| format!("cannot read {}: {}", path.display(), e))
}

/// One model to process and where its results go.
struct Job {
    input: PathBuf,
    output: PathBuf,
    description: Option<PathBuf>,
    /// False for `output.write_files=False`: run, write nothing.
    write: bool,
}

/// Inputs this large are read again for writing rather than kept parsed
/// through the run (peak memory of the largest entries).
const KEEP_PARSED_BELOW: usize = 32 << 20;

/// Read, process and write one model; returns its atom count.
fn process_one(job: &Job, ml: &monlib::MonLib, params: &Params, header: &str) -> Result<usize, String> {
    let text = read_model_text(&job.input)?;
    let lower = model_name(&job.input).to_string_lossy().to_ascii_lowercase();
    let is_cif = lower.ends_with(".cif") || lower.ends_with(".mmcif") || text.trim_start().starts_with("data_");
    let out_lower = job.output.to_string_lossy().to_ascii_lowercase();
    let write_pdb_format = out_lower.ends_with(".pdb") || out_lower.ends_with(".ent");
    // fixed mode writes mmCIF back into its source block: parse it once
    let keep_doc = is_cif && !write_pdb_format && !params.compat && text.len() < KEEP_PARSED_BELOW;
    model::mem_checkpoint("start");
    let source_block = |doc: &'_ reduce3::cif::Document<'_>| -> Result<usize, String> {
        doc.blocks.iter().position(|b| b.category("_atom_site").is_some()).ok_or_else(|| "no _atom_site loop found in the mmCIF file".to_string())
    };
    // the text is kept only when its parsed block is written back
    let (st, kept) = if keep_doc {
        (None, Some(text))
    } else {
        (Some(if is_cif { mmcif::read_mmcif(&text)? } else { pdbio::read_pdb(&text) }), None)
    };
    let doc = kept.as_deref().map(reduce3::cif::parse);
    let st = match (st, &doc) {
        (Some(st), _) => st,
        (None, Some(d)) => mmcif::structure_from_cif(&d.blocks[source_block(d)?])?,
        (None, None) => unreachable!("a model or its parsed text"),
    };
    model::mem_checkpoint("read");
    let n_atoms = st.atoms_size();
    if n_atoms == 0 {
        return Err(format!("no atoms found in {}", job.input.display()));
    }
    let out = pipeline::run(st, ml, params)?;
    model::mem_checkpoint("pipeline done");
    if !job.write {
        return Ok(n_atoms);
    }
    let model_text = if write_pdb_format && params.compat {
        pdbio::write_pdb(&out.structure, true)
    } else if write_pdb_format {
        // keep SSBOND, LINK and CONECT (fixed mode)
        pdbio::write_pdb_preserving(&out.structure)
    } else if is_cif && !params.compat {
        // keep every category of the input block (fixed mode)
        let write = |doc: &reduce3::cif::Document<'_>, size: usize| -> Result<String, String> {
            let block = &doc.blocks[source_block(doc)?];
            let mut w = reduce3::cifsource::CifText::with_capacity(size + out.structure.atoms_size() * 100);
            mmcif::write_cif_preserving(&out.structure, block, block.name, &mut w)?;
            Ok(w.out)
        };
        match &doc {
            Some(d) => write(d, kept.as_ref().map_or(0, |t| t.len()))?,
            None => {
                let text = read_model_text(&job.input)?;
                let d = reduce3::cif::parse(&text);
                write(&d, text.len())?
            }
        }
    } else {
        mmcif::write_mmcif(&out.structure)
    };
    std::fs::write(&job.output, model_text).map_err(|e| format!("cannot write {}: {}", job.output.display(), e))?;
    if let Some(description) = &job.description {
        let mut desc = String::with_capacity(header.len() + out.description.len());
        desc.push_str(header);
        desc.push_str(&out.description);
        std::fs::write(description, desc).map_err(|e| format!("cannot write {}: {}", description.display(), e))?;
    }
    Ok(n_atoms)
}

fn description_path(output: &Path) -> PathBuf {
    PathBuf::from(output.to_string_lossy().replace(".pdb", ".txt").replace(".cif", ".txt"))
}

/// The description file's first lines: version, time and command line.
fn description_header() -> String {
    let mut h = format!("reduce3 v.{}, run {}\n", env!("CARGO_PKG_VERSION"), utc_timestamp());
    for a in std::env::args() {
        h.push(' ');
        h.push_str(&a);
    }
    h.push('\n');
    h
}

fn run_cli(args: &[String]) -> Result<(), String> {
    let mut cli = parse_args(args)?;
    let say = |m: &str| {
        if !cli.quiet {
            eprintln!("{}", m);
        }
    };
    let root = monlib::MonLib::locate(cli.chem_data.as_deref()).ok_or_else(|| {
        "could not find chem_data (the monomer library). Pass --chem-data DIR or set REDUCE3_CHEM_DATA.".to_string()
    })?;
    if let Some(list) = &cli.batch_list {
        let text = if list == "-" {
            let mut s = String::new();
            std::io::Read::read_to_string(&mut std::io::stdin(), &mut s).map_err(|e| format!("cannot read the model list: {}", e))?;
            s
        } else {
            std::fs::read_to_string(list).map_err(|e| format!("cannot read {}: {}", list, e))?
        };
        cli.inputs.extend(text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')).map(PathBuf::from));
    }
    let batch = cli.inputs.len() > 1 || cli.batch_list.is_some() || cli.out_dir.is_some();
    if !batch {
        if let Some(n) = cli.threads {
            let _ = rayon::ThreadPoolBuilder::new().num_threads(n.max(1)).build_global();
        }
        let t0 = std::time::Instant::now();
        let input = cli.inputs[0].clone();
        let output = cli.output.clone().unwrap_or_else(|| default_output(&input, cli.params.opt.add_flip_movers));
        let description = (!cli.no_description).then(|| cli.description.clone().unwrap_or_else(|| description_path(&output)));
        let job = Job { input, output, description, write: cli.write_files };
        let ml = monlib::MonLib::load(&root)?;
        let n = process_one(&job, &ml, &cli.params, &description_header())?;
        say(&format!("Read {} ({} atoms)", job.input.display(), n));
        if job.write {
            let desc = job.description.as_ref().map(|d| format!(" and {}", d.display())).unwrap_or_default();
            say(&format!("Wrote {}{} ({:.3} s)", job.output.display(), desc, t0.elapsed().as_secs_f64()));
        }
        return Ok(());
    }
    if cli.output.is_some() || cli.description.is_some() {
        return Err("output.filename and output.description_file_name name one model's files; use --out-dir in batch mode".to_string());
    }
    let out_dir = cli.out_dir.clone().unwrap_or_else(|| PathBuf::from("."));
    std::fs::create_dir_all(&out_dir).map_err(|e| format!("cannot create {}: {}", out_dir.display(), e))?;
    let flips = cli.params.opt.add_flip_movers;
    let mut seen: std::collections::HashMap<PathBuf, PathBuf> = std::collections::HashMap::new();
    let mut jobs: Vec<Job> = Vec::with_capacity(cli.inputs.len());
    for input in &cli.inputs {
        let name = default_output(input, flips);
        let output = out_dir.join(&name);
        let description = (!cli.no_description).then(|| out_dir.join(description_path(&name)));
        for path in std::iter::once(&output).chain(description.as_ref()) {
            if let Some(other) = seen.insert(path.clone(), input.clone()) {
                return Err(format!(
                    "{} and {} would both write {} (run them separately, or use --no-description if only the reports collide)",
                    other.display(),
                    input.display(),
                    path.display()
                ));
            }
        }
        jobs.push(Job { input: input.clone(), output, description, write: cli.write_files });
    }
    let n_jobs = cli.jobs.unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)).max(1);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_jobs).build().map_err(|e| e.to_string())?;
    let t0 = std::time::Instant::now();
    let ml = monlib::MonLib::load(&root)?;
    let header = description_header();
    let failed = std::sync::atomic::AtomicUsize::new(0);
    pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter().with_max_len(1).for_each(|job| {
            if let Err(e) = reduce3::par::run_sequential(|| process_one(job, &ml, &cli.params, &header)) {
                failed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                eprintln!("reduce3: {}: error: {}", job.input.display(), e);
            }
        })
    });
    let failed = failed.into_inner();
    let secs = t0.elapsed().as_secs_f64();
    say(&format!(
        "{} models, {} failed, {:.1} s ({:.1} models/s, {} at once)",
        jobs.len(),
        failed,
        secs,
        jobs.len() as f64 / secs.max(1e-9),
        n_jobs
    ));
    if failed > 0 {
        return Err(format!("{} of {} models failed", failed, jobs.len()));
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    #[cfg(feature = "refcheck")]
    if args.len() > 2 && args[1] == "refcheck" {
        let fixed = args.iter().any(|a| a == "--fixed");
        for p in &args[2..] {
            if p.starts_with("--") {
                continue;
            }
            refcheck::run(p, !fixed);
        }
        return ExitCode::SUCCESS;
    }
    #[cfg(feature = "refcheck")]
    if args.len() > 4 && args[1] == "hcheck" {
        let fixed = args.iter().any(|a| a == "--fixed");
        refcheck::hcheck(&args[2], &args[3], &args[4], !fixed);
        return ExitCode::SUCCESS;
    }
    #[cfg(feature = "refcheck")]
    if args.len() > 2 && args[1] == "typecheck" {
        refcheck::typecheck(&args[2]);
        return ExitCode::SUCCESS;
    }
    #[cfg(feature = "refcheck")]
    if args.len() > 3 && args[1] == "ccdcheck" {
        refcheck::ccdcheck(&args[2], &args[3]);
        return ExitCode::SUCCESS;
    }
    #[cfg(feature = "refcheck")]
    if args.len() > 4 && args[1] == "wcheck" {
        refcheck::wcheck(&args[2], &args[3], &args[4]);
        return ExitCode::SUCCESS;
    }
    if args.len() < 2 {
        eprint!("{}", USAGE);
        return ExitCode::from(2);
    }
    match run_cli(&args[1..]) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.is_empty() => {
            print!("{}", USAGE);
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("reduce3: error: {}", e);
            ExitCode::FAILURE
        }
    }
}
