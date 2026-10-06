//! Reduce3: a fast Rust reimplementation of cctbx Reduce2 (hydrogen addition
//! and optimization of rotatable and flippable groups).

use reduce3::{hplace, mmcif, model, monlib, pdbio, pipeline};
#[cfg(feature = "refcheck")]
use reduce3::refcheck;
use pipeline::{Approach, Params};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

const USAGE: &str = "\
Reduce3: add hydrogens to a macromolecular model and optimize the rotatable
and flippable groups (a Rust reimplementation of cctbx Reduce2).

Usage: reduce3 [options] model.pdb [name=value ...]

Options:
  --compat              reproduce Reduce2 exactly, including its known bugs
  --chem-data DIR       chem_data directory (monomer library, CCD). Default:
                        $REDUCE3_CHEM_DATA, $CHEM_DATA, or the active conda env
  -o, --output FILE     output model (same as output.filename=FILE)
  --threads N           worker threads (default: all cores)
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
    input: PathBuf,
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
        input: PathBuf::new(),
        output: None,
        description: None,
        write_files: true,
        chem_data: None,
        threads: None,
        quiet: false,
        params: Params::default(),
    };
    let mut input: Option<PathBuf> = None;
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
            _ => {
                if input.is_some() {
                    return Err(format!("more than one input model given ({})", a));
                }
                input = Some(PathBuf::from(a));
            }
        }
        i += 1;
    }
    cli.input = input.ok_or_else(|| "no input model given".to_string())?;
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

fn default_output(input: &Path, flips: bool) -> PathBuf {
    let stem = input.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "model".into());
    let ext = input.extension().map(|s| format!(".{}", s.to_string_lossy())).unwrap_or_default();
    PathBuf::from(format!("{}{}{}", stem, if flips { "FH" } else { "H" }, ext))
}

fn run_cli(args: &[String]) -> Result<(), String> {
    let cli = parse_args(args)?;
    if let Some(n) = cli.threads {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n.max(1)).build_global();
    }
    let say = |m: &str| {
        if !cli.quiet {
            eprintln!("{}", m);
        }
    };
    let root = monlib::MonLib::locate(cli.chem_data.as_deref()).ok_or_else(|| {
        "could not find chem_data (the monomer library). Pass --chem-data DIR or set REDUCE3_CHEM_DATA.".to_string()
    })?;
    let t0 = std::time::Instant::now();
    let text = std::fs::read_to_string(&cli.input).map_err(|e| format!("cannot read {}: {}", cli.input.display(), e))?;
    let lower = cli.input.to_string_lossy().to_ascii_lowercase();
    let is_cif = lower.ends_with(".cif") || lower.ends_with(".mmcif") || text.trim_start().starts_with("data_");
    model::mem_checkpoint("start");
    let st = if is_cif { mmcif::read_mmcif(&text)? } else { pdbio::read_pdb(&text) };
    drop(text);
    model::mem_checkpoint("read");
    if st.atoms_size() == 0 {
        return Err(format!("no atoms found in {}", cli.input.display()));
    }
    let ml = monlib::MonLib::load(&root)?;
    say(&format!("Read {} ({} atoms)", cli.input.display(), st.atoms_size()));
    let out = pipeline::run(st, &ml, &cli.params)?;
    model::mem_checkpoint("pipeline done");
    let output = cli.output.clone().unwrap_or_else(|| default_output(&cli.input, cli.params.opt.add_flip_movers));
    let out_lower = output.to_string_lossy().to_ascii_lowercase();
    let write_pdb_format = out_lower.ends_with(".pdb") || out_lower.ends_with(".ent");
    let description = cli.description.clone().unwrap_or_else(|| {
        let s = output.to_string_lossy().replace(".pdb", ".txt").replace(".cif", ".txt");
        PathBuf::from(s)
    });
    if cli.write_files {
        let mut desc = format!("reduce3 v.{}, run {}\n", env!("CARGO_PKG_VERSION"), utc_timestamp());
        for a in std::env::args() {
            desc.push(' ');
            desc.push_str(&a);
        }
        desc.push('\n');
        desc.push_str(&out.description);
        let model_text = if write_pdb_format && cli.params.compat {
            pdbio::write_pdb(&out.structure, true)
        } else if write_pdb_format {
            // keep SSBOND, LINK and CONECT (fixed mode)
            pdbio::write_pdb_preserving(&out.structure)
        } else if is_cif && !cli.params.compat {
            // keep every category of the input block (fixed mode)
            let text = std::fs::read_to_string(&cli.input).map_err(|e| format!("cannot read {}: {}", cli.input.display(), e))?;
            let doc = reduce3::cif::parse(&text);
            let block = doc
                .blocks
                .iter()
                .find(|b| b.category("_atom_site").is_some())
                .ok_or_else(|| "no _atom_site loop found in the mmCIF file".to_string())?;
            let mut w = reduce3::cifsource::CifText::with_capacity(text.len() + out.structure.atoms_size() * 100);
            mmcif::write_cif_preserving(&out.structure, block, block.name, &mut w)?;
            w.out
        } else {
            mmcif::write_mmcif(&out.structure)
        };
        std::fs::write(&output, model_text).map_err(|e| format!("cannot write {}: {}", output.display(), e))?;
        std::fs::write(&description, desc).map_err(|e| format!("cannot write {}: {}", description.display(), e))?;
        say(&format!("Wrote {} and {} ({:.3} s)", output.display(), description.display(), t0.elapsed().as_secs_f64()));
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
