//! Golden outputs of fixed mode, the default, for the structures in
//! tests/golden (see its README.md): any change to a hydrogen's coordinates
//! as written (to 1e-3 Å), to a heavy atom's, or to a flip, rotor or His
//! decision in the report fails `fixed_mode_matches_golden_outputs`.
//!
//! The outputs depend on the platform (fused multiply-adds on aarch64 only,
//! the platform's libm), so each platform has its own file,
//! `fixed-<arch>-<os>.tsv`. They also depend on chem_data, so the test first
//! checks the files the set reads against the pinned version
//! (`chem_data.tsv`; `tools/golden/fetch_chem_data.py` fetches them).
//!
//! - `REDUCE3_GOLDEN_BLESS=1` writes this platform's file instead of
//!   comparing (after a deliberate change; review its diff).
//! - `REDUCE3_GOLDEN_DUMP=<dir>` also writes each structure's hydrogens,
//!   heavy atoms (exact and as written) and report to `<dir>/<id>.txt`.
//!
//! `determinism_across_thread_counts` is slow in debug builds, where it is
//! ignored: run `cargo test --release --test golden`.

#![cfg(feature = "gzip")]

#[path = "common/chem_data.rs"]
mod chem_data;

use std::fmt::Write as _;
use std::io::Read;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use reduce3::{mmcif, par, pipeline, MonLib, Output, Params};
use sha2::{Digest, Sha256};

const BLESS: &str = "REDUCE3_GOLDEN_BLESS";
const DUMP: &str = "REDUCE3_GOLDEN_DUMP";

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn golden_file() -> PathBuf {
    golden_dir().join(format!("fixed-{}-{}.tsv", std::env::consts::ARCH, std::env::consts::OS))
}

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes).iter().fold(String::with_capacity(64), |mut s, b| {
        let _ = write!(s, "{b:02x}");
        s
    })
}

fn env_flag(name: &str) -> bool {
    std::env::var_os(name).is_some_and(|v| !v.is_empty() && v != "0")
}

/// Rows of a tab-separated file with a header row, as maps from column name.
fn read_tsv(file: &Path) -> Vec<std::collections::HashMap<String, String>> {
    let text = std::fs::read_to_string(file).unwrap_or_else(|e| panic!("cannot read {}: {e}", file.display()));
    let mut lines = text.lines().filter(|l| !l.starts_with('#') && !l.is_empty());
    let head: Vec<&str> = lines.next().expect("a header row").split('\t').collect();
    lines
        .map(|l| {
            let cells: Vec<&str> = l.split('\t').collect();
            assert_eq!(cells.len(), head.len(), "{}: {l}", file.display());
            head.iter().zip(cells).map(|(h, c)| (h.to_string(), c.to_string())).collect()
        })
        .collect()
}

/// The chem_data root, checked against the pinned version, or None to skip.
fn chem_data() -> Option<MonLib> {
    let root = chem_data::locate()?;
    let pins = golden_dir().join("chem_data.tsv");
    let mut wrong = Vec::new();
    let rows = read_tsv(&pins);
    for row in &rows {
        let file = root.join(&row["path"]);
        let got = std::fs::read(&file).map_or_else(|_| "missing".to_string(), |b| sha256(&b));
        if got != row["sha256"] {
            wrong.push(format!("{}: {got}", row["path"]));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} of the {} chem_data files in {} differ from the pinned version ({}); fetch it with \
         tools/golden/fetch_chem_data.py. First: {:?}",
        wrong.len(),
        rows.len(),
        root.display(),
        pins.display(),
        &wrong[..wrong.len().min(5)]
    );
    Some(MonLib::load(&root).unwrap())
}

struct Input {
    id: String,
    text: String,
}

fn inputs() -> Vec<Input> {
    let dir = golden_dir();
    read_tsv(&dir.join("structures.tsv"))
        .into_par_iter()
        .map(|row| {
            let file = dir.join("inputs").join(format!("{}.cif.gz", row["id"]));
            let gz = std::fs::read(&file).unwrap_or_else(|e| panic!("cannot read {}: {e}", file.display()));
            let mut text = String::new();
            flate2::read::GzDecoder::new(&gz[..]).read_to_string(&mut text).unwrap();
            assert_eq!(sha256(text.as_bytes()), row["input_sha256"], "{} changed", file.display());
            Input { id: row["id"].clone(), text }
        })
        .collect()
}

/// One structure's fixed-mode result: what the golden file records, and the
/// exact output for comparing runs bit for bit.
struct Outcome {
    /// `id`, `status`, `hydrogens`, `hydrogen_sha256`, `movers`,
    /// `decision_sha256`, `heavy_sha256`.
    row: [String; 7],
    /// sha256 of every atom's coordinates, occupancy and B to the bit, and
    /// of the description but its timings.
    exact: String,
    decisions: String,
    dump: String,
}

const COLUMNS: [&str; 7] = ["id", "status", "hydrogens", "hydrogen_sha256", "movers", "decision_sha256", "heavy_sha256"];

fn run(input: &Input, ml: &MonLib) -> Outcome {
    let out = mmcif::read_mmcif(&input.text).and_then(|st| pipeline::run(st, ml, &Params::default()));
    match out {
        Ok(out) => describe(&input.id, &out),
        Err(e) => Outcome {
            row: [input.id.clone(), "error".into(), "0".into(), sha256(e.as_bytes()), "0".into(), "-".into(), "-".into()],
            exact: sha256(e.as_bytes()),
            decisions: String::new(),
            dump: format!("# error\n{e}\n"),
        },
    }
}

fn describe(id: &str, out: &Output) -> Outcome {
    let st = &out.structure;
    let (mut hydrogens, mut heavy, mut exact) = (String::new(), String::new(), String::new());
    let (mut dump_h, mut dump_heavy) = (String::new(), String::new());
    let mut n_h = 0;
    for m in &st.models {
        for c in &m.chains {
            for rg in &c.residue_groups {
                for ag in &rg.atom_groups {
                    for a in &ag.atoms {
                        let label = format!(
                            "{}\t{}\t{}{}\t{}\t{}\t{}",
                            m.id, c.id, rg.resseq.trim(), rg.icode.trim(), ag.altloc.trim(), ag.resname.trim(), a.name.trim()
                        );
                        let [x, y, z] = [a.xyz.x, a.xyz.y, a.xyz.z];
                        // as written: `fastfmt::push_fixed` is `{:.3}`
                        let written = format!("{label}\t{x:.3}\t{y:.3}\t{z:.3}\n");
                        let full = format!("{label}\t{x:.3}\t{y:.3}\t{z:.3}\t{x:?}\t{y:?}\t{z:?}\n");
                        let _ = writeln!(exact, "{label}\t{:x}\t{:x}\t{:x}\t{:x}\t{:x}\t{}", x.to_bits(), y.to_bits(), z.to_bits(),
                            a.occ.to_bits(), a.b.to_bits(), a.element.as_str());
                        if a.is_hydrogen() {
                            n_h += 1;
                            hydrogens += &written;
                            dump_h += &full;
                        } else {
                            heavy += &written;
                            dump_heavy += &full;
                        }
                    }
                }
            }
        }
    }
    // The report's decisions: each Mover and its final pose, without scores.
    let (mut decisions, mut report, mut in_report, mut movers) = (String::new(), String::new(), false, 0);
    for line in out.description.lines() {
        if line.trim_start().starts_with("BEGIN REPORT") {
            in_report = true;
            decisions += line.trim();
            decisions.push('\n');
        }
        if in_report {
            report += line;
            report.push('\n');
            if let Some((mover, pose)) = line.split_once(" pose ") {
                let mover = mover.trim();
                let mover = mover.split(" Initial score:").next().unwrap();
                let mover = mover.split(" final score:").next().unwrap();
                let _ = writeln!(decisions, "{mover}\t{pose}");
                movers += 1;
            }
        }
        if line.trim_start().starts_with("END REPORT") {
            in_report = false;
        }
    }
    // Everything but the timings, for comparing runs.
    for line in out.description.lines().filter(|l| !l.trim_start().starts_with("Time to ")) {
        exact += line;
        exact.push('\n');
    }
    Outcome {
        row: [
            id.to_string(),
            "ok".into(),
            n_h.to_string(),
            sha256(hydrogens.as_bytes()),
            movers.to_string(),
            sha256(decisions.as_bytes()),
            sha256(heavy.as_bytes()),
        ],
        exact: sha256(exact.as_bytes()),
        decisions,
        dump: format!("# hydrogens\n{dump_h}# heavy atoms\n{dump_heavy}# report\n{report}"),
    }
}

/// Every structure, run as the batch mode runs them: several at once, each
/// one sequentially (`determinism_across_thread_counts` checks that this
/// gives what one model run on its own gives).
fn run_all(inputs: &[Input], ml: &MonLib) -> Vec<Outcome> {
    inputs.par_iter().map(|i| par::run_sequential(|| run(i, ml))).collect()
}

#[test]
fn fixed_mode_matches_golden_outputs() {
    let Some(ml) = chem_data() else { return };
    let inputs = inputs();
    let results = run_all(&inputs, &ml);
    if let Some(dir) = std::env::var_os(DUMP) {
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).unwrap();
        for r in &results {
            std::fs::write(dir.join(format!("{}.txt", r.row[0])), &r.dump).unwrap();
        }
    }
    // the set exercises every kind of Mover
    for kind in ["SingleHydrogenRotator", "NH3Rotator", "AromaticMethylRotator", "AmideFlip", "HisFlip"] {
        assert!(results.iter().any(|r| r.decisions.contains(&format!("{kind} at "))), "no {kind} decision in the set");
    }
    let mut table = COLUMNS.join("\t") + "\n";
    for r in &results {
        table += &r.row.join("\t");
        table.push('\n');
    }
    let file = golden_file();
    if env_flag(BLESS) {
        std::fs::write(&file, table).unwrap();
        eprintln!("wrote {}", file.display());
        return;
    }
    let golden = read_tsv(&file);
    assert_eq!(golden.len(), results.len(), "{}: the structures changed; bless it again", file.display());
    let mut changed = Vec::new();
    for (g, r) in golden.iter().zip(&results) {
        assert_eq!(g["id"], r.row[0], "{}: the structures changed; bless it again", file.display());
        let what: Vec<&str> = COLUMNS.iter().zip(&r.row).filter(|(c, v)| g[**c] != **v).map(|(c, _)| *c).collect();
        if !what.is_empty() {
            changed.push(format!("{} ({})", r.row[0], what.join(", ")));
        }
    }
    assert!(
        changed.is_empty(),
        "{} of {} structures differ from {}:\n  {}\nSet {DUMP}=<dir> to see the outputs, or {BLESS}=1 to accept them.",
        changed.len(),
        results.len(),
        file.display(),
        changed.join("\n  ")
    );
}

/// One rayon thread and several give the same output, bit for bit, as does
/// running each model sequentially among others (the golden test's way).
#[test]
#[cfg_attr(debug_assertions, ignore = "slow in debug builds: cargo test --release --test golden")]
fn determinism_across_thread_counts() {
    let Some(ml) = chem_data() else { return };
    let inputs = inputs();
    let many = std::thread::available_parallelism().map_or(4, |n| n.get().max(4));
    let with = |threads: usize| -> Vec<String> {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        pool.install(|| inputs.iter().map(|i| run(i, &ml).exact).collect())
    };
    let one = with(1);
    let several = with(many);
    let batch: Vec<String> = run_all(&inputs, &ml).into_iter().map(|r| r.exact).collect();
    let differ = |other: &[String]| -> Vec<&str> {
        inputs.iter().zip(one.iter().zip(other)).filter(|(_, (a, b))| a != b).map(|(i, _)| i.id.as_str()).collect()
    };
    let (threads, batched) = (differ(&several), differ(&batch));
    assert!(threads.is_empty(), "1 thread and {many} threads differ on {threads:?}");
    assert!(batched.is_empty(), "1 thread and batch mode differ on {batched:?}");
}
