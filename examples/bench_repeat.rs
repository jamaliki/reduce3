//! Run the pipeline on one model many times in one process (warm caches).
//! Usage: bench_repeat CHEM_DATA MODEL N
use reduce3::{monlib::MonLib, pipeline, Params};
use std::time::Instant;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let ml = MonLib::load(std::path::Path::new(&a[1])).unwrap();
    let text = std::fs::read_to_string(&a[2]).unwrap();
    let n: usize = a[3].parse().unwrap();
    let mut params = Params::default();
    params.opt.add_flip_movers = true;
    let mut times = Vec::new();
    let (mut tr, mut tp, mut tw, mut add, mut opt) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let field = |d: &str, key: &str| -> f64 {
        d.lines().find(|l| l.contains(key)).and_then(|l| l.split('=').nth(1)).and_then(|v| v.trim().trim_end_matches("sec").trim().parse::<f64>().ok()).unwrap_or(0.0) * 1e3
    };
    let is_cif = a[2].ends_with(".cif");
    for _ in 0..n {
        let t = Instant::now();
        // as the command line does: mmCIF is parsed once and written back into its block
        let doc = is_cif.then(|| reduce3::cif::parse(&text));
        let block = doc.as_ref().map(|d| d.blocks.iter().find(|b| b.category("_atom_site").is_some()).unwrap());
        let st = match block { Some(b) => reduce3::mmcif::structure_from_cif(b).unwrap(), None => reduce3::pdbio::read_pdb(&text) };
        let t1 = Instant::now();
        let out = reduce3::par::run_sequential(|| pipeline::run(st, &ml, &params)).unwrap();
        let t2 = Instant::now();
        let s = match block {
            Some(b) => {
                let mut w = reduce3::cifsource::CifText::with_capacity(text.len() * 2);
                reduce3::mmcif::write_cif_preserving(&out.structure, b, b.name, &mut w).unwrap();
                w.out
            }
            None => reduce3::pdbio::write_pdb_preserving(&out.structure),
        };
        std::hint::black_box(s);
        tr.push((t1 - t).as_secs_f64() * 1e3);
        tp.push((t2 - t1).as_secs_f64() * 1e3);
        tw.push(t2.elapsed().as_secs_f64() * 1e3);
        add.push(field(&out.description, "Time to Add Hydrogen"));
        opt.push(field(&out.description, "Time to Optimize"));
        times.push(t.elapsed().as_secs_f64() * 1e3);
    }
    let med = |v: &mut Vec<f64>| { v.sort_by(|x, y| x.partial_cmp(y).unwrap()); v[v.len() / 2] };
    eprintln!("median read {:.2}  pipeline {:.2} (add H {:.1}, optimize {:.1})  write {:.2} ms", med(&mut tr), med(&mut tp), med(&mut add), med(&mut opt), med(&mut tw));
    times.sort_by(|x, y| x.partial_cmp(y).unwrap());
    println!("first {:.1} ms, median {:.1} ms, min {:.1} ms", times.iter().cloned().fold(0.0, f64::max), times[n / 2], times[0]);
}
