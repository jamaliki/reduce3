//! Replay captured, typed Probe2 inputs without chemical interpretation.
//! Run with `cargo run --features refcheck --example probe_fixture -- FILE`.

#[cfg(not(feature = "refcheck"))]
fn main() {
    eprintln!("This fixture tool requires the refcheck feature");
}

#[cfg(feature = "refcheck")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use reduce3::probe::AtomInfo;
    use reduce3::world::{AtomLabels, World};
    use serde_json::{Value, json};
    use std::collections::HashMap;
    use std::sync::Arc;

    let filename = std::env::args()
        .nth(1)
        .ok_or("Provide a reference capture JSON")?;
    let input: Value = serde_json::from_slice(&std::fs::read(filename)?)?;
    let captured = &input["probe2"]["initial_atom_info_calls"][0];
    let atoms = captured["atoms"]
        .as_array()
        .ok_or("Missing initial Probe2 atoms")?;
    let infos = captured["extra_atom_info"]
        .as_array()
        .ok_or("Missing initial Probe2 atom info")?;
    if atoms.len() != infos.len() {
        return Err("Incomplete atom typing capture".into());
    }
    let text = |v: &Value, key: &str| v[key].as_str().unwrap_or("").to_string();
    let number = |v: &Value, key: &str| v[key].as_f64().unwrap_or(0.0);
    let boolean = |v: &Value, key: &str| v[key].as_bool().unwrap_or(false);
    let mut models = HashMap::new();
    let mut groups = HashMap::new();
    let mut world = World {
        labels: Vec::new(),
        pos: Vec::new(),
        occ: Vec::new(),
        b: Vec::new(),
        info: Vec::new(),
        bonded: vec![Vec::new(); atoms.len()],
        riding: vec![None; atoms.len()],
        n_real: atoms.len(),
    };
    for (atom, info) in atoms.iter().zip(infos) {
        let model = text(atom, "model");
        let chain = text(atom, "chain");
        let residue = text(atom, "resseq");
        let component = text(atom, "component");
        let icode = text(atom, "icode");
        let altloc = text(atom, "altloc");
        let name = text(atom, "name");
        let element = text(atom, "element");
        let group_key = format!("{model}|{chain}|{residue}|{icode}|{component}|{altloc}");
        let model_index = models.len() as u32;
        let model_index = *models.entry(model.clone()).or_insert(model_index);
        let group_index = groups.len() as u32;
        let group_index = *groups.entry(group_key).or_insert(group_index);
        world.labels.push(AtomLabels {
            name: name.clone(),
            raw_name: name,
            element,
            resname: Arc::from(component.as_str()),
            resname_raw: Arc::from(component.as_str()),
            chain: Arc::from(chain),
            resseq: residue.trim().parse()?,
            resseq_raw: Arc::from(residue),
            icode: Arc::from(icode),
            altloc: Arc::from(altloc),
            model_index,
            model_id: Arc::from(model),
            ag: group_index,
            rg: group_index,
            hetero: boolean(atom, "hetero"),
            is_water: reduce3::resclass::get_class(&component)
                == reduce3::resclass::ResClass::CommonWater,
            class: None,
        });
        let position = atom["xyz"].as_array().ok_or("Missing coordinates")?;
        world.pos.push(reduce3::geom::Vec3::from_array([
            position[0].as_f64().ok_or("Bad coordinate")?,
            position[1].as_f64().ok_or("Bad coordinate")?,
            position[2].as_f64().ok_or("Bad coordinate")?,
        ]));
        world.occ.push(number(atom, "occupancy"));
        world.b.push(number(atom, "b_iso"));
        world.info.push(AtomInfo {
            vdw_radius: number(info, "vdw_radius"),
            is_acceptor: boolean(info, "is_acceptor"),
            is_donor: boolean(info, "is_donor"),
            is_dummy_hydrogen: boolean(info, "is_dummy_hydrogen"),
            is_ion: boolean(info, "is_ion"),
            charge: info["charge"].as_i64().ok_or("Missing charge")? as i8,
            alt: text(info, "alt").bytes().next().unwrap_or(0),
        });
    }
    for bond in captured["bonds"]
        .as_array()
        .ok_or("Missing covalent bonds")?
    {
        let values = bond.as_array().ok_or("Bad bond pair")?;
        let i = values[0].as_u64().ok_or("Bad bond index")? as usize;
        let j = values[1].as_u64().ok_or("Bad bond index")? as usize;
        if i >= atoms.len() || j >= atoms.len() {
            return Err("Bond index is out of range".into());
        }
        world.bonded[i].push(j as u32);
        world.bonded[j].push(i as u32);
    }
    let mut options =
        reduce3::clashcheck::Options::for_clashscore(&world, boolean(&input, "nuclear"));
    options.report_vdws = boolean(&input, "report_vdws");
    let report = reduce3::clashcheck::run(world, &options)?;
    let key = |i: u32| {
        let a = &report.atoms[i as usize];
        json!([
            a.model_id.as_ref(),
            a.chain.as_ref(),
            a.resseq,
            a.icode.trim(),
            a.resname.as_ref(),
            a.name.trim(),
            a.altloc.trim()
        ])
    };
    let pairs: Vec<_> = report
        .clashes
        .iter()
        .map(|clash| {
            json!({
        "first":key(clash.first),"second":key(clash.second),"overlap":clash.overlap,
        "location":clash.location.as_array(),"max_b_factor":clash.max_b_factor })
        })
        .collect();
    let contacts: Vec<_> = report
        .contacts
        .iter()
        .map(|contact| {
            use reduce3::probe::InteractionType::*;
            let kind = match contact.interaction_type {
                WideContact => "wc",
                CloseContact => "cc",
                WeakHydrogenBond => "wh",
                SmallOverlap => "so",
                Bump | BadBump => "bo",
                StandardHydrogenBond => "hb",
                Invalid => "invalid",
            };
            json!({"source":key(contact.source),"target":key(contact.target),"type":kind,
            "dot_count":contact.dot_count,"min_dot_gap":contact.min_dot_gap,
            "score_over_density":contact.score_over_density})
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({"pairs":pairs,
        "clashscore":report.clashscore,"n_selected":report.n_selected_nonwater_atoms,
        "phantom_count":report.phantom_count,"contact_count":report.contacts.len(),
        "contacts":contacts}))?
    );
    Ok(())
}
