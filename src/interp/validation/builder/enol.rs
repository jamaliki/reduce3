use super::*;

pub(super) fn prepare(
    graph: &mut ValidationGraph,
    residues: &mut [ResInterp],
    sources: &[ConfResidue],
    flat: &FlatAtoms,
) {
    for k in 0..residues.len() {
        let Some(hydrogen) = sources[k]
            .atoms
            .iter()
            .copied()
            .find(|&i| flat.name[i as usize].trim() == "HNO")
        else {
            continue;
        };
        if !residues[k].is_peptide {
            continue;
        }
        let Some(base) = &residues[k].comp else {
            continue;
        };
        let mut comp = (**base).clone();
        if comp.atom("HNO").is_none() {
            add_hno(&mut comp);
        }
        if let Some(atom_type) = comp.atom("HNO").and_then(|a| a.type_energy.clone()) {
            graph.topology.etype[hydrogen as usize] = EType::Typed(atom_type.into());
        }
        residues[k].comp = Some(Arc::new(comp));
        residues[k].id_to_atom.insert("HNO".into(), hydrogen);
        if k + 1 < residues.len() && sources[k + 1].link_to_previous {
            remove_next_hydrogen(graph, &mut residues[k + 1]);
        }
    }
}

fn add_hno(comp: &mut Comp) {
    use crate::monlib::{CompAngle, CompAtom, CompBond};
    comp.atoms.push(CompAtom {
        id: "HNO".into(),
        type_symbol: "H".into(),
        type_energy: Some("HOH1".into()),
    });
    comp.bonds.push(CompBond {
        a1: "HNO".into(),
        a2: "O".into(),
        value_dist: Some(0.97),
        value_dist_neutron: Some(0.97),
        esd: Some(0.02),
        ..Default::default()
    });
    for bond in &mut comp.bonds {
        if (bond.a1 == "C" && bond.a2 == "O") || (bond.a1 == "O" && bond.a2 == "C") {
            bond.value_dist = Some(1.43);
            bond.value_dist_neutron = Some(1.43);
            bond.esd = Some(0.02);
        }
    }
    comp.angles.push(CompAngle {
        a1: "HNO".into(),
        a2: "O".into(),
        a3: "C".into(),
        value: Some(120.),
        esd: Some(3.),
    });
    comp.applied_mods.push("reference-enol-peptide".into());
}

fn remove_next_hydrogen(graph: &mut ValidationGraph, residue: &mut ResInterp) {
    if lookup(residue, "H").is_some() {
        graph.unsupported(ALL, "Enol peptide has a hydrogen on the following nitrogen");
    } else if let Some(base) = &residue.comp {
        let mut comp = (**base).clone();
        comp.atoms.retain(|a| a.id != "H");
        comp.bonds.retain(|b| b.a1 != "H" && b.a2 != "H");
        comp.angles
            .retain(|a| a.a1 != "H" && a.a2 != "H" && a.a3 != "H");
        comp.chirs
            .retain(|c| c.centre != "H" && c.a.iter().all(|a| a != "H"));
        comp.planes.retain(|p| p.atom != "H");
        residue.comp = Some(Arc::new(comp));
    }
}

pub(super) fn links(graph: &mut ValidationGraph, residues: &[ResInterp], sources: &[ConfResidue]) {
    for k in 1..residues.len() {
        if !sources[k].link_to_previous
            || lookup(&residues[k - 1], "HNO").is_none()
            || sources[k - 1..=k]
                .iter()
                .any(|r| resclass::get_class(&r.resname) != resclass::ResClass::CommonAminoAcid)
        {
            continue;
        }
        let (Some(c), Some(n), Some(ca0), Some(ca1)) = (
            lookup(&residues[k - 1], "C"),
            lookup(&residues[k], "N"),
            lookup(&residues[k - 1], "CA"),
            lookup(&residues[k], "CA"),
        ) else {
            continue;
        };
        if dihedral_model(&graph.sites, [ca0, c, n, ca1]).is_some_and(|omega| omega.abs() < 45.) {
            continue;
        }
        let key = [c.min(n), c.max(n)];
        if let Some(&index) = graph.bond_keys.get(&key) {
            graph.bonds[index].ideal = 1.27;
            graph
                .topology
                .upsert_bond(c, n, 1.27, graph.bonds[index].origin);
        }
    }
}
