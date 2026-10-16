//! Preserve source hydrogen names at modified peptide termini.

use super::*;

fn parents<'a>(comp: &'a Comp, hydrogen: &str) -> Vec<&'a str> {
    comp.bonds.iter().filter_map(|bond| {
        let other = if bond.a1 == hydrogen { &bond.a2 }
            else if bond.a2 == hydrogen { &bond.a1 } else { return None };
        (!comp.is_h(other)).then_some(other.as_str())
    }).collect()
}

fn collision(comp: &Comp) -> bool {
    comp.test_for_peptide() && ["H1", "H2", "H3", "HN", "HN1", "HN2"].iter()
        .any(|name| comp.is_h(name) && parents(comp, name).iter().any(|&parent| parent != "N"))
}

fn names(comp: &Comp) -> [String; 3] {
    let mut serial = 1;
    std::array::from_fn(|_| loop {
        let name = format!("HNT{serial}");
        serial += 1;
        if !comp.has_atom(&name) { break name; }
    })
}

/// Return the reserved terminal names when source H names have other parents.
pub fn hydrogen_parents(comp: &Comp) -> Vec<(String, String)> {
    if !collision(comp) { return Vec::new(); }
    names(comp).into_iter().map(|name| (name, "N".into())).collect()
}

fn dictionary(ml: &MonLib, base: &Comp, names: &[String; 3]) -> Option<Comp> {
    let glycine = ml.comp("GLY")?;
    let charged = ml.apply_mod(&glycine, "NH3")?;
    let mut result = base.clone();
    let rename = |name: &str| match name {
        "H1" => names[0].clone(), "H2" => names[1].clone(), "H3" => names[2].clone(),
        _ => name.to_owned(),
    };
    let terminal = |name: &str| ["H1", "H2", "H3"].contains(&name);
    for old in ["H1", "H2", "H3"] {
        let mut atom = charged.atom(old)?.clone();
        if parents(&charged, old) != ["N"] { return None; }
        atom.id = rename(old);
        result.atoms.push(atom);
    }
    for source in &charged.bonds {
        if !terminal(&source.a1) && !terminal(&source.a2) { continue; }
        let mut bond = source.clone();
        bond.a1 = rename(&bond.a1);
        bond.a2 = rename(&bond.a2);
        result.bonds.push(bond);
    }
    for source in &charged.angles {
        if ![&source.a1, &source.a2, &source.a3].iter().any(|name| terminal(name)) { continue; }
        let mut angle = source.clone();
        angle.a1 = rename(&angle.a1);
        angle.a2 = rename(&angle.a2);
        angle.a3 = rename(&angle.a3);
        result.angles.push(angle);
    }
    for source in &charged.tors {
        if !source.a.iter().any(|name| terminal(name)) { continue; }
        let mut torsion = source.clone();
        torsion.a = torsion.a.map(|name| rename(&name));
        result.tors.push(torsion);
    }
    Some(result)
}

/// Keep source H names when a peptide dictionary uses terminal aliases elsewhere.
pub(super) fn place(rg: &mut ResidueGroup, members: &[(usize, usize)], np: (usize, usize),
    cap: (usize, usize), cp: (usize, usize), ml: &MonLib,
    auto_comps: &mut FxHashMap<String, Arc<Comp>>) -> bool {
    let group = &rg.atom_groups[np.0];
    let heavy: Vec<_> = group.atoms.iter().filter(|atom| !atom.is_hydrogen()).map(|atom| atom.name.clone()).collect();
    let Some(base) = hydrogen_dictionary(ml, &group.resname, &heavy, false) else { return false };
    if !collision(&base) { return false; }
    let key = group.resname.trim().to_ascii_uppercase();
    let names = names(&base);
    let comp = if let Some(comp) = auto_comps.get(&key).filter(|comp| names.iter().all(|name| comp.has_atom(name))) {
        comp.clone()
    } else {
        let Some(comp) = dictionary(ml, &base, &names) else { return false };
        Arc::new(comp)
    };
    let nitrogen = rg.atom_groups[np.0].atoms[np.1].clone();
    let ca = rg.atom_groups[cap.0].atoms[cap.1].xyz;
    let carbon = rg.atom_groups[cp.0].atoms[cp.1].xyz;
    let positions = construct_xyz(nitrogen.xyz, 1.0, ca, 109.5, carbon, 120.0);
    let groups: FxHashSet<_> = members.iter().map(|&(group, _)| group).collect();
    for group in groups {
        rg.atom_groups[group].atoms.retain(|atom| !atom.is_hydrogen()
            || parents(&base, atom.name.trim()) != ["N"]);
    }
    let group = &mut rg.atom_groups[np.0];
    for (name, position) in names.iter().zip(positions) {
        let mut hydrogen = Atom::new(&pad_name(name), "H", position);
        hydrogen.occ = nitrogen.occ;
        hydrogen.b = nitrogen.b;
        hydrogen.hetero = nitrogen.hetero;
        group.atoms.push(hydrogen);
    }
    auto_comps.insert(key, comp);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn private_names_keep_internal_source_atoms_and_restraints() {
        let Some(root) = MonLib::locate(None) else { return };
        let ml = MonLib::load(&root).unwrap();
        let source = ml.comp("0G7").unwrap();
        assert!(collision(&source));
        let names = ["HNT1", "HNT2", "HNT3"].map(String::from);
        let working = dictionary(&ml, &source, &names).unwrap();
        assert_eq!(parents(&working, "HN2"), ["N2"]);
        assert_eq!(hydrogen_parents(&source), ["HNT1", "HNT2", "HNT3"].map(|name| (name.into(), "N".into())));
        assert_eq!(parents(&working, "H3"), ["C3"]);
        assert_eq!(working.atoms[..source.atoms.len()].iter().map(|a| (&a.id, &a.type_energy)).collect::<Vec<_>>(),
            source.atoms.iter().map(|a| (&a.id, &a.type_energy)).collect::<Vec<_>>());
        assert_eq!(working.bonds[..source.bonds.len()].iter().map(|b| (&b.a1, &b.a2, b.value_dist)).collect::<Vec<_>>(),
            source.bonds.iter().map(|b| (&b.a1, &b.a2, b.value_dist)).collect::<Vec<_>>());
        assert_eq!(working.angles[..source.angles.len()].iter().map(|a| (&a.a1, &a.a2, &a.a3, a.value)).collect::<Vec<_>>(),
            source.angles.iter().map(|a| (&a.a1, &a.a2, &a.a3, a.value)).collect::<Vec<_>>());
        for name in names { assert_eq!(parents(&working, &name), ["N"]); }
        assert!(!collision(&ml.comp("GLY").unwrap()));
        assert!(hydrogen_parents(&ml.comp("GLY").unwrap()).is_empty());
        let mut non_peptide = source.as_ref().clone();
        non_peptide.atoms.retain(|a| a.id != "CA");
        assert!(!collision(&non_peptide));
    }
}
