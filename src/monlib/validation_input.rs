use super::*;

impl MonLib {
    /// Build a private interpretation overlay without changing shared caches.
    pub fn with_validation_cif(&self, inputs: &[(&str, &Path)]) -> Result<Self, String> {
        let mut result = Self {
            root: self.root.clone(),
            comp_synonyms: self.comp_synonyms.clone(),
            atom_synonyms: self.atom_synonyms.clone(),
            links: self.links.clone(),
            link_index: self.link_index.clone(),
            link_list_order: self.link_list_order.clone(),
            mods: self.mods.clone(),
            ener: self.ener.clone(),
            comp_cache: Arc::clone(&self.comp_cache),
            ccd_cache: Arc::clone(&self.ccd_cache),
            ccd_comp_cache: Arc::clone(&self.ccd_comp_cache),
            variant_cache: Arc::clone(&self.variant_cache),
            sites_cache: Arc::clone(&self.sites_cache),
            user_comps: self.user_comps.clone(),
        };
        for &(text, source) in inputs {
            result.validation_dictionary(text, source)?;
        }
        Ok(result)
    }

    fn validation_dictionary(&mut self, text: &str, source: &Path) -> Result<(), String> {
        let doc = cif::parse(text);
        validate_metadata(&doc)?;
        let mut definitions = self.validation_components(&doc, source)?;
        for block in &doc.blocks {
            definitions += self.validation_energy(block)?;
            if let Some(id) = block.name.strip_prefix("link_").filter(|id| *id != "list") {
                let mut link = self.link(id).cloned().unwrap_or_else(|| ChemLink {
                    id: id.into(),
                    ..Default::default()
                });
                update_link_metadata(&doc, id, &mut link);
                parse_link_block(block, &mut link);
                let index = if let Some(&index) = self.link_index.get(id) {
                    index
                } else {
                    let index = self.links.len();
                    self.links.push(ChemLink::default());
                    self.link_index.insert(id.into(), index);
                    index
                };
                self.links[index] = link;
                if !self.link_list_order.contains(&index) {
                    self.link_list_order.push(index);
                }
                definitions += 1;
            }
            if let Some(id) = block.name.strip_prefix("mod_").filter(|id| *id != "list") {
                let mut modification = ChemMod {
                    id: id.into(),
                    ..Default::default()
                };
                parse_mod_block(block, &mut modification);
                self.mods.insert(id.into(), modification);
                definitions += 1;
            }
        }
        if definitions == 0 {
            return Err(
                "The dictionary contains no component, link or modification definitions".into(),
            );
        }
        Ok(())
    }
    fn validation_components(
        &mut self,
        doc: &cif::Document,
        source: &Path,
    ) -> Result<usize, String> {
        let comps = parse_comp_doc(doc, source);
        let definitions = comps.len();
        for comp in comps {
            if comp.atoms.is_empty() {
                return Err(format!("Component {} has no atom definitions", comp.id));
            }
            let mut names = FxHashSet::default();
            if comp
                .atoms
                .iter()
                .any(|a| a.id.is_empty() || !names.insert(a.id.clone()))
            {
                return Err(format!(
                    "Component {} has ambiguous atom definitions",
                    comp.id
                ));
            }
            self.user_comps
                .insert(comp.id.to_ascii_uppercase(), Arc::new(comp));
        }

        Ok(definitions)
    }

    fn validation_energy(&mut self, block: &cif::Block) -> Result<usize, String> {
        let mut definitions = 0;
        if let Some(table) = block.category("_lib_atom") {
            for row in 0..table.nrows() {
                let text = |tag| table.get_tag(row, tag).map(s).unwrap_or_default();
                let number = |tag| table.get_tag(row, tag).and_then(cif::parse_f64);
                let id = text("type");
                if id.is_empty() {
                    return Err("An energy atom definition lacks a type".into());
                }
                self.ener.insert(
                    id.clone(),
                    EnerAtom {
                        type_: id,
                        hb_type: text("hb_type"),
                        element: text("element").to_ascii_uppercase(),
                        vdw_radius: number("vdw_radius"),
                        vdwh_radius: number("vdwh_radius"),
                        ion_radius: number("ion_radius"),
                        vdw_radius_neutron: number("vdw_radius_neutron"),
                    },
                );
                definitions += 1;
            }
        }

        Ok(definitions)
    }
}

fn validate_metadata(doc: &cif::Document) -> Result<(), String> {
    for block in &doc.blocks {
        for table in &block.categories {
            if table.name.starts_with("_chem_")
                && ![
                    "_chem_comp",
                    "_chem_comp_atom",
                    "_chem_comp_bond",
                    "_chem_comp_angle",
                    "_chem_comp_chir",
                    "_chem_comp_plane_atom",
                    "_chem_comp_tor",
                    "_chem_comp_tree",
                    "_chem_link",
                    "_chem_link_bond",
                    "_chem_link_angle",
                    "_chem_link_chir",
                    "_chem_link_plane",
                    "_chem_link_tor",
                    "_chem_mod",
                    "_chem_mod_atom",
                    "_chem_mod_bond",
                    "_chem_mod_angle",
                    "_chem_mod_chir",
                    "_chem_mod_plane_atom",
                    "_chem_mod_tor",
                ]
                .contains(&table.name.as_str())
            {
                return Err(format!(
                    "Unsupported chemical dictionary category {}",
                    table.name
                ));
            }
            if table.name.starts_with("_lib_") && table.name != "_lib_atom" {
                return Err(format!(
                    "Unsupported energy dictionary category {}",
                    table.name
                ));
            }
            for row in 0..table.nrows() {
                for (column, tag) in table.tags.iter().enumerate() {
                    let value = table.get(row, column);
                    if table.name.starts_with("_chem_link_")
                        && tag.starts_with("atom_")
                        && tag.ends_with("comp_id")
                        && !matches!(value, "1" | "2")
                    {
                        return Err(format!(
                            "{} row {row} has invalid component side {value}",
                            table.name
                        ));
                    }
                    if table.name.starts_with("_chem_mod_")
                        && tag == "function"
                        && !matches!(value, "add" | "change" | "delete")
                    {
                        return Err(format!(
                            "{} row {row} has invalid function {value}",
                            table.name
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

fn update_link_metadata(doc: &cif::Document, id: &str, link: &mut ChemLink) {
    for block in &doc.blocks {
        let Some(table) = block.category("_chem_link") else {
            continue;
        };
        for row in 0..table.nrows() {
            if table.get_tag(row, "id").unwrap_or_default() != id {
                continue;
            }
            let get = |name| table.get_tag(row, name).map(s).unwrap_or_default();
            link.comp_id_1 = get("comp_id_1");
            link.mod_id_1 = get("mod_id_1");
            link.group_comp_1 = get("group_comp_1");
            link.comp_id_2 = get("comp_id_2");
            link.mod_id_2 = get("mod_id_2");
            link.group_comp_2 = get("group_comp_2");
            link.name = get("name");
            link.listed = true;
        }
    }
}
