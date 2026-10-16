use reduce3::{hplace, mmcif, MonLib};

fn structure() -> reduce3::model::Structure {
    mmcif::read_mmcif(concat!(
        "data_afdb_ser200\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n",
        "_atom_site.label_alt_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n",
        "_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "N N . SER A 200 -9.168000221252441 3.950000047683716 23.77199935913086 1 0\n",
        "C CA . SER A 200 -8.21399974822998 2.8529999256134033 24.030000686645508 1 0\n",
        "C C . SER A 200 -7.635000228881836 2.2309999465942383 22.75200080871582 1 0\n",
        "O O . SER A 200 -8.364999771118164 2.109999895095825 21.764999389648438 1 0\n",
        "C CB . SER A 200 -8.914999961853027 1.7419999837875366 24.819000244140625 1 0\n",
        "O OG . SER A 200 -9.987000465393066 1.2289999723434448 24.047000885009766 1 0\n",
    )).unwrap()
}

#[test]
fn ser200_places_ha_in_fixed_mode_and_keeps_compat_and_heavy_sites() {
    let library = MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap();
    for compat in [false, true] {
        let mut st = structure();
        let heavy: Vec<_> = st.atom_paths().into_iter().map(|p| (st.atom(p).name.clone(), st.atom(p).xyz)).collect();
        let result = hplace::place_hydrogens(&mut st, &library, &hplace::HPlaceParams {
            neutron: false, n_terminal_charge: hplace::NTermCharge::NoCharge,
            exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat, cell: None,
        });
        assert_eq!(result.site_labels_no_para.iter().any(|label| label.trim_start().starts_with("HA ")), compat);
        for (name, xyz) in heavy {
            let atom = st.atom_paths().into_iter().map(|p| st.atom(p)).find(|atom| atom.name == name).unwrap();
            assert_eq!(atom.xyz, xyz);
        }
        if !compat {
            let site = |name: &str| st.atom_paths().into_iter().map(|p| st.atom(p))
                .find(|atom| atom.name.trim() == name).unwrap().xyz;
            let (ca, ha, n, c, cb) = (site("CA"), site("HA"), site("N"), site("C"), site("CB"));
            assert!((ha.dist(ca) - 0.97).abs() < 0.03);
            let normal = (n - ca).cross(c - ca);
            assert!(normal.dot(ha - ca) * normal.dot(cb - ca) < 0.0);
        }
    }
}
