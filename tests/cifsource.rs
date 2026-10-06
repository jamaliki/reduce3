//! The in-memory CIF interface: a source that is not Reduce3's own parser must
//! build the same model, and the sink must receive the block `write_mmcif`
//! writes.

use reduce3::cifsource::{CifCell, CifSink, CifSource, CifTable};
use reduce3::{cif, mmcif};
use std::borrow::Cow;
use std::collections::HashMap;

const MODEL: &str = "\
data_test
_cell.length_a 40.000
_cell.length_b 50.000
_cell.length_c 60.000
_cell.angle_alpha 90.00
_cell.angle_beta 90.00
_cell.angle_gamma 90.00
_symmetry.space_group_name_H-M 'P 21 21 21'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1 N N   . SER A 1 ? 10.000 12.000 14.000 1.00 20.00 ? 1 SER A N   1
ATOM   2 C CA  . SER A 1 ? 11.400 12.100 14.200 1.00 21.50 ? 1 SER A CA  1
ATOM   3 C C   . SER A 1 ? 12.000 13.500 14.100 1.00 22.00 ? 1 SER A C   1
ATOM   4 O O   . SER A 1 ? 11.300 14.500 14.000 1.00 23.00 ? 1 SER A O   1
ATOM   5 C CB  A SER A 1 ? 11.900 11.300 15.400 0.60 24.00 ? 1 SER A CB  1
ATOM   6 O OG  A SER A 1 ? 11.500 9.950  15.300 0.60 25.00 ? 1 SER A OG  1
ATOM   7 C CB  B SER A 1 ? 11.800 11.200 15.500 0.40 24.00 ? 1 SER A CB  1
ATOM   8 O OG  B SER A 1 ? 13.200 11.100 15.600 0.40 25.00 ? 1 SER A OG  1
HETATM 9 O O   . HOH B . ? 15.000 15.000 15.000 1.00 30.00 ? 101 HOH A O 1
HETATM 10 ZN ZN . ZN C . ? 20.000 20.000 20.000 1.00 30.00 2 102 ZN A ZN 1
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1 0.2500 0.2600 0.2700 0.0100 -0.0200 0.0300
loop_
_atom_type.symbol
C
N
O
ZN
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.details
metalc1 metalc A 1 OG C ZN
;a text field
over two lines
;
_entity.id 1
_entity.pdbx_description \"the protein's 'name'\"
";

/// A source holding owned tables, as another parser might.
struct Tables(HashMap<String, Table>, Vec<String>);

#[derive(Default)]
struct Table {
    category: String,
    tags: Vec<String>,
    rows: Vec<Vec<String>>,
    is_loop: bool,
}

impl CifTable for &Table {
    fn row_count(&self) -> usize {
        self.rows.len()
    }
    fn column(&self, tag: &str) -> Option<usize> {
        self.tags.iter().position(|t| t.eq_ignore_ascii_case(tag))
    }
    fn cell(&self, row: usize, column: usize) -> Cow<'_, str> {
        Cow::Owned(self.rows[row][column].clone())
    }
    fn tags(&self) -> Vec<Cow<'_, str>> {
        self.tags.iter().map(|t| Cow::Owned(format!("_{}.{}", self.category, t))).collect()
    }
    fn is_loop(&self) -> bool {
        self.is_loop
    }
}

impl CifSource for Tables {
    type Table<'a> = &'a Table;
    fn table(&self, category: &str) -> Option<&Table> {
        self.0.get(&category.to_ascii_lowercase())
    }
    fn categories(&self) -> Vec<String> {
        self.1.clone()
    }
}

/// Copy Reduce3's parsed block into owned tables.
fn tables_from_text(text: &str) -> Tables {
    let doc = cif::parse(text);
    let mut out = HashMap::new();
    let mut order = Vec::new();
    for c in &doc.blocks[0].categories {
        let name = c.name.trim_start_matches('_').to_string();
        let t = Table {
            category: name.clone(),
            tags: c.tags.clone(),
            rows: (0..c.nrows()).map(|r| (0..c.ncols()).map(|k| c.get(r, k).to_string()).collect()).collect(),
            is_loop: c.is_loop,
        };
        order.push(name.clone());
        out.insert(name, t);
    }
    Tables(out, order)
}

/// A sink that collects the block into owned tables.
#[derive(Default)]
struct Collect {
    code: String,
    tables: HashMap<String, Table>,
    open: Option<String>,
}

fn split(tag: &str) -> (String, String) {
    let (c, t) = tag.trim_start_matches('_').split_once('.').unwrap();
    (c.to_ascii_lowercase(), t.to_ascii_lowercase())
}

fn text(v: CifCell<'_>) -> String {
    match v {
        CifCell::Text(s) => s.to_string(),
        CifCell::Unknown => "?".into(),
        CifCell::NotApplicable => ".".into(),
    }
}

impl CifSink for Collect {
    fn begin_block(&mut self, code: &str) {
        self.code = code.to_string();
    }
    fn item(&mut self, tag: &str, value: CifCell<'_>) {
        let (c, t) = split(tag);
        let table = self.tables.entry(c).or_default();
        table.tags.push(t);
        if table.rows.is_empty() {
            table.rows.push(Vec::new());
        }
        table.rows[0].push(text(value));
    }
    fn begin_loop(&mut self, tags: &[&str], rows: usize) {
        let (c, _) = split(tags[0]);
        let table =
            Table { category: c.clone(), tags: tags.iter().map(|t| split(t).1).collect(), rows: Vec::with_capacity(rows), is_loop: true };
        self.tables.insert(c.clone(), table);
        self.open = Some(c);
    }
    fn row(&mut self, values: &[CifCell<'_>]) {
        let c = self.open.as_ref().unwrap();
        self.tables.get_mut(c).unwrap().rows.push(values.iter().map(|&v| text(v)).collect());
    }
    fn end_loop(&mut self) {
        self.open = None;
    }
}

#[test]
fn foreign_source_builds_the_same_model() {
    let from_text = mmcif::read_mmcif(MODEL).unwrap();
    let from_tables = mmcif::structure_from_cif(&tables_from_text(MODEL)).unwrap();
    assert_eq!(from_text.atoms_size(), 10);
    assert_eq!(mmcif::write_mmcif(&from_text), mmcif::write_mmcif(&from_tables));
    let first = &from_tables.models[0].chains[0].residue_groups[0];
    assert_eq!(first.atom_groups.len(), 3, "blank, A and B conformers");
    assert!(first.atom_groups[0].atoms[0].uij.is_some());
}

#[test]
fn sink_receives_what_the_text_writer_writes() {
    let st = mmcif::read_mmcif(MODEL).unwrap();
    let text = mmcif::write_mmcif(&st);
    let mut sink = Collect::default();
    mmcif::write_cif(&st, &mut sink);
    assert_eq!(sink.code, "default");
    let atoms = &sink.tables["atom_site"];
    assert_eq!(atoms.rows.len(), 10);
    assert_eq!(sink.tables["atom_site_anisotrop"].rows.len(), 1);
    // reading the collected tables back gives the same model and text
    let again = mmcif::structure_from_cif(&Tables(sink.tables, Vec::new())).unwrap();
    assert_eq!(mmcif::write_mmcif(&again), text);
}

#[test]
fn null_values_and_quoting() {
    let st = mmcif::read_mmcif(MODEL).unwrap();
    let mut sink = Collect::default();
    mmcif::write_cif(&st, &mut sink);
    let atoms = &sink.tables["atom_site"];
    let col = |t: &str| atoms.tags.iter().position(|x| x == t).unwrap();
    assert_eq!(atoms.rows[0][col("label_alt_id")], ".");
    assert_eq!(atoms.rows[0][col("pdbx_pdb_ins_code")], "?");
    assert_eq!(atoms.rows[9][col("pdbx_formal_charge")], "2");
    assert_eq!(sink.tables["symmetry"].rows[0][0], "P 21 21 21");
    assert!(mmcif::write_mmcif(&st).contains("_symmetry.space_group_name_H-M    'P 21 21 21'\n"));
    assert!(reduce3::cifsource::needs_quotes("data_x") && !reduce3::cifsource::needs_quotes("x,y,z"));
}

/// Fixed mode writes the model back into its source block.
#[test]
fn preserving_writer_keeps_the_source_block() {
    use reduce3::cifsource::CifText;
    use reduce3::geom::v3;
    use reduce3::model::Atom;
    let mut st = mmcif::read_mmcif(MODEL).unwrap();
    // a new hydrogen on SER 1 and a moved OG of conformer A
    let rg = &mut st.models[0].chains[0].residue_groups[0];
    let mut h = Atom::new(" HA ", "H", v3(11.9, 11.6, 13.3));
    h.occ = 1.0;
    h.b = 21.5;
    rg.atom_groups[0].atoms.push(h);
    rg.atom_groups[1].atoms[1].xyz = v3(11.25, 9.9, 15.3);
    let doc = cif::parse(MODEL);
    let mut w = CifText::with_capacity(4096);
    mmcif::write_cif_preserving(&st, &doc.blocks[0], "test", &mut w).unwrap();
    let text = w.out;
    let out = cif::parse(&text);
    let b = &out.blocks[0];
    assert_eq!(b.name, "test");
    // every source category, in order, the untouched ones unchanged
    let names: Vec<&str> = b.categories.iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["_cell", "_symmetry", "_atom_site", "_atom_site_anisotrop", "_atom_type", "_struct_conn", "_entity"]);
    let src = &doc.blocks[0];
    for cat in ["_cell", "_symmetry", "_struct_conn", "_entity"] {
        assert_eq!(b.category(cat).unwrap().values, src.category(cat).unwrap().values, "{}", cat);
    }
    assert_eq!(b.item("_struct_conn.details"), Some("a text field\nover two lines"));
    assert_eq!(b.item("_entity.pdbx_description"), Some("the protein's 'name'"));
    let atoms = b.category("_atom_site").unwrap();
    assert_eq!(atoms.nrows(), 11);
    let ids: Vec<&str> = (0..11).map(|r| atoms.get_tag(r, "id").unwrap()).collect();
    assert_eq!(ids, ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]);
    // the new H after the blank-altloc atoms, with its residue's labels
    let row = (0..11).find(|&r| atoms.get_tag(r, "label_atom_id") == Some("HA")).unwrap();
    assert_eq!(row, 4);
    for (tag, value) in [
        ("type_symbol", "H"),
        ("label_alt_id", "."),
        ("label_comp_id", "SER"),
        ("label_asym_id", "A"),
        ("label_seq_id", "1"),
        ("auth_seq_id", "1"),
        ("pdbx_formal_charge", "?"),
        ("Cartn_x", "11.900"),
        ("B_iso_or_equiv", "21.50"),
    ] {
        assert_eq!(atoms.get_tag(row, tag), Some(value), "{}", tag);
    }
    // unchanged input values keep their spelling; a changed coordinate is rewritten
    assert_eq!(atoms.get_tag(0, "Cartn_x"), Some("10.000"));
    let og_a = (0..11).find(|&r| atoms.get_tag(r, "label_atom_id") == Some("OG") && atoms.get_tag(r, "label_alt_id") == Some("A")).unwrap();
    assert_eq!(atoms.get_tag(og_a, "Cartn_x"), Some("11.250"));
    assert_eq!(atoms.get_tag(og_a, "Cartn_z"), Some("15.300"));
    // anisotropic row follows its atom's new id; H joins the atom types
    assert_eq!(b.item("_atom_site_anisotrop.id"), Some("1"));
    let types = b.category("_atom_type").unwrap();
    assert_eq!((0..types.nrows()).map(|r| types.get(r, 0)).collect::<Vec<_>>(), ["C", "N", "O", "ZN", "H"]);
}
