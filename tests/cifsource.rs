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
";

/// A source holding owned tables, as another parser might.
struct Tables(HashMap<String, Table>);

#[derive(Default)]
struct Table {
    tags: Vec<String>,
    rows: Vec<Vec<String>>,
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
}

impl CifSource for Tables {
    type Table<'a> = &'a Table;
    fn table(&self, category: &str) -> Option<&Table> {
        self.0.get(&category.to_ascii_lowercase())
    }
}

/// Copy Reduce3's parsed block into owned tables.
fn tables_from_text(text: &str) -> Tables {
    let doc = cif::parse(text);
    let mut out = HashMap::new();
    for c in &doc.blocks[0].categories {
        let t = Table {
            tags: c.tags.clone(),
            rows: (0..c.nrows()).map(|r| (0..c.ncols()).map(|k| c.get(r, k).to_string()).collect()).collect(),
        };
        out.insert(c.name.trim_start_matches('_').to_string(), t);
    }
    Tables(out)
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
        let table = Table { tags: tags.iter().map(|t| split(t).1).collect(), rows: Vec::with_capacity(rows) };
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
    let again = mmcif::structure_from_cif(&Tables(sink.tables)).unwrap();
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
