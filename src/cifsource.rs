//! Working on CIF data another parser already holds in memory.
//!
//! [`CifSource`] is the read side: Reduce3 builds its structure from any parsed
//! document that can hand out categories as tables, without re-reading text.
//! [`CifSink`] is the write side: Reduce3 streams its output model as items and
//! loops, which a text writer formats and an in-memory document builder can
//! store directly. Reduce3's own mmCIF reader and writer use the same two traits,
//! so every path builds and emits exactly the same data.

use crate::cif;
use std::borrow::Cow;

/// One category of a parsed CIF block: a loop, or the category's scalar items
/// as a single row.
pub trait CifTable {
    /// Number of rows.
    fn row_count(&self) -> usize;
    /// Column of the item `tag`, given without the category (for example
    /// `Cartn_x`). Names compare case-insensitively.
    fn column(&self, tag: &str) -> Option<usize>;
    /// A value as logical CIF text without delimiters, with `?` for unknown and
    /// `.` for not applicable.
    fn cell(&self, row: usize, column: usize) -> Cow<'_, str>;
    /// A numeric value, ignoring a standard uncertainty such as `1.23(4)`;
    /// `None` for nulls and text that is not a number. Sources that already
    /// hold parsed numbers can return them directly; they must equal what
    /// Rust's `str::parse::<f64>` gives for the same text.
    fn number(&self, row: usize, column: usize) -> Option<f64> {
        cif::parse_f64(&self.cell(row, column))
    }
}

/// One parsed CIF data block.
pub trait CifSource {
    /// The table type the source hands out.
    type Table<'a>: CifTable
    where
        Self: 'a;
    /// The category `category`, given without the leading underscore (for
    /// example `atom_site`), or `None` when the block has no such category.
    fn table(&self, category: &str) -> Option<Self::Table<'_>>;
}

/// A value emitted by Reduce3.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CifCell<'a> {
    /// A present value, without CIF quoting. The text writer quotes it only
    /// when [`needs_quotes`] says so; a bare value keeps its usual CIF typing
    /// (numbers stay numbers), and a sink that types values should do the same.
    Text(&'a str),
    /// CIF unknown (`?`).
    Unknown,
    /// CIF not applicable (`.`).
    NotApplicable,
}

/// Receives an output data block in order.
pub trait CifSink {
    /// The data block code (`data_<code>`); called once, first.
    fn begin_block(&mut self, code: &str);
    /// A scalar item; `tag` is the full name with its leading underscore.
    fn item(&mut self, tag: &str, value: CifCell<'_>);
    /// The start of a loop with these full tags and the given number of rows.
    fn begin_loop(&mut self, tags: &[&str], rows: usize);
    /// One loop row, with one value per tag.
    fn row(&mut self, values: &[CifCell<'_>]);
    /// The end of the current loop.
    fn end_loop(&mut self);
}

/// True when the text writer puts `value` in quotes: it is empty, starts with
/// a character that would make it something else, contains blanks, or is a
/// reserved word.
pub fn needs_quotes(value: &str) -> bool {
    value.is_empty()
        || value.starts_with(['_', '#', '$', '\'', '"', '[', ']', ';'])
        || value.contains([' ', '\t'])
        || matches!(value.to_ascii_lowercase().as_str(), "loop_" | "stop_" | "global_")
        || value.get(..5).is_some_and(|p| p.eq_ignore_ascii_case("data_") || p.eq_ignore_ascii_case("save_"))
}

/// The `CifSink` behind `write_mmcif`: CIF text in the layout iotbx uses.
pub struct CifText {
    /// The text written so far.
    pub out: String,
}

impl CifText {
    /// An empty writer with room for `capacity` bytes.
    pub fn with_capacity(capacity: usize) -> CifText {
        CifText { out: String::with_capacity(capacity) }
    }

    fn value(&mut self, v: CifCell<'_>) {
        match v {
            CifCell::Unknown => self.out.push('?'),
            CifCell::NotApplicable => self.out.push('.'),
            CifCell::Text(s) if !needs_quotes(s) => self.out.push_str(s),
            CifCell::Text(s) => {
                let q = if s.contains('\'') { '"' } else { '\'' };
                self.out.push(q);
                self.out.push_str(s);
                self.out.push(q);
            }
        }
    }
}

impl CifSink for CifText {
    fn begin_block(&mut self, code: &str) {
        self.out.push_str("data_");
        self.out.push_str(code);
        self.out.push('\n');
    }
    fn item(&mut self, tag: &str, value: CifCell<'_>) {
        self.out.push_str(tag);
        for _ in tag.len()..34 {
            self.out.push(' ');
        }
        self.value(value);
        self.out.push('\n');
    }
    fn begin_loop(&mut self, tags: &[&str], _rows: usize) {
        self.out.push_str("loop_\n");
        for t in tags {
            self.out.push_str("  ");
            self.out.push_str(t);
            self.out.push('\n');
        }
    }
    fn row(&mut self, values: &[CifCell<'_>]) {
        self.out.push_str("  ");
        for &v in values {
            self.out.push(' ');
            self.value(v);
        }
        self.out.push('\n');
    }
    fn end_loop(&mut self) {
        self.out.push('\n');
    }
}

impl<'b> CifTable for &cif::Category<'b> {
    fn row_count(&self) -> usize {
        self.nrows()
    }
    fn column(&self, tag: &str) -> Option<usize> {
        self.col(tag)
    }
    fn cell(&self, row: usize, column: usize) -> Cow<'_, str> {
        Cow::Borrowed(self.get(row, column))
    }
}

impl<'b> CifSource for cif::Block<'b> {
    type Table<'a>
        = &'a cif::Category<'b>
    where
        Self: 'a;
    fn table(&self, category: &str) -> Option<&cif::Category<'b>> {
        self.category(&format!("_{}", category))
    }
}
