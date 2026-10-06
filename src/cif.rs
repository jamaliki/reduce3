//! Minimal, fast, zero-copy CIF (STAR) reader.
//!
//! Handles data blocks, single items, `loop_` tables, quoted strings and
//! semicolon text fields, which covers mmCIF model files, the CCP4/geostd
//! monomer library and the wwPDB chemical component dictionary. Category and
//! item names are matched case-insensitively.

use rustc_hash::FxHashMap;

#[derive(Debug)]
pub struct Category<'a> {
    /// Lower-case category name including the leading underscore, e.g. `_atom_site`.
    pub name: String,
    /// Lower-case item names (without the category prefix).
    pub tags: Vec<String>,
    /// Row-major values; `values.len() == rows * tags.len()`.
    pub values: Vec<&'a str>,
    /// The full tags as the source spells them (`_atom_site.Cartn_x`).
    pub spelling: Vec<&'a str>,
    /// Written as a `loop_` (otherwise as scalar items).
    pub is_loop: bool,
}

impl<'a> Category<'a> {
    #[inline]
    pub fn ncols(&self) -> usize {
        self.tags.len()
    }
    #[inline]
    pub fn nrows(&self) -> usize {
        if self.tags.is_empty() { 0 } else { self.values.len() / self.tags.len() }
    }
    pub fn col(&self, tag: &str) -> Option<usize> {
        let t = tag.to_ascii_lowercase();
        self.tags.iter().position(|x| *x == t)
    }
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> &'a str {
        self.values[row * self.tags.len() + col]
    }
    /// Value at (row, tag) or None when the column is absent.
    pub fn get_tag(&self, row: usize, tag: &str) -> Option<&'a str> {
        self.col(tag).map(|c| self.get(row, c))
    }
}

#[derive(Debug)]
pub struct Block<'a> {
    pub name: &'a str,
    pub categories: Vec<Category<'a>>,
    index: FxHashMap<String, usize>,
}

impl<'a> Block<'a> {
    pub fn category(&self, name: &str) -> Option<&Category<'a>> {
        let n = name.to_ascii_lowercase();
        self.index.get(&n).map(|&i| &self.categories[i])
    }
    /// A single item value such as `_cell.length_a`.
    pub fn item(&self, full: &str) -> Option<&'a str> {
        let (cat, tag) = split_tag(full);
        let c = self.category(&cat)?;
        let col = c.col(&tag)?;
        if c.nrows() == 0 { None } else { Some(c.get(0, col)) }
    }
}

#[derive(Debug, Default)]
pub struct Document<'a> {
    pub blocks: Vec<Block<'a>>,
}

impl<'a> Document<'a> {
    pub fn block(&self, name: &str) -> Option<&Block<'a>> {
        self.blocks.iter().find(|b| b.name.eq_ignore_ascii_case(name))
    }
}

/// True for CIF null markers.
#[inline]
pub fn is_null(v: &str) -> bool {
    v == "?" || v == "."
}

/// Split `_cat.item` into (`_cat`, `item`), lower-cased. Tags without a dot (old
/// style) are treated as a category with a single unnamed item.
fn split_tag(full: &str) -> (String, String) {
    let lower = full.to_ascii_lowercase();
    match lower.find('.') {
        Some(p) => (lower[..p].to_string(), lower[p + 1..].to_string()),
        None => (lower, String::new()),
    }
}

#[derive(Debug)]
enum Tok<'a> {
    Value(&'a str),
    /// A value written in quotes or as a text field; never a keyword.
    Quoted(&'a str),
    Tag(&'a str),
    Loop,
    Data(&'a str),
    Save,
}

struct Lexer<'a> {
    s: &'a str,
    b: &'a [u8],
    pos: usize,
}

impl<'a> Lexer<'a> {
    fn new(s: &'a str) -> Self {
        Lexer { s, b: s.as_bytes(), pos: 0 }
    }

    #[inline]
    fn at_line_start(&self, p: usize) -> bool {
        p == 0 || self.b[p - 1] == b'\n' || self.b[p - 1] == b'\r'
    }

    fn next(&mut self) -> Option<Tok<'a>> {
        let b = self.b;
        let n = b.len();
        loop {
            // skip whitespace
            while self.pos < n && (b[self.pos] as char).is_ascii_whitespace() {
                self.pos += 1;
            }
            if self.pos >= n {
                return None;
            }
            let c = b[self.pos];
            if c == b'#' {
                while self.pos < n && b[self.pos] != b'\n' {
                    self.pos += 1;
                }
                continue;
            }
            if c == b';' && self.at_line_start(self.pos) {
                // text field: runs until a line starting with ';'
                let start = self.pos + 1;
                let mut p = start;
                loop {
                    match memchr_nl(&b[p..]) {
                        Some(off) => {
                            p += off + 1;
                            if p < n && b[p] == b';' {
                                let mut end = p - 1; // the newline before ';'
                                if end > start && b[end - 1] == b'\r' {
                                    end -= 1;
                                }
                                self.pos = p + 1;
                                return Some(Tok::Quoted(&self.s[start..end.max(start)]));
                            }
                        }
                        None => {
                            self.pos = n;
                            return Some(Tok::Quoted(&self.s[start..n]));
                        }
                    }
                }
            }
            if c == b'\'' || c == b'"' {
                let start = self.pos + 1;
                let mut p = start;
                while p < n {
                    if b[p] == c && (p + 1 >= n || (b[p + 1] as char).is_ascii_whitespace()) {
                        self.pos = p + 1;
                        return Some(Tok::Quoted(&self.s[start..p]));
                    }
                    if b[p] == b'\n' {
                        // unterminated quote: treat the rest of the line as the value
                        self.pos = p;
                        return Some(Tok::Quoted(&self.s[start..p]));
                    }
                    p += 1;
                }
                self.pos = n;
                return Some(Tok::Quoted(&self.s[start..n]));
            }
            let start = self.pos;
            while self.pos < n && !(b[self.pos] as char).is_ascii_whitespace() {
                self.pos += 1;
            }
            let t = &self.s[start..self.pos];
            if c == b'_' {
                return Some(Tok::Tag(t));
            }
            if t.len() >= 5 {
                let lower5 = t[..5].to_ascii_lowercase();
                if lower5 == "data_" {
                    return Some(Tok::Data(&t[5..]));
                }
                if lower5 == "loop_" && t.len() == 5 {
                    return Some(Tok::Loop);
                }
                if lower5 == "save_" {
                    return Some(Tok::Save);
                }
            }
            return Some(Tok::Value(t));
        }
    }
}

#[inline]
fn memchr_nl(b: &[u8]) -> Option<usize> {
    b.iter().position(|&c| c == b'\n')
}

struct BlockBuilder<'a> {
    name: &'a str,
    categories: Vec<Category<'a>>,
    index: FxHashMap<String, usize>,
}

impl<'a> BlockBuilder<'a> {
    fn new(name: &'a str) -> Self {
        BlockBuilder { name, categories: Vec::new(), index: FxHashMap::default() }
    }
    fn add_item(&mut self, tag: &'a str, value: &'a str) {
        let (cat, item) = split_tag(tag);
        match self.index.get(&cat) {
            Some(&i) => {
                let c = &mut self.categories[i];
                if c.nrows() <= 1 {
                    if let Some(col) = c.tags.iter().position(|x| *x == item) {
                        c.values[col] = value;
                    } else {
                        c.tags.push(item);
                        c.values.push(value);
                        c.spelling.push(tag);
                    }
                }
            }
            None => {
                self.index.insert(cat.clone(), self.categories.len());
                self.categories.push(Category { name: cat, tags: vec![item], values: vec![value], spelling: vec![tag], is_loop: false });
            }
        }
    }
    fn add_loop(&mut self, tags: Vec<&'a str>, values: Vec<&'a str>) {
        if tags.is_empty() {
            return;
        }
        let (cat, _) = split_tag(tags[0]);
        let items: Vec<String> = tags.iter().map(|t| split_tag(t).1).collect();
        let ncol = items.len();
        let mut values = values;
        let rem = values.len() % ncol;
        if rem != 0 {
            values.truncate(values.len() - rem);
        }
        match self.index.get(&cat) {
            Some(&i) if self.categories[i].tags == items => {
                self.categories[i].values.extend(values);
            }
            _ => {
                self.index.insert(cat.clone(), self.categories.len());
                self.categories.push(Category { name: cat, tags: items, values, spelling: tags, is_loop: true });
            }
        }
    }
    fn finish(self) -> Block<'a> {
        Block { name: self.name, categories: self.categories, index: self.index }
    }
}

/// Parse a CIF document. Malformed constructs are skipped rather than fatal.
pub fn parse(text: &str) -> Document<'_> {
    let mut lx = Lexer::new(text);
    let mut doc = Document::default();
    let mut cur: Option<BlockBuilder> = None;
    let mut pending: Option<Tok> = None;
    loop {
        let tok = match pending.take() {
            Some(t) => t,
            None => match lx.next() {
                Some(t) => t,
                None => break,
            },
        };
        match tok {
            Tok::Data(name) => {
                if let Some(b) = cur.take() {
                    doc.blocks.push(b.finish());
                }
                cur = Some(BlockBuilder::new(name));
            }
            Tok::Tag(tag) => {
                let v = match lx.next() {
                    Some(Tok::Value(v)) | Some(Tok::Quoted(v)) => v,
                    other => {
                        pending = other;
                        continue;
                    }
                };
                cur.get_or_insert_with(|| BlockBuilder::new("")).add_item(tag, v);
            }
            Tok::Loop => {
                let mut tags = Vec::new();
                let mut values = Vec::new();
                loop {
                    match lx.next() {
                        Some(Tok::Tag(t)) if values.is_empty() => tags.push(t),
                        Some(Tok::Value(v)) | Some(Tok::Quoted(v)) => values.push(v),
                        other => {
                            pending = other;
                            break;
                        }
                    }
                }
                cur.get_or_insert_with(|| BlockBuilder::new("")).add_loop(tags, values);
            }
            Tok::Save | Tok::Value(_) | Tok::Quoted(_) => {}
        }
    }
    if let Some(b) = cur.take() {
        doc.blocks.push(b.finish());
    }
    doc
}

/// Parse a CIF number, ignoring a trailing standard-uncertainty suffix like `1.23(4)`.
pub fn parse_f64(v: &str) -> Option<f64> {
    if is_null(v) {
        return None;
    }
    let v = match v.find('(') {
        Some(p) => &v[..p],
        None => v,
    };
    v.trim().parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_loops_and_items() {
        let txt = "data_x\n_cell.length_a 10.5\n_cell.length_b 'a b'\nloop_\n_atom.id\n_atom.name\n1 'C A'\n2 \"O\"\n;text\nfield\n;\nZ\n";
        let d = parse(txt);
        let b = &d.blocks[0];
        assert_eq!(b.name, "x");
        assert_eq!(b.item("_cell.length_a"), Some("10.5"));
        assert_eq!(b.item("_cell.length_b"), Some("a b"));
        let c = b.category("_atom").unwrap();
        assert_eq!(c.nrows(), 3);
        assert_eq!(c.get(0, 1), "C A");
        assert_eq!(c.get(2, 0), "text\nfield");
        assert_eq!(c.get(2, 1), "Z");
    }

    #[test]
    fn quote_inside_token() {
        let d = parse("data_y\nloop_\n_a.n\n_a.m\nO5' x\n\"C1'\" y\n");
        let c = d.blocks[0].category("_a").unwrap();
        assert_eq!(c.get(0, 0), "O5'");
        assert_eq!(c.get(1, 0), "C1'");
    }
}
