//! Unit cells exactly as cctbx computes them (`uctbx::unit_cell`,
//! `sgtbx::space_group::average_unit_cell`, `iotbx.pdb.cryst1_interpretation`).
//!
//! Reduce2's model processing moves every site through fractional space and
//! back (`site_symmetry_table.apply_symmetry_sites`, twice), which perturbs
//! coordinates in the last bits. Compat mode repeats that so that exact score
//! ties break the same way. On arm64 the cctbx build fuses the multiply-adds
//! of `fractionalize`/`orthogonalize`; the same fused forms are used here.

use crate::geom::{v3, Vec3};
use crate::model::Structure;
use crate::spacegroup_data::{ROT_SETS, SG_INFO, SYMBOLS};

const PI_180: f64 = std::f64::consts::PI / 180.0;

#[derive(Clone, Debug)]
pub struct UnitCell {
    pub params: [f64; 6],
    pub orth: [f64; 9],
    pub frac: [f64; 9],
    pub metr: [f64; 6],
    pub volume: f64,
}

#[inline]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    if cfg!(target_arch = "aarch64") { a.mul_add(b, c) } else { a * b + c }
}

impl UnitCell {
    /// `unit_cell(parameters)`.
    pub fn new(p: [f64; 6]) -> Option<UnitCell> {
        if p.iter().any(|&x| !(x > 0.0)) {
            return None;
        }
        let mut cos_ang = [0.0; 3];
        let mut sin_ang = [0.0; 3];
        for i in 3..6 {
            if p[i] >= 180.0 {
                return None;
            }
            let a = p[i] * PI_180;
            cos_ang[i - 3] = a.cos();
            sin_ang[i - 3] = a.sin();
            if sin_ang[i - 3] == 0.0 {
                return None;
            }
        }
        let mut d = 1.0;
        for i in 0..3 {
            d -= cos_ang[i] * cos_ang[i];
        }
        d += 2.0 * cos_ang[0] * cos_ang[1] * cos_ang[2];
        if d < 0.0 {
            return None;
        }
        let mut r_cos = [0.0; 3];
        for i in 0..3 {
            let denom = sin_ang[(i + 1) % 3] * sin_ang[(i + 2) % 3];
            if denom == 0.0 {
                return None;
            }
            r_cos[i] = (cos_ang[(i + 1) % 3] * cos_ang[(i + 2) % 3] - cos_ang[i]) / denom;
        }
        for c in r_cos.iter_mut() {
            if *c < -1.0 || *c > 1.0 {
                return None;
            }
            *c = c.acos().cos();
        }
        let s1rca2 = (1.0 - r_cos[0] * r_cos[0]).sqrt();
        if s1rca2 == 0.0 {
            return None;
        }
        let orth = [
            p[0],
            cos_ang[2] * p[1],
            cos_ang[1] * p[2],
            0.0,
            sin_ang[2] * p[1],
            -sin_ang[1] * r_cos[0] * p[2],
            0.0,
            0.0,
            sin_ang[1] * p[2] * s1rca2,
        ];
        let frac = [
            1.0 / p[0],
            -cos_ang[2] / (sin_ang[2] * p[0]),
            -(cos_ang[2] * sin_ang[1] * r_cos[0] + cos_ang[1] * sin_ang[2]) / (sin_ang[1] * s1rca2 * sin_ang[2] * p[0]),
            0.0,
            1.0 / (sin_ang[2] * p[1]),
            r_cos[0] / (s1rca2 * sin_ang[2] * p[1]),
            0.0,
            0.0,
            1.0 / (sin_ang[1] * s1rca2 * p[2]),
        ];
        let metr = [
            p[0] * p[0],
            p[1] * p[1],
            p[2] * p[2],
            p[0] * p[1] * cos_ang[2],
            p[0] * p[2] * cos_ang[1],
            p[1] * p[2] * cos_ang[0],
        ];
        let volume = p[0] * p[1] * p[2] * d.sqrt();
        Some(UnitCell { params: p, orth, frac, metr, volume })
    }

    /// `unit_cell(metrical_matrix)`.
    pub fn from_metrical(g: [f64; 6]) -> Option<UnitCell> {
        if g[0] <= 0.0 || g[1] <= 0.0 || g[2] <= 0.0 {
            return None;
        }
        let a = g[0].sqrt();
        let b = g[1].sqrt();
        let c = g[2].sqrt();
        let acos_deg = |x: f64| x.acos() / PI_180;
        UnitCell::new([a, b, c, acos_deg(g[5] / b / c), acos_deg(g[4] / c / a), acos_deg(g[3] / a / b)])
    }

    /// `space_group.average_unit_cell(self)` for the given rotation parts.
    pub fn averaged(&self, rots: &[i8]) -> Option<UnitCell> {
        let t = self.metr;
        let n = rots.len() / 9;
        let mut r = [0.0f64; 6];
        for k in 0..n {
            let c: Vec<f64> = rots[9 * k..9 * k + 9].iter().map(|&x| x as f64).collect();
            // sym_mat3::tensor_transpose_transform(c): c^T * t * c
            let ctt = [
                c[0] * t[0] + c[3] * t[3] + c[6] * t[4],
                c[3] * t[1] + c[0] * t[3] + c[6] * t[5],
                c[6] * t[2] + c[0] * t[4] + c[3] * t[5],
                c[1] * t[0] + c[4] * t[3] + c[7] * t[4],
                c[4] * t[1] + c[1] * t[3] + c[7] * t[5],
                c[7] * t[2] + c[1] * t[4] + c[4] * t[5],
                c[2] * t[0] + c[5] * t[3] + c[8] * t[4],
                c[5] * t[1] + c[2] * t[3] + c[8] * t[5],
                c[8] * t[2] + c[2] * t[4] + c[5] * t[5],
            ];
            let s = [
                ctt[0] * c[0] + ctt[1] * c[3] + ctt[2] * c[6],
                ctt[3] * c[1] + ctt[4] * c[4] + ctt[5] * c[7],
                ctt[6] * c[2] + ctt[7] * c[5] + ctt[8] * c[8],
                ctt[0] * c[1] + ctt[1] * c[4] + ctt[2] * c[7],
                ctt[0] * c[2] + ctt[1] * c[5] + ctt[2] * c[8],
                ctt[3] * c[2] + ctt[4] * c[5] + ctt[5] * c[8],
            ];
            for i in 0..6 {
                r[i] += s[i];
            }
        }
        for x in r.iter_mut() {
            *x /= n as f64;
        }
        UnitCell::from_metrical(r)
    }

    /// `unit_cell.fractionalize(sites_cart)` (upper-triangular form).
    #[inline]
    pub fn fractionalize(&self, x: Vec3) -> Vec3 {
        let f = &self.frac;
        v3(fma(f[2], x.z, fma(f[0], x.x, f[1] * x.y)), fma(f[4], x.y, f[5] * x.z), f[8] * x.z)
    }

    /// `unit_cell.orthogonalize(sites_frac)` (upper-triangular form).
    #[inline]
    pub fn orthogonalize(&self, x: Vec3) -> Vec3 {
        let o = &self.orth;
        v3(fma(o[2], x.z, fma(o[0], x.x, o[1] * x.y)), fma(o[4], x.y, o[5] * x.z), o[8] * x.z)
    }
}

impl UnitCell {
    /// `sym_mat3::tensor_transform(c)`: c * t * c^T (via `mat3 * sym_mat3`).
    fn tensor_transform(t: &[f64; 6], c: &[f64; 9]) -> [f64; 6] {
        let f3 = |a: f64, b: f64, c_: f64, d: f64, e: f64, f: f64| fma(e, f, fma(a, b, c_ * d));
        let (l, r) = (c, t);
        let ct = [
            f3(l[0], r[0], l[1], r[3], l[2], r[4]),
            f3(l[1], r[1], l[0], r[3], l[2], r[5]),
            f3(l[2], r[2], l[0], r[4], l[1], r[5]),
            f3(l[3], r[0], l[4], r[3], l[5], r[4]),
            f3(l[4], r[1], l[3], r[3], l[5], r[5]),
            f3(l[5], r[2], l[3], r[4], l[4], r[5]),
            f3(l[6], r[0], l[7], r[3], l[8], r[4]),
            f3(l[7], r[1], l[6], r[3], l[8], r[5]),
            f3(l[8], r[2], l[6], r[4], l[7], r[5]),
        ];
        [
            f3(ct[0], c[0], ct[1], c[1], ct[2], c[2]),
            f3(ct[3], c[3], ct[4], c[4], ct[5], c[5]),
            f3(ct[6], c[6], ct[7], c[7], ct[8], c[8]),
            f3(ct[0], c[3], ct[1], c[4], ct[2], c[5]),
            f3(ct[0], c[6], ct[1], c[7], ct[2], c[8]),
            f3(ct[3], c[6], ct[4], c[7], ct[5], c[8]),
        ]
    }

    /// U_cart -> U* -> U_cart (`adptbx.u_cart_as_u_star`, `u_star_as_u_cart`),
    /// which the X-ray structure round trip applies to anisotropic atoms.
    pub fn u_cart_roundtrip(&self, u: [f64; 6]) -> [f64; 6] {
        let us = Self::tensor_transform(&u, &self.frac);
        Self::tensor_transform(&us, &self.orth)
    }
}

/// What `crystal.symmetry.as_cif_block` writes for a space group.
pub struct SpaceGroupInfo {
    pub hall: &'static str,
    pub number: u16,
    pub crystal_system: &'static str,
    pub hm: &'static str,
    pub ops: &'static [&'static str],
}

/// Rotation parts of the space group named by a CRYST1 symbol, following
/// `iotbx.pdb.cryst1_interpretation` (rhombohedral and short monoclinic
/// symbols resolved with the cell). None when cctbx would not know it.
pub fn cryst1_rotations(symbol: &str, p: &[f64; 6]) -> Option<&'static [i8]> {
    lookup(symbol, p).map(|(r, _)| ROT_SETS[r])
}

/// The space group named by a CRYST1 symbol (see `cryst1_rotations`).
pub fn space_group_info(symbol: &str, p: &[f64; 6]) -> Option<SpaceGroupInfo> {
    lookup(symbol, p).map(|(_, i)| {
        let (hall, number, crystal_system, hm, ops) = SG_INFO[i];
        SpaceGroupInfo { hall, number, crystal_system, hm, ops }
    })
}

fn lookup(symbol: &str, p: &[f64; 6]) -> Option<(usize, usize)> {
    let s = symbol.trim().replace(' ', "").to_ascii_uppercase();
    let is90 = |a: f64| (a - 90.0).abs() < 0.01;
    let is120 = |a: f64| (a - 120.0).abs() < 0.01;
    let equiv = |r: f64, s: f64, t: f64| {
        let m = (r + s + t) / 3.0;
        (r - m).abs() < 0.01 && (s - m).abs() < 0.01 && (t - m).abs() < 0.01
    };
    let rhombo = [
        ("R3", "R3"), ("H3", "R3"), ("R-3", "R-3"), ("H-3", "R-3"), ("R32", "R32"), ("H32", "R32"),
        ("R3M", "R3M"), ("H3M", "R3M"), ("R3C", "R3C"), ("H3C", "R3C"), ("R-3M", "R-3M"), ("H-3M", "R-3M"),
        ("R-3C", "R-3C"), ("H-3C", "R-3C"),
    ];
    let key = if let Some(&(_, r)) = rhombo.iter().find(|(k, _)| *k == s) {
        let (a, b, c, al, be, ga) = (p[0], p[1], p[2], p[3], p[4], p[5]);
        if (a - b).abs() <= 0.01 && is90(al) && is90(be) && is120(ga) {
            format!("{}:H", r)
        } else if equiv(a, b, c) && equiv(al, be, ga) {
            format!("{}:R", r)
        } else {
            return None;
        }
    } else if ["P2", "P21", "C2", "A2", "B2", "I2"].contains(&s.as_str()) {
        let (z, t) = s.split_at(1);
        if is90(p[3]) && is90(p[5]) {
            if z == "B" {
                return None;
            }
            format!("{}1{}1", z, t)
        } else if is90(p[3]) && is90(p[4]) {
            if z == "C" {
                return None;
            }
            format!("{}11{}", z, t)
        } else {
            return None;
        }
    } else {
        s
    };
    SYMBOLS.binary_search_by(|(k, _, _)| (*k).cmp(key.as_str())).ok().map(|i| (SYMBOLS[i].1 as usize, SYMBOLS[i].2 as usize))
}

/// The cell cctbx processes the model with: the symmetry-averaged CRYST1
/// cell, or, without a usable one, the P1 box Reduce2 puts around the model
/// (`add_crystal_symmetry_if_necessary(box_cushion=5)`).
pub fn processing_cell(st: &Structure) -> Option<UnitCell> {
    if let Some(cs) = &st.crystal {
        let p = cs.cell;
        let dummy_len = p[..3].iter().all(|&v| v == 0.0 || v == 1.0);
        let dummy_ang = p[3..].iter().all(|&v| v == 0.0 || v == 90.0);
        let sg = cs.space_group.replace(' ', "");
        let dummy = dummy_len && dummy_ang && (sg.is_empty() || sg == "P1");
        if !dummy {
            if let Some(rots) = cryst1_rotations(&cs.space_group, &p) {
                if let Some(uc) = UnitCell::new(p).and_then(|u| u.averaged(rots)) {
                    return Some(uc);
                }
            }
        }
    }
    // box around all atoms
    let mut mn = v3(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut mx = v3(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for m in &st.models {
        for c in &m.chains {
            for rg in &c.residue_groups {
                for ag in &rg.atom_groups {
                    for a in &ag.atoms {
                        any = true;
                        mn = v3(mn.x.min(a.xyz.x), mn.y.min(a.xyz.y), mn.z.min(a.xyz.z));
                        mx = v3(mx.x.max(a.xyz.x), mx.y.max(a.xyz.y), mx.z.max(a.xyz.z));
                    }
                }
            }
        }
    }
    if !any {
        return None;
    }
    let p = [(mx.x - mn.x) + 10.0, (mx.y - mn.y) + 10.0, (mx.z - mn.z) + 10.0, 90.0, 90.0, 90.0];
    let p1 = cryst1_rotations("P 1", &p)?;
    UnitCell::new(p).and_then(|u| u.averaged(p1))
}

/// Move every site through fractional coordinates and back `times` times.
pub fn roundtrip_sites(pos: &mut [Vec3], uc: &UnitCell, times: usize) {
    for p in pos.iter_mut() {
        let mut x = *p;
        for _ in 0..times {
            x = uc.orthogonalize(uc.fractionalize(x));
        }
        *p = x;
    }
}

/// A crystallographic symmetry operator in fractional coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymOp {
    pub r: [[f64; 3]; 3],
    pub t: [f64; 3],
}

impl SymOp {
    /// Parse an operator such as `-x+1/2,y,-z` or `x-y,x,z+1/6`.
    pub fn parse(s: &str) -> Option<SymOp> {
        let mut op = SymOp { r: [[0.0; 3]; 3], t: [0.0; 3] };
        let rows: Vec<&str> = s.split(',').collect();
        if rows.len() != 3 {
            return None;
        }
        for (k, row) in rows.iter().enumerate() {
            let b = row.trim().as_bytes();
            let mut i = 0;
            while i < b.len() {
                let sign = match b[i] {
                    b'-' => {
                        i += 1;
                        -1.0
                    }
                    b'+' => {
                        i += 1;
                        1.0
                    }
                    _ => 1.0,
                };
                let start = i;
                while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'/' || b[i] == b'.') {
                    i += 1;
                }
                let number = &row.trim()[start..i];
                if i < b.len() && matches!(b[i], b'x' | b'y' | b'z' | b'X' | b'Y' | b'Z') {
                    let c = if number.is_empty() { 1.0 } else { number.parse::<f64>().ok()? };
                    let axis = (b[i].to_ascii_lowercase() - b'x') as usize;
                    op.r[k][axis] += sign * c;
                    i += 1;
                } else if !number.is_empty() {
                    let v = match number.split_once('/') {
                        Some((p, q)) => p.parse::<f64>().ok()? / q.parse::<f64>().ok()?,
                        None => number.parse::<f64>().ok()?,
                    };
                    op.t[k] += sign * v;
                } else {
                    return None;
                }
            }
        }
        Some(op)
    }

    /// The image of a fractional site.
    pub fn apply(&self, f: Vec3) -> Vec3 {
        let a = [f.x, f.y, f.z];
        let row = |k: usize| self.r[k][0] * a[0] + self.r[k][1] * a[1] + self.r[k][2] * a[2] + self.t[k];
        v3(row(0), row(1), row(2))
    }

    pub fn is_identity(&self) -> bool {
        self.r == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] && self.t == [0.0; 3]
    }
}

/// The model's crystal symmetry: the symmetry-averaged cell and the space
/// group operators. None without a usable CRYST1 (cctbx then boxes the model
/// in P1, which has no contacts between copies).
pub fn crystal_operators(st: &Structure) -> Option<(UnitCell, Vec<SymOp>)> {
    let cs = st.crystal.as_ref()?;
    let p = cs.cell;
    let dummy_len = p[..3].iter().all(|&v| v == 0.0 || v == 1.0);
    let dummy_ang = p[3..].iter().all(|&v| v == 0.0 || v == 90.0);
    let sg = cs.space_group.replace(' ', "");
    if dummy_len && dummy_ang && (sg.is_empty() || sg == "P1") {
        return None;
    }
    let info = space_group_info(&cs.space_group, &p)?;
    let uc = UnitCell::new(p).and_then(|u| u.averaged(cryst1_rotations(&cs.space_group, &p)?))?;
    let ops: Option<Vec<SymOp>> = info.ops.iter().map(|s| SymOp::parse(s)).collect();
    Some((uc, ops?))
}

#[cfg(test)]
mod symop_tests {
    use super::*;

    #[test]
    fn parses_operators() {
        let op = SymOp::parse("-x+y,-x,z+2/3").unwrap();
        assert_eq!(op.r, [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(op.t, [0.0, 0.0, 2.0 / 3.0]);
        assert!(SymOp::parse("x,y,z").unwrap().is_identity());
        let f = SymOp::parse("-x+1/2,y,-z").unwrap().apply(v3(0.1, 0.2, 0.3));
        assert!((f.x - 0.4).abs() < 1e-12 && (f.y - 0.2).abs() < 1e-12 && (f.z + 0.3).abs() < 1e-12);
    }
}
