//! Small 3D vector toolkit matching the scitbx conventions used by Reduce2.
//!
//! * Dihedral angles follow the IUPAC sign convention (same as
//!   `scitbx.math.dihedral_angle`).
//! * Rotations about an axis are right-handed (same as
//!   `scitbx::math::rotate_point_around_axis`).

use std::ops::{Add, AddAssign, Div, Index, Mul, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[inline(always)]
pub const fn v3(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3 { x, y, z }
}

impl Vec3 {
    pub const ZERO: Vec3 = v3(0.0, 0.0, 0.0);

    #[inline(always)]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    #[inline(always)]
    pub fn cross(self, o: Vec3) -> Vec3 {
        v3(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    #[inline(always)]
    pub fn length_sq(self) -> f64 {
        self.dot(self)
    }
    #[inline(always)]
    pub fn length(self) -> f64 {
        self.length_sq().sqrt()
    }
    /// Unit vector; returns the zero vector for a zero-length input instead of NaNs.
    #[inline(always)]
    pub fn normalize(self) -> Vec3 {
        let l = self.length();
        if l > 0.0 { self / l } else { Vec3::ZERO }
    }
    #[inline(always)]
    pub fn dist(self, o: Vec3) -> f64 {
        (self - o).length()
    }
    #[inline(always)]
    pub fn dist_sq(self, o: Vec3) -> f64 {
        (self - o).length_sq()
    }
    /// Angle between two vectors in radians, clamped so rounding never yields NaN.
    pub fn angle(self, o: Vec3) -> f64 {
        let d = self.length() * o.length();
        if d == 0.0 {
            return 0.0;
        }
        (self.dot(o) / d).clamp(-1.0, 1.0).acos()
    }
    #[inline(always)]
    pub fn as_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    #[inline(always)]
    pub fn from_array(a: [f64; 3]) -> Vec3 {
        v3(a[0], a[1], a[2])
    }
    /// A vector perpendicular to this one (not normalized), exactly as
    /// `scitbx.matrix.col.ortho()` chooses it.
    pub fn ortho(self) -> Vec3 {
        let (a, b, c) = (self.x.abs(), self.y.abs(), self.z.abs());
        if c <= a && c <= b {
            v3(-self.y, self.x, 0.0)
        } else if b <= a && b <= c {
            v3(-self.z, 0.0, self.x)
        } else {
            v3(0.0, -self.z, self.y)
        }
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline(always)]
    fn add(self, o: Vec3) -> Vec3 {
        v3(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl AddAssign for Vec3 {
    #[inline(always)]
    fn add_assign(&mut self, o: Vec3) {
        self.x += o.x;
        self.y += o.y;
        self.z += o.z;
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    #[inline(always)]
    fn sub(self, o: Vec3) -> Vec3 {
        v3(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl SubAssign for Vec3 {
    #[inline(always)]
    fn sub_assign(&mut self, o: Vec3) {
        self.x -= o.x;
        self.y -= o.y;
        self.z -= o.z;
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline(always)]
    fn mul(self, s: f64) -> Vec3 {
        v3(self.x * s, self.y * s, self.z * s)
    }
}
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline(always)]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}
impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline(always)]
    fn div(self, s: f64) -> Vec3 {
        v3(self.x / s, self.y / s, self.z / s)
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    #[inline(always)]
    fn neg(self) -> Vec3 {
        v3(-self.x, -self.y, -self.z)
    }
}
impl Index<usize> for Vec3 {
    type Output = f64;
    #[inline(always)]
    fn index(&self, i: usize) -> &f64 {
        match i {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }
}

/// Dihedral angle in degrees with the IUPAC sign convention, or `None` when the
/// geometry is degenerate (collinear atoms), like `scitbx.math.dihedral_angle`.
pub fn dihedral_deg(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> Option<f64> {
    // scitbx divides by pi/180 (not multiply by 180/pi); keep the same rounding.
    dihedral_rad(p0, p1, p2, p3).map(|r| r / PI_180)
}

pub fn dihedral_rad(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> Option<f64> {
    // Same arithmetic as scitbx::math::dihedral::angle().
    let d_01 = p0 - p1;
    let d_21 = p2 - p1;
    let d_23 = p2 - p3;
    let n_0121 = d_01.cross(d_21);
    let n0 = n_0121.length_sq();
    let n_2123 = d_21.cross(d_23);
    let n1 = n_2123.length_sq();
    if n0 == 0.0 || n1 == 0.0 {
        return None;
    }
    let cos_angle = (n_0121.dot(n_2123) / (n0 * n1).sqrt()).clamp(-1.0, 1.0);
    let mut result = cos_angle.acos();
    if d_21.dot(n_0121.cross(n_2123)) < 0.0 {
        result = -result;
    }
    Some(result)
}

/// `scitbx::math::dihedral::angle()` as the arm64 cctbx build evaluates it
/// (fused dot and cross products). Python callers of
/// `scitbx.math.dihedral_angle` get this one.
pub fn dihedral_rad_cpp(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> Option<f64> {
    use crate::riding::cx;
    let d_01 = p0 - p1;
    let d_21 = p2 - p1;
    let d_23 = p2 - p3;
    let n_0121 = cx::cross(d_01, d_21);
    let n0 = cx::dot(n_0121, n_0121);
    let n_2123 = cx::cross(d_21, d_23);
    let n1 = cx::dot(n_2123, n_2123);
    if n0 == 0.0 || n1 == 0.0 {
        return None;
    }
    let cos_angle = (cx::dot(n_0121, n_2123) / (n0 * n1).sqrt()).clamp(-1.0, 1.0);
    let mut result = cos_angle.acos();
    if cx::dot(d_21, cx::cross(n_0121, n_2123)) < 0.0 {
        result = -result;
    }
    Some(result)
}

/// `scitbx::math::rotate_point_around_axis` with identical arithmetic.
#[inline]
pub fn rotate_point_around_axis(a1: Vec3, a2: Vec3, point: Vec3, angle_rad: f64) -> Vec3 {
    let (xa, ya, za) = (a1.x, a1.y, a1.z);
    let xl = a2.x - xa;
    let yl = a2.y - ya;
    let zl = a2.z - za;
    let xlsq = xl * xl;
    let ylsq = yl * yl;
    let zlsq = zl * zl;
    let dlsq = xlsq + ylsq + zlsq;
    let dl = dlsq.sqrt();
    let ca = angle_rad.cos();
    let dsa = angle_rad.sin() / dl;
    let oca = (1.0 - ca) / dlsq;
    let xlylo = xl * yl * oca;
    let xlzlo = xl * zl * oca;
    let ylzlo = yl * zl * oca;
    let xma = point.x - xa;
    let yma = point.y - ya;
    let zma = point.z - za;
    v3(
        xma * (xlsq * oca + ca) + yma * (xlylo - zl * dsa) + zma * (xlzlo + yl * dsa) + xa,
        xma * (xlylo + zl * dsa) + yma * (ylsq * oca + ca) + zma * (ylzlo - xl * dsa) + ya,
        xma * (xlzlo - yl * dsa) + yma * (ylzlo + xl * dsa) + zma * (zlsq * oca + ca) + za,
    )
}

/// Reduce2's `RotatePointDegreesAroundAxisDir`: axis through `origin` with
/// direction `dir` (second axis point is origin + dir), right-handed.
#[inline]
pub fn rotate_deg_axis_dir(origin: Vec3, dir: Vec3, point: Vec3, degrees: f64) -> Vec3 {
    rotate_point_around_axis(origin, origin + dir, point, degrees * (PI_180))
}

const PI_180: f64 = std::f64::consts::PI / 180.0;

/// `scitbx.matrix.col.rotate_around_origin` (radians).
#[inline]
pub fn rotate_around_origin_scitbx(x: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let n = axis / axis.length();
    let (c, s) = (angle.cos(), angle.sin());
    x * c + n * n.dot(x) * (1.0 - c) + n.cross(x) * s
}

/// Rotate `point` about the axis through `origin` with direction `dir`
/// (need not be normalized) by `degrees`, right-handed.
#[inline]
pub fn rotate_about_axis_dir(origin: Vec3, dir: Vec3, point: Vec3, degrees: f64) -> Vec3 {
    let (s, c) = degrees.to_radians().sin_cos();
    rotate_about_axis_sc(origin, dir.normalize(), point, s, c)
}

/// Rodrigues rotation with a precomputed unit axis and sin/cos of the angle.
#[inline(always)]
pub fn rotate_about_axis_sc(origin: Vec3, unit: Vec3, point: Vec3, s: f64, c: f64) -> Vec3 {
    let v = point - origin;
    let rotated = v * c + unit.cross(v) * s + unit * (unit.dot(v) * (1.0 - c));
    origin + rotated
}

/// Rotate a vector about an axis through the origin (radians), like
/// `scitbx.matrix.col.rotate_around_origin`.
#[inline]
pub fn rotate_around_origin(v: Vec3, axis: Vec3, radians: f64) -> Vec3 {
    let (s, c) = radians.sin_cos();
    rotate_about_axis_sc(Vec3::ZERO, axis.normalize(), v, s, c)
}

/// Clamped arc-cosine in degrees, guarding against values just outside [-1, 1]
/// produced by rounding (the original code could raise a math domain error).
#[inline]
pub fn acos_deg_clamped(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).acos().to_degrees()
}

/// Determinant of a 3x3 matrix given in row-major order.
#[inline]
pub fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dihedral_sign_matches_scitbx() {
        let d = dihedral_deg(v3(1., 0., 0.), v3(0., 0., 0.), v3(0., 0., 1.), v3(0., 1., 1.)).unwrap();
        assert!((d - 90.0).abs() < 1e-12);
        let d = dihedral_deg(v3(1., 0., 0.), v3(0., 0., 0.), v3(0., 0., 1.), v3(1., 0., 1.)).unwrap();
        assert!(d.abs() < 1e-12);
    }

    #[test]
    fn rotation_is_right_handed() {
        let p = rotate_about_axis_dir(v3(0., 0., 0.), v3(0., 0., 1.), v3(1., 0., 0.), 90.0);
        assert!((p - v3(0., 1., 0.)).length() < 1e-12);
        let p = rotate_about_axis_dir(v3(1., 2., 3.), v3(0., 0., 2.), v3(2., 2., 3.), 90.0);
        assert!((p - v3(1., 3., 3.)).length() < 1e-12);
    }
}
