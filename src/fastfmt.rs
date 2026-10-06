//! Fixed-point formatting of floats, the same text as `format!("{:w.p}")`
//! (round half to even on the exact binary value, a sign for every negative
//! value, -0 included), without the general formatting machinery.

/// Append `x` as `format!("{:>width$.prec$}")` would write it.
#[inline]
pub fn push_fixed(out: &mut String, x: f64, prec: u32, width: usize) {
    if !x.is_finite() || prec > 9 {
        use std::fmt::Write;
        let _ = write!(out, "{:>w$.p$}", x, w = width, p = prec as usize);
        return;
    }
    let bits = x.to_bits();
    let neg = bits >> 63 == 1;
    let exp = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1u64 << 52) - 1);
    let (m, e) = if exp == 0 { (frac, -1074) } else { (frac | (1u64 << 52), exp - 1075) };
    let scale = 10u128.pow(prec);
    // n = round_half_even(m * 2^e * 10^prec), exactly
    let n: u128 = if e >= 0 {
        if e > 40 {
            use std::fmt::Write;
            let _ = write!(out, "{:>w$.p$}", x, w = width, p = prec as usize);
            return;
        }
        ((m as u128) * scale) << e
    } else {
        let sh = (-e) as u32;
        let p = (m as u128) * scale;
        if sh >= 127 {
            0
        } else {
            let q = p >> sh;
            let r = p & ((1u128 << sh) - 1);
            let half = 1u128 << (sh - 1);
            if r > half || (r == half && q & 1 == 1) { q + 1 } else { q }
        }
    };
    let (int, fr) = (n / scale, n % scale);
    let mut buf = [0u8; 64];
    let mut k = buf.len();
    // fraction digits
    let mut f = fr;
    for _ in 0..prec {
        k -= 1;
        buf[k] = b'0' + (f % 10) as u8;
        f /= 10;
    }
    if prec > 0 {
        k -= 1;
        buf[k] = b'.';
    }
    let mut i = int;
    loop {
        k -= 1;
        buf[k] = b'0' + (i % 10) as u8;
        i /= 10;
        if i == 0 {
            break;
        }
    }
    if neg {
        k -= 1;
        buf[k] = b'-';
    }
    let len = buf.len() - k;
    for _ in len..width {
        out.push(' ');
    }
    // the digits are ASCII
    out.push_str(std::str::from_utf8(&buf[k..]).unwrap());
}

#[cfg(test)]
mod tests {
    use super::push_fixed;

    fn check(x: f64, prec: u32, width: usize) {
        let mut a = String::new();
        push_fixed(&mut a, x, prec, width);
        let b = format!("{:>w$.p$}", x, w = width, p = prec as usize);
        assert_eq!(a, b, "x = {:e} ({:#x}), prec {}, width {}", x, x.to_bits(), prec, width);
    }

    #[test]
    fn matches_format_on_random_ties_and_edges() {
        let mut s: u64 = 0x9e3779b97f4a7c15;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let specials = [0.0, -0.0, 1e-300, -1e-300, 5e-324, 0.5, -0.5, 1.5, 2.5, 0.0625, -0.0625, 0.1875, 1.0625, 0.0005, 0.00049999, 99999.995, 1e7, -1e7, 1e8, 123456.789];
        for prec in [0u32, 1, 2, 3, 5] {
            for width in [0usize, 6, 7, 8] {
                for &x in &specials {
                    check(x, prec, width);
                }
                for _ in 0..200_000 {
                    let r = next();
                    // coordinates and B factors: a random value, a value on a binary grid (ties)
                    let a = ((r >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2.0e4;
                    check(a, prec, width);
                    let t = ((r % 4_000_000) as f64 - 2_000_000.0) / 2f64.powi((next() % 12) as i32);
                    check(t, prec, width);
                    check(f64::from_bits(r) % 1e6, prec, width);
                }
            }
        }
    }
}
