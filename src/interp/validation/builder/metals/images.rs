use super::*;
use crate::cell::{SymOp, crystal_operators};

pub(super) fn for_each(
    st: &Structure,
    flat: &FlatAtoms,
    anchors: &[usize],
    visit: &mut impl FnMut(usize, Vec3, SymOp),
) -> Result<(), String> {
    if let Some((cell, ops)) = crystal_operators(st) {
        let fractional: Vec<_> = flat.pos.iter().map(|&p| cell.fractionalize(p)).collect();
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for position in &fractional {
            for (axis, value) in position.as_array().into_iter().enumerate() {
                lo[axis] = lo[axis].min(value);
                hi[axis] = hi[axis].max(value);
            }
        }
        let margin: Vec<_> = cell
            .frac
            .chunks_exact(3)
            .map(|row| 3.5 * row.iter().map(|v| v * v).sum::<f64>().sqrt())
            .collect();
        for &anchor in anchors {
            for op in &ops {
                let position = op.apply(fractional[anchor]).as_array();
                let (low, high) = shift_bounds(lo, hi, &margin, position)?;
                for x in low[0]..=high[0] {
                    for y in low[1]..=high[1] {
                        for z in low[2]..=high[2] {
                            let mut image = *op;
                            for (t, shift) in image.t.iter_mut().zip([x, y, z]) {
                                *t += shift as f64;
                            }
                            if !image.is_identity() {
                                visit(
                                    anchor,
                                    cell.orthogonalize(image.apply(fractional[anchor])),
                                    image,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn shift_bounds(
    lo: [f64; 3],
    hi: [f64; 3],
    margin: &[f64],
    position: [f64; 3],
) -> Result<([i64; 3], [i64; 3]), String> {
    let mut low = [0; 3];
    let mut high = [0; 3];
    for axis in 0..3 {
        let lower = (lo[axis] - margin[axis] - position[axis]).ceil();
        let upper = (hi[axis] + margin[axis] - position[axis]).floor();
        if [lower, upper]
            .iter()
            .any(|value| !value.is_finite() || value.abs() >= (i64::MAX / 4) as f64)
        {
            return Err("Metal symmetry shifts exceed the integer range".into());
        }
        low[axis] = lower as i64;
        high[axis] = upper as i64;
    }
    let images = (0..3).try_fold(1u64, |n, axis| {
        n.checked_mul((high[axis] - low[axis] + 1).max(0) as u64)
    });
    if images.is_none_or(|n| n > 1_000_000) {
        return Err("Metal symmetry search exceeds its image bound".into());
    }
    Ok((low, high))
}

#[cfg(test)]
mod tests {
    use super::shift_bounds;

    #[test]
    fn extreme_cells_cannot_overflow_shift_bounds() {
        assert!(shift_bounds([0.; 3], [1.; 3], &[f64::INFINITY; 3], [0.; 3]).is_err());
        assert!(shift_bounds([0.; 3], [1e30; 3], &[1.; 3], [0.; 3]).is_err());
        assert!(shift_bounds([0.; 3], [1000.; 3], &[1.; 3], [0.; 3]).is_err());
        assert_eq!(
            shift_bounds([0.; 3], [1.; 3], &[1.; 3], [0.; 3]).unwrap(),
            ([-1; 3], [2; 3])
        );
    }
}
