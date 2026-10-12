use crate::geom::Vec3;
use crate::probe::SpatialGrid;
use rustc_hash::FxHashMap;

pub(crate) enum Query {
    Dense(SpatialGrid),
    Sparse {
        cells: FxHashMap<[i64; 3], Vec<u32>>,
        positions: Vec<Vec3>,
    },
}

fn bucket(position: Vec3) -> [i64; 3] {
    [
        (position.x / 4.0).floor() as i64,
        (position.y / 4.0).floor() as i64,
        (position.z / 4.0).floor() as i64,
    ]
}

impl Query {
    pub fn new(positions: &[Vec3]) -> Result<Self, String> {
        let mut lower = [f64::INFINITY; 3];
        let mut upper = [f64::NEG_INFINITY; 3];
        for p in positions {
            for (axis, value) in p.as_array().into_iter().enumerate() {
                if !value.is_finite() || value.abs() > i64::MAX as f64 {
                    return Err("Probe coordinates exceed the spatial-index range".into());
                }
                lower[axis] = lower[axis].min(value);
                upper[axis] = upper[axis].max(value);
            }
        }
        let cells = (0..3).try_fold(1usize, |count, axis| {
            let side = ((upper[axis] - lower[axis]) / 4.0).floor() as usize;
            count.checked_mul(side.saturating_add(2))
        });
        if cells.is_some_and(|count| count <= positions.len().saturating_mul(32).max(4096)) {
            let points: Vec<_> = positions
                .iter()
                .enumerate()
                .map(|(i, &p)| (i as u32, p))
                .collect();
            return Ok(Self::Dense(SpatialGrid::new(&points, 4.0)));
        }
        let mut cells: FxHashMap<[i64; 3], Vec<u32>> = FxHashMap::default();
        for (i, &position) in positions.iter().enumerate() {
            cells.entry(bucket(position)).or_default().push(i as u32);
        }
        Ok(Self::Sparse {
            cells,
            positions: positions.to_vec(),
        })
    }

    pub fn for_each_within<F: FnMut(u32, Vec3, f64)>(
        &self,
        p: Vec3,
        min: f64,
        max: f64,
        mut visit: F,
    ) {
        match self {
            Self::Dense(grid) => grid.for_each_within(p, min, max, visit),
            Self::Sparse { cells, positions } => {
                let low = bucket(p - crate::geom::v3(max, max, max));
                let high = bucket(p + crate::geom::v3(max, max, max));
                for z in low[2]..=high[2] {
                    for y in low[1]..=high[1] {
                        for x in low[0]..=high[0] {
                            if let Some(items) = cells.get(&[x, y, z]) {
                                for &i in items {
                                    let position = positions[i as usize];
                                    let distance = p.dist_sq(position);
                                    if distance >= min * min && distance <= max * max {
                                        visit(i, position, distance);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
