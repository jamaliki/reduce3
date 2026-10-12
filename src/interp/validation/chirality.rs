use super::*;

/// CCTBX's target is the unit-cell volume from three bond lengths and angles.
pub(super) fn target(comp: &Comp, chir: &crate::monlib::CompChir) -> Option<f64> {
    let bond = |name: &str| {
        comp.bonds
            .iter()
            .find(|b| {
                (b.a1 == chir.centre && b.a2 == name) || (b.a2 == chir.centre && b.a1 == name)
            })
            .and_then(|b| b.value_dist)
    };
    let angle = |a: &str, b: &str| {
        comp.angles
            .iter()
            .find(|p| p.a2 == chir.centre && ((p.a1 == a && p.a3 == b) || (p.a1 == b && p.a3 == a)))
            .and_then(|p| p.value)
            .map(f64::to_radians)
    };
    let lengths = [bond(&chir.a[0])?, bond(&chir.a[1])?, bond(&chir.a[2])?];
    let cosines = [
        angle(&chir.a[1], &chir.a[2])?.cos(),
        angle(&chir.a[0], &chir.a[2])?.cos(),
        angle(&chir.a[0], &chir.a[1])?.cos(),
    ];
    let [a, b, c] = cosines;
    let determinant = 1.0 + 2.0 * a * b * c - a * a - b * b - c * c;
    if determinant <= 0.0 || !determinant.is_finite() {
        return None;
    }
    let sign = if chir.volume_sign.to_ascii_lowercase().starts_with("neg") {
        -1.0
    } else {
        1.0
    };
    Some(sign * lengths.iter().product::<f64>() * determinant.sqrt())
}

pub(super) fn pseudo(residue: &str, name: &str) -> bool {
    name == "P"
        || matches!(residue, "SF4" | "F3S")
        || (residue == "VAL" && name == "CB")
        || (residue == "LEU" && name == "CG")
}

pub(super) fn link_target(
    previous: &Comp,
    current: &Comp,
    link: &crate::monlib::ChemLink,
    chir: &crate::monlib::LinkChir,
) -> Option<f64> {
    let key = |side: u8, name: &str| format!("{side}:{name}");
    let mut merged = Comp::default();
    for (side, comp) in [(1, previous), (2, current)] {
        merged
            .bonds
            .extend(comp.bonds.iter().cloned().map(|mut bond| {
                bond.a1 = key(side, &bond.a1);
                bond.a2 = key(side, &bond.a2);
                bond
            }));
        merged
            .angles
            .extend(comp.angles.iter().cloned().map(|mut angle| {
                angle.a1 = key(side, &angle.a1);
                angle.a2 = key(side, &angle.a2);
                angle.a3 = key(side, &angle.a3);
                angle
            }));
    }
    merged.bonds.extend(
        link.bonds
            .iter()
            .filter(|bond| bond.c1 != bond.c2)
            .map(|bond| crate::monlib::CompBond {
                a1: key(bond.c1, &bond.a1),
                a2: key(bond.c2, &bond.a2),
                value_dist: bond.value_dist,
                ..Default::default()
            }),
    );
    merged.angles.extend(
        link.angles
            .iter()
            .filter(|angle| angle.c.iter().any(|&side| side != angle.c[0]))
            .map(|angle| crate::monlib::CompAngle {
                a1: key(angle.c[0], &angle.a[0]),
                a2: key(angle.c[1], &angle.a[1]),
                a3: key(angle.c[2], &angle.a[2]),
                value: angle.value,
                ..Default::default()
            }),
    );
    target(
        &merged,
        &crate::monlib::CompChir {
            centre: key(chir.c[0], &chir.a[0]),
            a: std::array::from_fn(|i| key(chir.c[i + 1], &chir.a[i + 1])),
            volume_sign: chir.volume_sign.clone(),
            ..Default::default()
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::monlib::{CompAngle, CompBond, CompChir};

    #[test]
    fn target_uses_dictionary_geometry_and_sign() {
        let mut comp = Comp::default();
        for a in ["A", "B", "C"] {
            comp.bonds.push(CompBond {
                a1: "O".into(),
                a2: a.into(),
                value_dist: Some(2.),
                ..Default::default()
            });
        }
        for (a, b) in [("A", "B"), ("A", "C"), ("B", "C")] {
            comp.angles.push(CompAngle {
                a1: a.into(),
                a2: "O".into(),
                a3: b.into(),
                value: Some(90.),
                ..Default::default()
            });
        }
        let mut chir = CompChir {
            centre: "O".into(),
            a: ["A".into(), "B".into(), "C".into()],
            volume_sign: "positive".into(),
            ..Default::default()
        };
        assert_eq!(target(&comp, &chir), Some(8.));
        chir.volume_sign = "negative".into();
        assert_eq!(target(&comp, &chir), Some(-8.));
    }
}
