use reduce3::cell::SymOp;
use reduce3::geom::v3;

#[test]
fn inverse_preserves_translations_and_nonorthogonal_fractional_rotations() {
    for operation in ["x-1,y,z", "-y,x-y,z+1/3", "-x,y+1/2,-z"] {
        let forward = SymOp::parse(operation).unwrap();
        let reverse = forward.inverse().unwrap();
        for input in [v3(0., 0., 0.), v3(0.13, -2.3, 4.5)] {
            let output = reverse.apply(forward.apply(input));
            assert!(input.dist(output) < 1e-12);
        }
    }
}

#[test]
fn singular_operator_cannot_supply_a_reverse_symmetry_bond() {
    let operation = SymOp {
        r: [[0.; 3]; 3],
        t: [0.; 3],
    };
    assert!(operation.inverse().is_none());
}
