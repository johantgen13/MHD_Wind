// This is a simple non-relativisetic MHD code.
//
// Author: Brayden JoHantgen
// Last Update: 5/22/2026

pub mod math_func;

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: f64 = 10.0;
const DISCON: f64 = 0.5;
const ADIABATIC: f64 = 2.0;

////////////////////
// Usage Functions
////////////////////
fn init_prim() -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive = Vec::new();
    for i in 0..(CELL_NUM as u8) {
        if i < ((CELL_NUM * DISCON) as u8) {
            init_primitive.push((1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0));
        } else {
            init_primitive.push((0.1, 0.125, 0.0, 0.0, 0.0, 0.5, -1.0, 0.0));
        }
    }
    init_primitive
}

///////////////
// Simulation
///////////////
fn main() {
    let init_prims = init_prim();
    let con_0 = math_func::prim_to_cons(init_prims[0], ADIABATIC);
    let test_val = math_func::flux(init_prims[0], ADIABATIC);
    println!("{:?}", test_val);
}
