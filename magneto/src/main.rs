// This is the rust file containing the 1D relativistic mhd simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/1/2026

//use std::fs;
//use std::env;

/////////////////
/// Dataclasses
/////////////////
//struct Config {
//    num_zones: u8,
//    dr: f64,
//    discon: f64,
//    p: f64,
//    rho: f64,
//    v: f64,
//}

//struct Physics {
//}

//struct Driver {
//}

///////////////////
// Math Functions
///////////////////
fn specific_energy_gas(prim: Vec<f64>, a_index: f64) -> f64 {
    let e = prim[0] / ((a_index - 1.0) * prim[1]);
    e
}

fn sound_speed(prim: Vec<f64>, a_index: f64) -> f64 {
    let cs = ((a_index * prim[0]) / prim[1]).sqrt();
    cs
}

fn tot_energy_density(prim: Vec<f64>, a_index: f64) -> f64 {
    let e = specific_energy_gas(prim.clone(), a_index);
    let energy = prim[1] * e + (0.5) * prim[1] * (prim[2] * prim[2]);
    energy
}

/////////////////////
/// Usage Functions
/////////////////////
fn init_prim() -> Vec<Vec<f64>> {
    let mut init_primitive = Vec::new();
    for i in 0..20 {
        if i < 10 {
            let vec_fill = vec![1.0, 1.0, 0.0];
            init_primitive.push(vec_fill);
        } else {
            let vec_fill = vec![0.125, 0.1, 0.0];
            init_primitive.push(vec_fill);
        }
    }
    init_primitive
}

////////////////
/// Simulation
////////////////
fn main() {
    let prim_vec = init_prim(); 
//    let test = tot_energy_density(prim_vec, 1.4);
    println!("{:?}", prim_vec);  
}