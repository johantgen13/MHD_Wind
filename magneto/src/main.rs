// This is the rust file containing the 1D relativistic mhd simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/4/2026

//use std::fs;
//use std::env;

#![allow(dead_code)]

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: f64 = 20.0;
const DISCON: f64 = 0.5;
const DR: f64 = 1.0 / CELL_NUM;

////////////////
// Dataclasses
////////////////
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
fn tuple_max(tup: (f64, f64, f64)) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut max_check: f64 = 0.0;
    for i in 0..3 {
        if arr[i] > max_check {
            max_check = arr[i];
        }
    }
    max_check
}

fn specific_energy_gas(prim: (f64, f64, f64), a_index: f64) -> f64 {
    let e = prim.0 / ((a_index - 1.0) * prim.1);
    e
}

fn sound_speed(prim: (f64, f64, f64), a_index: f64) -> f64 {
    let cs = ((a_index * prim.0) / prim.1).sqrt();
    cs
}

fn tot_energy_density(prim: (f64, f64, f64), a_index: f64) -> f64 {
    let e = specific_energy_gas(prim.clone(), a_index);
    let energy = prim.1 * e + (0.5) * prim.1 * (prim.2 * prim.2);
    energy
}

fn inverse_energy_density(energy: f64, rho: f64, v: f64, a_index: f64) -> f64 {
    let pressure = (a_index - 1.0) * (energy - 0.5 * rho * v * v);
    pressure
}

fn prim_to_cons(prim: (f64, f64, f64), a_index: f64) -> (f64, f64, f64) {
    let energy_density = tot_energy_density(prim.clone(), a_index);
    let cons = (prim.1, (prim.1 * prim.2), energy_density);
    cons
}

fn cons_to_prim(cons: (f64, f64, f64), a_index: f64) -> (f64, f64, f64) {
    let v = cons.1 / cons.0;
    let pressure = inverse_energy_density(cons.2, cons.0, v, a_index);
    let prim = (pressure, cons.0, v);
    prim
}

fn flux(prim: (f64, f64, f64), a_index: f64) -> (f64, f64, f64) {
    let energy = tot_energy_density(prim.clone(), a_index);
    let cell_flux = ((prim.1 * prim.2), ((prim.1 * prim.2 * prim.2) + prim.0), ((energy + prim.0) * prim.2));
    cell_flux
}

fn p_eigen(prim: (f64, f64, f64), a_index: f64) -> f64 {
    let cs = sound_speed(prim.clone(), a_index);
    let max_eigen = prim.2 + cs;
    max_eigen
}

fn m_eigen(prim: (f64, f64, f64), a_index: f64) -> f64 {
    let cs = sound_speed(prim.clone(), a_index);
    let min_eigen = prim.2 - cs;
    min_eigen
}

fn hll_flux(prim_l: (f64, f64, f64), prim_r: (f64, f64, f64), a_index: f64) -> (f64, f64, f64) {
    let plus_l = p_eigen(prim_l.clone(), a_index);
    let minus_l = m_eigen(prim_l.clone(), a_index);
    let u_l = prim_to_cons(prim_l.clone(), a_index);
    let f_l = flux(prim_l.clone(), a_index);

    let plus_r = p_eigen(prim_r.clone(), a_index);
    let minus_r = m_eigen(prim_r.clone(), a_index);
    let u_r = prim_to_cons(prim_r.clone(), a_index);
    let f_r = flux(prim_r.clone(), a_index);

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2);
    hll
}

////////////////////
// Usage Functions
////////////////////
fn init_prim() -> Vec<(f64, f64, f64)> {
    let mut init_primitive = Vec::new(); 
    for i in 0..(CELL_NUM as u8) {
        if i < ((CELL_NUM * DISCON) as u8) {
            init_primitive.push((1.0, 1.0, 0.0))
        } else {
            init_primitive.push((0.125, 0.1, 0.0))
        }
    }
    init_primitive
}

///////////////
// Simulation
///////////////
fn main() {
    let prim_vec = init_prim(); 
    let hll_test = hll_flux(prim_vec[9], prim_vec[10], 1.4);
    println!("{:?}", hll_test);  
}