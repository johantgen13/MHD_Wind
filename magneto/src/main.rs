// This is the rust file containing the 1D relativistic mhd simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/5/2026

//use std::fs;
//use std::env;

#![allow(dead_code)]

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: f64 = 10.0;
const DISCON: f64 = 0.5;
const DR: f64 = 1.0 / CELL_NUM;
const ADIABATIC: f64 = 1.4;
const HIGH_ORDER: bool = false;
const T_FINAL: f64 = 0.401;

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

/// Input: 
///     tup: A tuple that contains three elements of type f64.
/// Output: 
///     The maximum value of the elements in the tuple, a float.
/// Description: 
///     This function converts the tuple to an array and then iterates over the elements of the array. 
///     If an array/tuple element is larger than zero it is saved and compared to the rest of the 
///     elements. The largest element is returned. This function will fail if all elements are negative.
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

fn godonov(prims_vec: Vec<(f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64)> {
    let mut go_vec = Vec::new();
    go_vec.push(hll_flux(prims_vec[0], prims_vec[1], a_index));
    for i in 0..((CELL_NUM - 1.0) as u32) {
        let index_1: usize = (i).try_into().unwrap();
        let index_2: usize = (i+1).try_into().unwrap();
        let go_fill = hll_flux(prims_vec[index_1], prims_vec[index_2], a_index);
        go_vec.push(go_fill);
    }
    let index_a: usize = ((CELL_NUM - 2.0) as u32).try_into().unwrap();
    let index_b: usize = ((CELL_NUM - 1.0) as u32).try_into().unwrap();
    go_vec.push(hll_flux(prims_vec[index_a], prims_vec[index_b], a_index));
    go_vec
}

fn compute_time_step(prim_l: (f64, f64, f64), prim_r: (f64, f64, f64), a_index: f64) -> f64 {
    let plus_l = p_eigen(prim_l.clone(), a_index);
    let minus_l = m_eigen(prim_l.clone(), a_index);

    let plus_r = p_eigen(prim_r.clone(), a_index);
    let minus_r = m_eigen(prim_r.clone(), a_index);

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let mut dt: f64 = 0.0;

    if a_minus > a_plus {
        dt += DR / a_minus;
    } else {
        dt += DR / a_plus;
    }    
    dt
}

fn l_function(prims_vec: Vec<(f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64)>, dt: f64) -> Vec<(f64, f64, f64)> {
    let go_vec = godonov(prims_vec, ADIABATIC);
    let mut new_cons_vec = Vec::new();
    for i in 0..(CELL_NUM as u8) {
        let index_1: usize = (i).try_into().unwrap();
        let index_2: usize = (i+1).try_into().unwrap();
        let new_0 = cons_vec[index_1].0 - (go_vec[index_2].0 - go_vec[index_1].0) * dt / DR;
        let new_1 = cons_vec[index_1].1 - (go_vec[index_2].1 - go_vec[index_1].1) * dt / DR;
        let new_2 = cons_vec[index_1].2 - (go_vec[index_2].2 - go_vec[index_1].2) * dt / DR;
        let new_fill = (new_0, new_1, new_2);
        new_cons_vec.push(new_fill);
    }
    new_cons_vec
}

////////////////////
// Usage Functions
////////////////////
fn init_prim() -> Vec<(f64, f64, f64)> {
    let mut init_primitive = Vec::new(); 
    for i in 0..(CELL_NUM as u8) {
        if i < ((CELL_NUM * DISCON) as u8) {
            init_primitive.push((1.0, 1.0, 0.0));
        } else {
            init_primitive.push((0.125, 0.1, 0.0));
        }
    }
    init_primitive
}

fn cons_vec_from_prim(prims: Vec<(f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64)> {
    let mut cons_vec = Vec::new();
    for i in 0..(CELL_NUM as u8) {
        let index: usize = (i).try_into().unwrap();
        cons_vec.push(prim_to_cons(prims[index], a_index));
    }
    cons_vec
}

fn prim_vec_from_cons(cons: Vec<(f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64)> {
    let mut prims_vec = Vec::new();
    for i in 0..(CELL_NUM as u8) {
        let index: usize = (i).try_into().unwrap();
        prims_vec.push(cons_to_prim(cons[index], a_index));
    }
    prims_vec
}

//fn rk_time_step(cons_vec: Vec<(f64, f64, f64)>, dt: f64) -> Vec<(f64, f64, f64)> {

//}

///////////////
// Simulation
///////////////
fn main() {
    let mut t: f64 = 0.0;
    let initial_primitives = init_prim();
    let mut conserved_vec = cons_vec_from_prim(initial_primitives.clone(), ADIABATIC);

    while t < T_FINAL {
        let primitives = prim_vec_from_cons(conserved_vec, ADIABATIC);
        let conserve = cons_vec_from_prim(primitives.clone(), ADIABATIC);

        let mut dt = 1.0;
        for i in 0..((CELL_NUM - 1.0) as u8) {
            let index_1: usize = (i).try_into().unwrap();
            let index_2: usize = (i+1).try_into().unwrap();
            let dt_check = compute_time_step(primitives[index_1], primitives[index_2], ADIABATIC);
            if dt_check < dt {
                dt = dt_check;
            }
        }
        conserved_vec = l_function(primitives.clone(), conserve, dt);

        println!("First Cell: {:?}", primitives[0]);
        println!("Middle Cell: {:?}", primitives[5]);
        println!("Last Cell: {:?}", primitives[9]);

        t += dt;
    }
    
} 