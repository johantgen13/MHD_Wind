// This is a simple non-relativisetic 1D MHD code.
//
// Author: Brayden JoHantgen
// Last Update: 6/7/2026

use std::fs;
//use std::str;
use std::io::{BufWriter, Write, BufRead, BufReader, Read};
use std::path::Path;
use std::time::Instant;
//use std::collections::HashMap;

pub mod math_func;

////////////
// Structs
////////////
struct Physics {
    adiabatic_index: f64,
    p: (f64, f64),
    rho: (f64, f64),
    vx: (f64, f64),
    vy: (f64, f64),
    vz: (f64, f64),
    bx: (f64, f64),
    by: (f64, f64),
    bz: (f64, f64),
}

struct Driver {
    cfl: f64,
    tfinal: f64,
    checkpoint: f64,
    num_zones_x: usize,
    num_zones_y: usize,
    discontinuity: f64,

    dimensionality: String,
    plm: bool,
    grid_type: String,
}

////////////////////
// Usage Functions
////////////////////

// Input:
// Output:
// Description:
//fn read_config(file_path: String) -> Result<(), Box<dyn std::error::Error>> { //Physics, Driver {
//    let input_file_path = Path::new(&file_path);
//    let file = fs::read_to_string(input_file_path).map_err(|e| e.to_string())?;
    
//    for l in file.lines() {
//        let line = l.trim();
//    }

//    if line.is_empty() {
//        continue;
//    }

//    let parts: Vec<_> = line.split_whitespace().collect();
//    println!("{:?}", parts);
//    Ok(())   
//}

/// Input:
/// Output:
/// Description:
fn init_prims_1D(phys: &Physics, zones: usize, discon: f64, plm: bool) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
    if plm == false {
        for i in 0..zones {
            if i < (((zones as f64) * discon) as usize) {
                init_primitive.push((phys.p.0, phys.rho.0, phys.vx.0, phys.vy.0, phys.vz.0, phys.bx.0, phys.by.0, phys.bz.0));
            } else {
                init_primitive.push((phys.p.1, phys.rho.1, phys.vx.1, phys.vy.1, phys.vz.1, phys.bx.1, phys.by.1, phys.bz.1));
            }
        }
    } else {
        for i in 0..(zones+2) {
            if i < (((zones as f64) * discon + 1.0) as usize) {
                init_primitive.push((phys.p.0, phys.rho.0, phys.vx.0, phys.vy.0, phys.vz.0, phys.bx.0, phys.by.0, phys.bz.0));
            } else {
                init_primitive.push((phys.p.1, phys.rho.1, phys.vx.1, phys.vy.1, phys.vz.1, phys.bx.1, phys.by.1, phys.bz.1));
            }
        }
    }
    init_primitive
}

/// Input:
/// Output:
/// Description:
fn cons_vec_from_prim(prims: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut cons_vec = Vec::new();
    for i in 0..(prims.len() as usize) {
        cons_vec.push(math_func::prim_to_cons(prims[i], a_index));
    }
    cons_vec
}

/// Input:
/// Output:
/// Description:
fn prim_vec_from_cons(cons: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut prims_vec = Vec::new();
    for i in 0..(cons.len() as usize) {
        prims_vec.push(math_func::cons_to_prim(cons[i], a_index));
    }
    prims_vec
}

/// Input:
/// Output:
/// Description:
fn hll_flux_x(prim_l: (f64, f64, f64, f64, f64, f64, f64, f64), prim_r: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let plus_l = math_func::max_eigen(prim_l.clone(), a_index);
    let minus_l = math_func::min_eigen(prim_l.clone(), a_index);
    let u_l = math_func::prim_to_cons(prim_l.clone(), a_index);
    let f_l = math_func::flux_x(prim_l.clone(), a_index);

    let plus_r = math_func::max_eigen(prim_r.clone(), a_index);
    let minus_r = math_func::min_eigen(prim_r.clone(), a_index);
    let u_r = math_func::prim_to_cons(prim_r.clone(), a_index);
    let f_r = math_func::flux_x(prim_r.clone(), a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll_7 = ((a_plus * f_l.7) + (a_minus * f_r.7) - (a_plus * a_minus * (u_r.7 - u_l.7))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6, hll_7);
    hll
}

/// Input:
/// Output:
/// Description:
fn hll_flux_y(prim_l: (f64, f64, f64, f64, f64, f64, f64, f64), prim_r: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let plus_l = math_func::max_eigen(prim_l.clone(), a_index);
    let minus_l = math_func::min_eigen(prim_l.clone(), a_index);
    let u_l = math_func::prim_to_cons(prim_l.clone(), a_index);
    let f_l = math_func::flux_y(prim_l.clone(), a_index);

    let plus_r = math_func::max_eigen(prim_r.clone(), a_index);
    let minus_r = math_func::min_eigen(prim_r.clone(), a_index);
    let u_r = math_func::prim_to_cons(prim_r.clone(), a_index);
    let f_r = math_func::flux_y(prim_r.clone(), a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll_7 = ((a_plus * f_l.7) + (a_minus * f_r.7) - (a_plus * a_minus * (u_r.7 - u_l.7))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6, hll_7);
    hll
}

/// Input:
/// Output:
/// Description:
fn hll_flux_x_plm(prim_1: (f64, f64, f64, f64, f64, f64, f64, f64), prim_2: (f64, f64, f64, f64, f64, f64, f64, f64), prim_3: (f64, f64, f64, f64, f64, f64, f64, f64), prim_4: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let p_l = math_func::left_reconstruction(prim_1.0, prim_2.0, prim_3.0);
    let rho_l = math_func::left_reconstruction(prim_1.1, prim_2.1, prim_3.1);
    let vx_l = math_func::left_reconstruction(prim_1.2, prim_2.2, prim_3.2);
    let vy_l = math_func::left_reconstruction(prim_1.3, prim_2.3, prim_3.3);
    let vz_l = math_func::left_reconstruction(prim_1.4, prim_2.4, prim_3.4);
    let bx_l = math_func::left_reconstruction(prim_1.5, prim_2.5, prim_3.5);
    let by_l = math_func::left_reconstruction(prim_1.6, prim_2.6, prim_3.6);
    let bz_l = math_func::left_reconstruction(prim_1.7, prim_2.7, prim_3.7);

    let p_r = math_func::right_reconstruction(prim_2.0, prim_3.0, prim_4.0);
    let rho_r = math_func::right_reconstruction(prim_2.1, prim_3.1, prim_4.1);
    let vx_r = math_func::right_reconstruction(prim_2.2, prim_3.2, prim_4.2);
    let vy_r = math_func::right_reconstruction(prim_2.3, prim_3.3, prim_4.3);
    let vz_r = math_func::right_reconstruction(prim_2.4, prim_3.4, prim_4.4);
    let bx_r = math_func::right_reconstruction(prim_2.5, prim_3.5, prim_4.5);
    let by_r = math_func::right_reconstruction(prim_2.6, prim_3.6, prim_4.6);
    let bz_r = math_func::right_reconstruction(prim_2.7, prim_3.7, prim_4.7);

    let prim_l = (p_l, rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l);
    let prim_r = (p_r, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r);

    let plus_l = math_func::max_eigen(prim_l.clone(), a_index);
    let minus_l = math_func::min_eigen(prim_l.clone(), a_index);
    let u_l = math_func::prim_to_cons(prim_l.clone(), a_index);
    let f_l = math_func::flux_x(prim_l.clone(), a_index);

    let plus_r = math_func::max_eigen(prim_r.clone(), a_index);
    let minus_r = math_func::min_eigen(prim_r.clone(), a_index);
    let u_r = math_func::prim_to_cons(prim_r.clone(), a_index);
    let f_r = math_func::flux_x(prim_r.clone(), a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll_7 = ((a_plus * f_l.7) + (a_minus * f_r.7) - (a_plus * a_minus * (u_r.7 - u_l.7))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6, hll_7);
    hll
}

/// Input:
/// Output:
/// Description:
fn hll_flux_y_plm(prim_1: (f64, f64, f64, f64, f64, f64, f64, f64), prim_2: (f64, f64, f64, f64, f64, f64, f64, f64), prim_3: (f64, f64, f64, f64, f64, f64, f64, f64), prim_4: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let p_l = math_func::left_reconstruction(prim_1.0, prim_2.0, prim_3.0);
    let rho_l = math_func::left_reconstruction(prim_1.1, prim_2.1, prim_3.1);
    let vx_l = math_func::left_reconstruction(prim_1.2, prim_2.2, prim_3.2);
    let vy_l = math_func::left_reconstruction(prim_1.3, prim_2.3, prim_3.3);
    let vz_l = math_func::left_reconstruction(prim_1.4, prim_2.4, prim_3.4);
    let bx_l = math_func::left_reconstruction(prim_1.5, prim_2.5, prim_3.5);
    let by_l = math_func::left_reconstruction(prim_1.6, prim_2.6, prim_3.6);
    let bz_l = math_func::left_reconstruction(prim_1.7, prim_2.7, prim_3.7);

    let p_r = math_func::right_reconstruction(prim_2.0, prim_3.0, prim_4.0);
    let rho_r = math_func::right_reconstruction(prim_2.1, prim_3.1, prim_4.1);
    let vx_r = math_func::right_reconstruction(prim_2.2, prim_3.2, prim_4.2);
    let vy_r = math_func::right_reconstruction(prim_2.3, prim_3.3, prim_4.3);
    let vz_r = math_func::right_reconstruction(prim_2.4, prim_3.4, prim_4.4);
    let bx_r = math_func::right_reconstruction(prim_2.5, prim_3.5, prim_4.5);
    let by_r = math_func::right_reconstruction(prim_2.6, prim_3.6, prim_4.6);
    let bz_r = math_func::right_reconstruction(prim_2.7, prim_3.7, prim_4.7);

    let prim_l = (p_l, rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l);
    let prim_r = (p_r, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r);

    let plus_l = math_func::max_eigen(prim_l.clone(), a_index);
    let minus_l = math_func::min_eigen(prim_l.clone(), a_index);
    let u_l = math_func::prim_to_cons(prim_l.clone(), a_index);
    let f_l = math_func::flux_y(prim_l.clone(), a_index);

    let plus_r = math_func::max_eigen(prim_r.clone(), a_index);
    let minus_r = math_func::min_eigen(prim_r.clone(), a_index);
    let u_r = math_func::prim_to_cons(prim_r.clone(), a_index);
    let f_r = math_func::flux_y(prim_r.clone(), a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll_7 = ((a_plus * f_l.7) + (a_minus * f_r.7) - (a_plus * a_minus * (u_r.7 - u_l.7))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6, hll_7);
    hll
}

/// Input:
/// Output:
/// Description:
fn godonov_x(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64, plm: bool) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut go_vec = Vec::new();
    if plm == false {
        go_vec.push(hll_flux_x(prims_vec[0], prims_vec[1], a_index));
        for i in 1..((prims_vec.len()) as usize) {
            let go_fill = hll_flux_x(prims_vec[i-1], prims_vec[i], a_index);
            go_vec.push(go_fill);
        }
        go_vec.push(hll_flux_x(prims_vec[(prims_vec.len() - 2) as usize], prims_vec[(prims_vec.len() - 1) as usize], a_index));
    } else {
        go_vec.push(hll_flux_x_plm(prims_vec[0], prims_vec[1], prims_vec[2], prims_vec[3], a_index));
        for i in 1..((prims_vec.len() - 2) as usize) {
            let go_fill = hll_flux_x_plm(prims_vec[i-1], prims_vec[i], prims_vec[i+1], prims_vec[i+2], a_index);
            go_vec.push(go_fill);
        }
        go_vec.push(hll_flux_x_plm(prims_vec[(prims_vec.len() - 4) as usize], prims_vec[(prims_vec.len() - 3) as usize], prims_vec[(prims_vec.len() - 2) as usize], prims_vec[(prims_vec.len() - 1) as usize], a_index));
    }
    go_vec
}


/// Input:
/// Output:
/// Description:
fn l_function_1D(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64, plm: bool) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut new_cons_vec = Vec::new();
    let go_vec = godonov_x(prims_vec.clone(), a_index, plm);
    if plm == false {
        let dx: f64 = 1.0 / (prims_vec.len() as f64);
        for i in 0..((prims_vec.len()) as usize) {
            let new_0 = - (go_vec[i+1].0 - go_vec[i].0) / dx;
            let new_1 = - (go_vec[i+1].1 - go_vec[i].1) / dx;
            let new_2 = - (go_vec[i+1].2 - go_vec[i].2) / dx;
            let new_3 = - (go_vec[i+1].3 - go_vec[i].3) / dx;
            let new_4 = - (go_vec[i+1].4 - go_vec[i].4) / dx;
            let new_5 = - (go_vec[i+1].5 - go_vec[i].5) / dx;
            let new_6 = - (go_vec[i+1].6 - go_vec[i].6) / dx;
            let new_7 = - (go_vec[i+1].7 - go_vec[i].7) / dx;
            let new_fill = (new_0, new_1, new_2, new_3, new_4, new_5, new_6, new_7);
            new_cons_vec.push(new_fill);
        }
    } else {
        let dx: f64 = 1.0 / ((prims_vec.len() - 2) as f64);
        new_cons_vec.push(cons_vec[0]);
        for i in 1..((prims_vec.len() + 1) as usize) {
            let new_0 = - (go_vec[i].0 - go_vec[i-1].0) / dx;
            let new_1 = - (go_vec[i].1 - go_vec[i-1].1) / dx;
            let new_2 = - (go_vec[i].2 - go_vec[i-1].2) / dx;
            let new_3 = - (go_vec[i].3 - go_vec[i-1].3) / dx;
            let new_4 = - (go_vec[i].4 - go_vec[i-1].4) / dx;
            let new_5 = - (go_vec[i].5 - go_vec[i-1].5) / dx;
            let new_6 = - (go_vec[i].6 - go_vec[i-1].6) / dx;
            let new_7 = - (go_vec[i].7 - go_vec[i-1].7) / dx;
            let new_fill = (new_0, new_1, new_2, new_3, new_4, new_5, new_6, new_7);
            new_cons_vec.push(new_fill);
        }
        new_cons_vec.push(cons_vec[(prims_vec.len() + 1) as usize])
    }
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
fn euler_timestep(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64, dt: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let l_cons = l_function_1D(prims_vec.clone(), cons_vec.clone(), a_index, false);
    let mut new_cons_vec = Vec::new();
    for i in 0..(prims_vec.len() as usize) {
        let fill_0 = cons_vec[i].0 + dt * l_cons[i].0;
        let fill_1 = cons_vec[i].1 + dt * l_cons[i].1;
        let fill_2 = cons_vec[i].2 + dt * l_cons[i].2;
        let fill_3 = cons_vec[i].3 + dt * l_cons[i].3;
        let fill_4 = cons_vec[i].4 + dt * l_cons[i].4;
        let fill_5 = cons_vec[i].5 + dt * l_cons[i].5;
        let fill_6 = cons_vec[i].6 + dt * l_cons[i].6;
        let fill_7 = cons_vec[i].7 + dt * l_cons[i].7;
        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6, fill_7);
        new_cons_vec.push(fill);
    }
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
//fn rk4_step(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64, f64, f64, f64, f64)>, dt: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
//    let l_cons = l_function(prims_vec.clone(), cons_vec.clone());
//    let mut cons_1 = Vec::new();
//    cons_1.push(cons_vec[0]);
//    for i in 1..((CELL_NUM + 1.0) as u64) {
//        let index_a: usize = (i).try_into().unwrap();
//        let fill_0 = cons_vec[index_a].0 + dt * l_cons[index_a].0;
//        let fill_1 = cons_vec[index_a].1 + dt * l_cons[index_a].1;
//        let fill_2 = cons_vec[index_a].2 + dt * l_cons[index_a].2;
//        let fill_3 = cons_vec[index_a].3 + dt * l_cons[index_a].3;
//        let fill_4 = cons_vec[index_a].4 + dt * l_cons[index_a].4;
//        let fill_5 = cons_vec[index_a].5 + dt * l_cons[index_a].5;
//        let fill_6 = cons_vec[index_a].6 + dt * l_cons[index_a].6;
//        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
//        cons_1.push(fill);
//    }
//    let b = (CELL_NUM + 1.0) as u64;
//    let index_b: usize = (b).try_into().unwrap();
//    cons_1.push(cons_vec[index_b]);

//    let prims_1 = prim_vec_from_cons(cons_1.clone(), ADIABATIC, BX);
//    let l_cons_1 = l_function(prims_1.clone(), cons_1.clone());
//    let mut cons_2 = Vec::new();
//    cons_2.push(cons_vec[0]);
//    for i in 1..((CELL_NUM + 1.0) as u64) {
//        let index_c: usize = (i).try_into().unwrap();
//        let fill_0 = 0.75 * cons_vec[index_c].0 + 0.25 * cons_1[index_c].0 + 0.25 * dt * l_cons_1[index_c].0;
//        let fill_1 = 0.75 * cons_vec[index_c].1 + 0.25 * cons_1[index_c].1 + 0.25 * dt * l_cons_1[index_c].1;
//        let fill_2 = 0.75 * cons_vec[index_c].2 + 0.25 * cons_1[index_c].2 + 0.25 * dt * l_cons_1[index_c].2;
//        let fill_3 = 0.75 * cons_vec[index_c].3 + 0.25 * cons_1[index_c].3 + 0.25 * dt * l_cons_1[index_c].3;
//        let fill_4 = 0.75 * cons_vec[index_c].4 + 0.25 * cons_1[index_c].4 + 0.25 * dt * l_cons_1[index_c].4;
//        let fill_5 = 0.75 * cons_vec[index_c].5 + 0.25 * cons_1[index_c].5 + 0.25 * dt * l_cons_1[index_c].5;
//        let fill_6 = 0.75 * cons_vec[index_c].6 + 0.25 * cons_1[index_c].6 + 0.25 * dt * l_cons_1[index_c].6;
//        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
//        cons_2.push(fill);
//    }
//    cons_2.push(cons_vec[index_b]);

//    let prims_2 = prim_vec_from_cons(cons_1.clone(), ADIABATIC, BX);
//    let l_cons_2 = l_function(prims_2.clone(), cons_2.clone());
//    let mut new_cons = Vec::new();
//    new_cons.push(cons_vec[0]);
//    for i in 1..((CELL_NUM + 1.0) as u64) {
//        let index_d: usize = (i).try_into().unwrap();
//        let fill_0 = 0.33 * cons_vec[index_d].0 + 0.67 * cons_2[index_d].0 + 0.67 * dt * l_cons_2[index_d].0;
//        let fill_1 = 0.33 * cons_vec[index_d].1 + 0.67 * cons_2[index_d].1 + 0.67 * dt * l_cons_2[index_d].1;
//        let fill_2 = 0.33 * cons_vec[index_d].2 + 0.67 * cons_2[index_d].2 + 0.67 * dt * l_cons_2[index_d].2;
//        let fill_3 = 0.33 * cons_vec[index_d].3 + 0.67 * cons_2[index_d].3 + 0.67 * dt * l_cons_2[index_d].3;
//        let fill_4 = 0.33 * cons_vec[index_d].4 + 0.67 * cons_2[index_d].4 + 0.67 * dt * l_cons_2[index_d].4;
//        let fill_5 = 0.33 * cons_vec[index_d].5 + 0.67 * cons_2[index_d].5 + 0.67 * dt * l_cons_2[index_d].5;
//        let fill_6 = 0.33 * cons_vec[index_d].6 + 0.67 * cons_2[index_d].6 + 0.67 * dt * l_cons_2[index_d].6;
//        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
//        new_cons.push(fill);
//    }
//    new_cons.push(cons_vec[index_b]);
//    new_cons
//}

/// Input:
/// Output:
/// Description:
fn write_checkpoint(prims: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, t: f64, check_count: i8, plm: bool) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);

    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));

    let num_fill = format!("{} {}", "cell_num:".to_string(), &(prims.len().to_string()+&" ".to_string()));

    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut vx_string = "vx: ".to_string();
    let mut vy_string = "vy: ".to_string();
    let mut vz_string = "vz: ".to_string();
    let mut bx_string = "Bx: ".to_string();
    let mut by_string = "By: ".to_string();
    let mut bz_string = "Bz: ".to_string();

    let mut prims_fill = Vec::new();
    if plm == true{
        for i in 1..((prims.len()+1) as usize) {
            prims_fill.push(prims[i]);
        }
    } else {
        for i in 0..((prims.len()-1) as usize) {
            prims_fill.push(prims[i]);
        }
    }
    

    for i in prims_fill {
        let p_fill = i.0.to_string() + &" ".to_string();
        p_string.push_str(&p_fill);

        let rho_fill = i.1.to_string() + &" ".to_string();
        rho_string.push_str(&rho_fill);

        let vx_fill = i.2.to_string() + &" ".to_string();
        vx_string.push_str(&vx_fill);

        let vy_fill = i.3.to_string() + &" ".to_string();
        vy_string.push_str(&vy_fill);

        let vz_fill = i.4.to_string() + &" ".to_string();
        vz_string.push_str(&vz_fill);

        let bx_fill = i.5.to_string() + &" ".to_string();
        bx_string.push_str(&bx_fill);

        let by_fill = i.6.to_string() + &" ".to_string();
        by_string.push_str(&by_fill);

        let bz_fill = i.7.to_string() + &" ".to_string();
        bz_string.push_str(&bz_fill);
    }
 
    let file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(output_file_path)?;

    let mut writer = BufWriter::new(file);

    writeln!(writer, "{}", num_fill)?;
    writeln!(writer, "{}", t_fill)?;
    writeln!(writer, "{}", p_string)?;
    writeln!(writer, "{}", rho_string)?;
    writeln!(writer, "{}", vx_string)?;
    writeln!(writer, "{}", vy_string)?;
    writeln!(writer, "{}", vz_string)?;
    writeln!(writer, "{}", bx_string)?;
    writeln!(writer, "{}", by_string)?;
    writeln!(writer, "{}", bz_string)?;
    writer.flush()?;

    Ok(())
}

///////////////
// Simulation
///////////////
fn main() {
    let phys = Physics{adiabatic_index: 2.0, p: (1.0, 0.1), rho: (1.0, 0.125), vx: (0.0, 0.0), vy: (0.0, 0.0), vz: (0.0, 0.0), bx: (0.75, 0.75), by: (1.0, -1.0), bz: (0.0, 0.0)};
    let drive = Driver{cfl: 0.8, tfinal: 0.401, checkpoint: 0.0125, num_zones_x: 800, num_zones_y: 10, discontinuity: 0.5, dimensionality: "1D".to_string(), plm: false, grid_type: "Cartesian".to_string()};

    let before = Instant::now();

    //let path = "/Users/toogan13/Desktop/MHD_Wind/magneto_config.txt";
    //let config = read_config(path.to_string());
    //println!("{:?}", config);

    let mut t: f64 = 0.0;
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i8 = 0;

    let initial_primitives = init_prims_1D(&phys, drive.num_zones_x, drive.discontinuity, drive.plm);
    let mut conserved_vec = cons_vec_from_prim(initial_primitives.clone(), phys.adiabatic_index);
    //let test_val = l_function_1D(initial_primitives, conserved_vec, phys.adiabatic_index, drive.plm);
    //println!("{:?}", test_val);

    while t < drive.tfinal {
        let primitives = prim_vec_from_cons(conserved_vec.clone(), phys.adiabatic_index);
        let conserve = cons_vec_from_prim(primitives.clone(), phys.adiabatic_index);
        
        let mut dt = 1.0;
        for i in 0..((initial_primitives.len() - 1) as usize) {
            let dt_check = math_func::compute_time_step(primitives[i], primitives[i+1], phys.adiabatic_index, 1.0 / (initial_primitives.len() as f64));
            if dt_check < dt {
                dt = dt_check; 
                }
            }
        dt = drive.cfl * dt;
        
        if drive.plm == true {
            //conserved_vec = rk4_step(primitives.clone(), conserve, dt);
        } else {
            conserved_vec = euler_timestep(primitives.clone(), conserve, phys.adiabatic_index, dt);
        }

        if t >= t_checkpoint {
            let _ = write_checkpoint(primitives.clone(), t, check_count, drive.plm);
            t_checkpoint += drive.checkpoint;
            check_count += 1;
            }
        t += dt;
        time_step_count += 1.0;
    }
    let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
    let performance = (time_step_count * (initial_primitives.len() as f64)) / runtime;
    println!("The number of timesteps: {:?}", time_step_count);
    println!("The runtime in seconds: {:?}", runtime);
    println!("The performance in zones per second: {:.2?}", performance);
}