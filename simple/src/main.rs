// This is a simple non-relativisetic 1D MHD code.
//
// Author: Brayden JoHantgen
// Last Update: 5/27/2026

use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

pub mod math_func;

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: f64 = 800.0;
const DISCON: f64 = 0.5;
const ADIABATIC: f64 = 2.0;
const DR: f64 = 1.0;// / CELL_NUM;
const T_FINAL: f64 = 0.401;
const CHECK_INTERVAL: f64 = 0.025;
const CFL: f64 = 0.8;
const BX: f64 = 0.75;

////////////////////
// Usage Functions
////////////////////

/// Input:
/// Output:
/// Description:
fn init_prim() -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive = Vec::new();
    for i in 0..(CELL_NUM as u32) {
        if i < ((CELL_NUM * DISCON) as u32) {
            init_primitive.push((1.0, 1.0, 0.0, 0.0, 0.0, BX, 1.0, 0.0));
        } else {
            init_primitive.push((0.1, 0.125, 0.0, 0.0, 0.0, BX, -1.0, 0.0));
        }
    }
    init_primitive
}

/// Input:
/// Output:
/// Description:
fn cons_vec_from_prim(prims: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
    let mut cons_vec = Vec::new();
    for i in 0..(CELL_NUM as u32) {
        let index: usize = (i).try_into().unwrap();
        cons_vec.push(math_func::prim_to_cons(prims[index], a_index));
    }
    cons_vec
}

/// Input:
/// Output:
/// Description:
fn prim_vec_from_cons(cons: Vec<(f64, f64, f64, f64, f64, f64, f64)>, a_index: f64, bx: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut prims_vec = Vec::new();
    for i in 0..(CELL_NUM as u32) {
        let index: usize = (i).try_into().unwrap();
        prims_vec.push(math_func::cons_to_prim(cons[index], a_index, bx));
    }
    prims_vec
}

/// Input:
/// Output:
/// Description:
fn hll_flux(prim_l: (f64, f64, f64, f64, f64, f64, f64, f64), prim_r: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
    let plus_l = math_func::max_eigen(prim_l.clone(), a_index);
    let minus_l = math_func::min_eigen(prim_l.clone(), a_index);
    let u_l = math_func::prim_to_cons(prim_l.clone(), a_index);
    let f_l = math_func::flux(prim_l.clone(), a_index);

    let plus_r = math_func::max_eigen(prim_r.clone(), a_index);
    let minus_r = math_func::min_eigen(prim_r.clone(), a_index);
    let u_r = math_func::prim_to_cons(prim_r.clone(), a_index);
    let f_r = math_func::flux(prim_r.clone(), a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6);
    hll
}

/// Input:
/// Output:
/// Description:
fn godonov(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
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

/// Input:
/// Output:
/// Description:
fn l_function(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64, f64, f64, f64, f64)>, dt: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
    let go_vec = godonov(prims_vec, ADIABATIC);
    let mut new_cons_vec = Vec::new();
    for i in 0..(CELL_NUM as u32) {
        let index_1: usize = (i).try_into().unwrap();
        let index_2: usize = (i+1).try_into().unwrap();
        let new_0 = cons_vec[index_1].0 - (go_vec[index_2].0 - go_vec[index_1].0) * dt / DR;
        let new_1 = cons_vec[index_1].1 - (go_vec[index_2].1 - go_vec[index_1].1) * dt / DR;
        let new_2 = cons_vec[index_1].2 - (go_vec[index_2].2 - go_vec[index_1].2) * dt / DR;
        let new_3 = cons_vec[index_1].3 - (go_vec[index_2].3 - go_vec[index_1].3) * dt / DR;
        let new_4 = cons_vec[index_1].4 - (go_vec[index_2].4 - go_vec[index_1].4) * dt / DR;
        let new_5 = cons_vec[index_1].5 - (go_vec[index_2].5 - go_vec[index_1].5) * dt / DR;
        let new_6 = cons_vec[index_1].6 - (go_vec[index_2].6 - go_vec[index_1].6) * dt / DR;
        let new_fill = (new_0, new_1, new_2, new_3, new_4, new_5, new_6);
        new_cons_vec.push(new_fill);
    }
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
fn write_checkpoint(prims: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, t: f64, check_count: i8) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);

    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));

    let num_fill = format!("{} {}", "cell_num:".to_string(), &(CELL_NUM.to_string()+&" ".to_string()));

    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut vx_string = "vx: ".to_string();
    let mut vy_string = "vy: ".to_string();
    let mut vz_string = "vz: ".to_string();
    let mut bx_string = "Bx: ".to_string();
    let mut by_string = "By: ".to_string();
    let mut bz_string = "Bz: ".to_string();

    for i in prims {
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
    let mut t: f64 = 0.0;
    let mut t_checkpoint = CHECK_INTERVAL;
    let mut check_count: i8 = 0;

    let initial_primitives = init_prim();
    let mut conserved_vec = cons_vec_from_prim(initial_primitives.clone(), ADIABATIC);

    while t < T_FINAL {
        let primitives = prim_vec_from_cons(conserved_vec.clone(), ADIABATIC, BX);
        let conserve = cons_vec_from_prim(primitives.clone(), ADIABATIC);
        
        let mut dt = 1.0;
        for i in 0..((CELL_NUM - 1.0) as u8) {
            let index_1: usize = (i).try_into().unwrap();
            let index_2: usize = (i+1).try_into().unwrap();
            let dt_check = math_func::compute_time_step(primitives[index_1], primitives[index_2], ADIABATIC, DR);
            if dt_check < dt {
                dt = dt_check;
                }
            }
        dt = CFL * dt;
        
         conserved_vec = l_function(primitives.clone(), conserve, dt);

        if t >= t_checkpoint {
            let _ = write_checkpoint(primitives.clone(), t, check_count);
            t_checkpoint += CHECK_INTERVAL;
            check_count += 1;
            }
        t += dt;
    }
}