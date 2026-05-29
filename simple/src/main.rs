// This is a simple non-relativisetic 1D MHD code.
//
// Author: Brayden JoHantgen
// Last Update: 5/29/2026

use std::fs;
//use std::str;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

pub mod math_func;

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: usize = 100;
const DISCON: f64 = 0.5;
const ADIABATIC: f64 = 2.0;
const DR: f64 = 1.0 / (CELL_NUM as f64);
const T_FINAL: f64 = 1.26;
const CHECK_INTERVAL: f64 = 0.025;
const CFL: f64 = 0.8;
const BX: f64 = 0.75;

////////////////////
// Usage Functions
////////////////////

/// Input:
/// Output:
/// Description:
//fn read_config(file_path: String) -> Vec<u8> {
//    let input_file_path = Path::new(&file_path);
//    let data_bytes = fs::read(input_file_path);
//    let byte_vec = data_bytes.expect("Something went wrong!");
//    let index_vec = math_func::vector_index(byte_vec, 10);
//    index_vec
    //for i in byte_vec {
    //    println!("{:?}", i);
    //}
//}

/// Input:
/// Output:
/// Description:
fn init_prim() -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive = Vec::new();
    for i in 0..(CELL_NUM + 2) {
        if i < (((CELL_NUM as f64) * DISCON + 1.0) as usize) {
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
fn cons_vec_from_prim(prims: &Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
    let mut cons_vec = Vec::new();
    for i in 0..(CELL_NUM+2) {
        cons_vec.push(math_func::prim_to_cons(prims[i], a_index));
    }
    cons_vec
}

/// Input:
/// Output:
/// Description:
fn prim_vec_from_cons(cons: &Vec<(f64, f64, f64, f64, f64, f64, f64)>, a_index: f64, bx: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut prims_vec = Vec::new();
    for i in 0..(CELL_NUM+2) {
        prims_vec.push(math_func::cons_to_prim(cons[i], a_index, bx));
    }
    prims_vec
}

/// Input:
/// Output:
/// Description:
fn hll_flux(prim_1: (f64, f64, f64, f64, f64, f64, f64, f64), prim_2: (f64, f64, f64, f64, f64, f64, f64, f64), prim_3: (f64, f64, f64, f64, f64, f64, f64, f64), prim_4: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
    let p_l = math_func::left_reconstruction(prim_1.0, prim_2.0, prim_3.0);
    let rho_l = math_func::left_reconstruction(prim_1.1, prim_2.1, prim_3.1);
    let vx_l = math_func::left_reconstruction(prim_1.2, prim_2.2, prim_3.2);
    let vy_l = math_func::left_reconstruction(prim_1.3, prim_2.3, prim_3.3);
    let vz_l = math_func::left_reconstruction(prim_1.4, prim_2.4, prim_3.4);
    let bx_l = prim_2.5;
    let by_l = math_func::left_reconstruction(prim_1.6, prim_2.6, prim_3.6);
    let bz_l = math_func::left_reconstruction(prim_1.7, prim_2.7, prim_3.7);

    let p_r = math_func::right_reconstruction(prim_2.0, prim_3.0, prim_4.0);
    let rho_r = math_func::right_reconstruction(prim_2.1, prim_3.1, prim_4.1);
    let vx_r = math_func::right_reconstruction(prim_2.2, prim_3.2, prim_4.2);
    let vy_r = math_func::right_reconstruction(prim_2.3, prim_3.3, prim_4.3);
    let vz_r = math_func::right_reconstruction(prim_2.4, prim_3.4, prim_4.4);
    let bx_r = prim_3.5;
    let by_r = math_func::right_reconstruction(prim_2.6, prim_3.6, prim_4.6);
    let bz_r = math_func::right_reconstruction(prim_2.7, prim_3.7, prim_4.7);

    let prim_l = (p_l, rho_l, vx_l, vy_l, vz_l, bx_l, by_l, bz_l);
    let prim_r = (p_r, rho_r, vx_r, vy_r, vz_r, bx_r, by_r, bz_r);

    let plus_l = math_func::max_eigen(prim_l, a_index);
    let minus_l = math_func::min_eigen(prim_l, a_index);
    let u_l = math_func::prim_to_cons(prim_l, a_index);
    let f_l = math_func::flux(prim_l, a_index);

    let plus_r = math_func::max_eigen(prim_r, a_index);
    let minus_r = math_func::min_eigen(prim_r, a_index);
    let u_r = math_func::prim_to_cons(prim_r, a_index);
    let f_r = math_func::flux(prim_r, a_index);

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
    go_vec.push(hll_flux(prims_vec[0], prims_vec[1], prims_vec[2], prims_vec[3], a_index));
    for i in 1..CELL_NUM {
        let go_fill = hll_flux(prims_vec[i-1], prims_vec[i], prims_vec[i+1], prims_vec[i+2], a_index);
        go_vec.push(go_fill);
    }
    go_vec.push(hll_flux(prims_vec[CELL_NUM-2], prims_vec[CELL_NUM-1], prims_vec[CELL_NUM], prims_vec[CELL_NUM+1], a_index));
    go_vec
}

/// Input:
/// Output:
/// Description:
fn l_function(prims_vec: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: &Vec<(f64, f64, f64, f64, f64, f64, f64)>) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
    let go_vec = godonov(prims_vec, ADIABATIC);
    let mut new_cons_vec = Vec::new();
    new_cons_vec.push(cons_vec[0]);
    for i in 1..(CELL_NUM+1) {
        let new_0 = - (go_vec[i].0 - go_vec[i-1].0) / DR;
        let new_1 = - (go_vec[i].1 - go_vec[i-1].1) / DR;
        let new_2 = - (go_vec[i].2 - go_vec[i-1].2) / DR;
        let new_3 = - (go_vec[i].3 - go_vec[i-1].3) / DR;
        let new_4 = - (go_vec[i].4 - go_vec[i-1].4) / DR;
        let new_5 = - (go_vec[i].5 - go_vec[i-1].5) / DR;
        let new_6 = - (go_vec[i].6 - go_vec[i-1].6) / DR;
        let new_fill = (new_0, new_1, new_2, new_3, new_4, new_5, new_6);
        new_cons_vec.push(new_fill);
    }
    new_cons_vec.push(cons_vec[CELL_NUM+1]);
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
fn rk4_step(prims_vec: &Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, cons_vec: Vec<(f64, f64, f64, f64, f64, f64, f64)>, dt: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64)> {
    let l_cons = l_function(prims_vec.to_vec(), &cons_vec);
    let mut cons_1 = Vec::new();
    cons_1.push(cons_vec[0]);
    for i in 1..(CELL_NUM+1) {
        let fill_0 = cons_vec[i].0 + dt * l_cons[i].0;
        let fill_1 = cons_vec[i].1 + dt * l_cons[i].1;
        let fill_2 = cons_vec[i].2 + dt * l_cons[i].2;
        let fill_3 = cons_vec[i].3 + dt * l_cons[i].3;
        let fill_4 = cons_vec[i].4 + dt * l_cons[i].4;
        let fill_5 = cons_vec[i].5 + dt * l_cons[i].5;
        let fill_6 = cons_vec[i].6 + dt * l_cons[i].6;
        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
        cons_1.push(fill);
    }
    cons_1.push(cons_vec[CELL_NUM+1]);

    let prims_1 = prim_vec_from_cons(&cons_1, ADIABATIC, BX);
    let l_cons_1 = l_function(prims_1, &cons_1);
    let mut cons_2 = Vec::new();
    cons_2.push(cons_vec[0]);
    for i in 1..(CELL_NUM+1) {
        let fill_0 = 0.75 * cons_vec[i].0 + 0.25 * cons_1[i].0 + 0.25 * dt * l_cons_1[i].0;
        let fill_1 = 0.75 * cons_vec[i].1 + 0.25 * cons_1[i].1 + 0.25 * dt * l_cons_1[i].1;
        let fill_2 = 0.75 * cons_vec[i].2 + 0.25 * cons_1[i].2 + 0.25 * dt * l_cons_1[i].2;
        let fill_3 = 0.75 * cons_vec[i].3 + 0.25 * cons_1[i].3 + 0.25 * dt * l_cons_1[i].3;
        let fill_4 = 0.75 * cons_vec[i].4 + 0.25 * cons_1[i].4 + 0.25 * dt * l_cons_1[i].4;
        let fill_5 = 0.75 * cons_vec[i].5 + 0.25 * cons_1[i].5 + 0.25 * dt * l_cons_1[i].5;
        let fill_6 = 0.75 * cons_vec[i].6 + 0.25 * cons_1[i].6 + 0.25 * dt * l_cons_1[i].6;
        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
        cons_2.push(fill);
    }
    cons_2.push(cons_vec[CELL_NUM+1]);

    let prims_2 = prim_vec_from_cons(&cons_1, ADIABATIC, BX);
    let l_cons_2 = l_function(prims_2, &cons_2);
    let mut new_cons = Vec::new();
    new_cons.push(cons_vec[0]);
    for i in 1..(CELL_NUM+1) {
        let fill_0 = 0.33 * cons_vec[i].0 + 0.67 * cons_2[i].0 + 0.67 * dt * l_cons_2[i].0;
        let fill_1 = 0.33 * cons_vec[i].1 + 0.67 * cons_2[i].1 + 0.67 * dt * l_cons_2[i].1;
        let fill_2 = 0.33 * cons_vec[i].2 + 0.67 * cons_2[i].2 + 0.67 * dt * l_cons_2[i].2;
        let fill_3 = 0.33 * cons_vec[i].3 + 0.67 * cons_2[i].3 + 0.67 * dt * l_cons_2[i].3;
        let fill_4 = 0.33 * cons_vec[i].4 + 0.67 * cons_2[i].4 + 0.67 * dt * l_cons_2[i].4;
        let fill_5 = 0.33 * cons_vec[i].5 + 0.67 * cons_2[i].5 + 0.67 * dt * l_cons_2[i].5;
        let fill_6 = 0.33 * cons_vec[i].6 + 0.67 * cons_2[i].6 + 0.67 * dt * l_cons_2[i].6;
        let fill = (fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6);
        new_cons.push(fill);
    }
    new_cons.push(cons_vec[CELL_NUM+1]);
    new_cons
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

    let mut prims_fill = Vec::new();
    for i in 1..(CELL_NUM+1) {
        prims_fill.push(prims[i]);
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
    let before = Instant::now();

    //let path = "/Users/toogan13/Desktop/MHD_Wind/magneto_config.txt";
    //let bytes = read_config(path.to_string());
    //println!("{:?}", bytes);

    let mut t: f64 = 0.0;
    let mut t_checkpoint = CHECK_INTERVAL;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i8 = 0;

    let initial_primitives = init_prim();
    let mut conserved_vec = cons_vec_from_prim(&initial_primitives, ADIABATIC);
    //println!("{:?}", initial_primitives);

    while t < T_FINAL {
        let primitives = prim_vec_from_cons(&conserved_vec, ADIABATIC, BX);
        let conserve = cons_vec_from_prim(&primitives, ADIABATIC);
        
        let mut dt = 1.0;
        for i in 0..(CELL_NUM-1) {
            let dt_check = math_func::compute_time_step(primitives[i], primitives[i+1], ADIABATIC, DR);
            if dt_check < dt {
                dt = dt_check; 
                }
            }
        dt = CFL * dt;
        
        conserved_vec = rk4_step(&primitives, conserve, dt);

        if t >= t_checkpoint {
            let _ = write_checkpoint(primitives, t, check_count);
            t_checkpoint += CHECK_INTERVAL;
            check_count += 1;
            }
        t += dt;
        time_step_count += 1.0;
    }
    let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
    let performance = (time_step_count * (CELL_NUM as f64)) / runtime;
    println!("The number of timesteps: {:?}", time_step_count);
    println!("The runtime in seconds: {:?}", runtime);
    println!("The performance in zones per second: {:.2?}", performance);
}