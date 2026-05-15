// This is the rust file containing the 1D relativistic mhd simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/15/2026

//#![allow(dead_code)]

// use std::fs;
//use std::io::{BufWriter, Write};
//use std::path::Path;
//use std::env;

pub mod math_functions;

/////////////////////
// Useful Variables
/////////////////////
const CELL_NUM: f64 = 10.0;
const DISCON: f64 = 0.5;
//const DR: f64 = 1.0 / CELL_NUM;
//const ADIABATIC: f64 = 1.4;
//const T_FINAL: f64 = 0.401;
//const CHECK_INTERVAL: f64 = 0.025;
//const CFL: f64 = 0.5;

////////////////
// Dataclasses
////////////////

///////////////////
// Math Functions
///////////////////

//fn p_eigen(prim: (f64, f64, f64), a_index: f64) -> f64 {
//    let cs = sound_speed(prim.clone(), a_index);
//    let max_eigen = prim.2 + cs;
//    max_eigen
//}

//fn m_eigen(prim: (f64, f64, f64), a_index: f64) -> f64 {
//    let cs = sound_speed(prim.clone(), a_index);
//   let min_eigen = prim.2 - cs;
//    min_eigen
//}

//fn sgn(num: f64) -> f64 {
//    let mut sign: f64;
//    if num > 0.0 {
//        sign = 1.0;
//    } else if num < 0.0 {
//        sign = -1.0;
//    }
//    else {
//        sign = 0.0;
//    }
//    sign
//}

//fn minmod(x: f64, y: f64, z: f64) -> f64 {
//    let mm_1 = (sgn(x) + sgn(y)).abs();
//    let mm_2 = sgn(x) + sgn(z);
//    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
//    let mm = 0.25 * mm_1 * mm_2 * mm_3;
//    mm
//}

//fn compute_time_step(prim_l: (f64, f64, f64), prim_r: (f64, f64, f64), a_index: f64) -> f64 {
//    let plus_l = p_eigen(prim_l.clone(), a_index);
//    let minus_l = m_eigen(prim_l.clone(), a_index);

//    let plus_r = p_eigen(prim_r.clone(), a_index);
//    let minus_r = m_eigen(prim_r.clone(), a_index);

//    let a_plus = tuple_max((0.0, plus_l, plus_r));
//    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

//    let mut dt: f64 = 0.0;

//    if a_minus > a_plus {
//        dt += DR / a_minus;
//    } else {
//        dt += DR / a_plus;
//    }    
//    dt
//}

////////////////////
// Usage Functions
////////////////////

/// Input:
///     There is no input for this function.
/// Output:
///     init_primitive: this is a vector with six f64 components which
///         correspont to P, rho, v_x, v_y, B_x, B_y.
/// Description:
///     This function uses the values for the cell numbers and the location
///     of the discontinuity to write the primitive variables for the defined
///     grid. This function depends on the read_config function.
fn init_prim() -> Vec<(f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive = Vec::new();
    for i in 0..(CELL_NUM as u8) {
        if i < ((CELL_NUM * DISCON) as u8) {
            init_primitive.push((1.0, 1.0, 0.0, 0.0, 0.5, 1.0));
        } else {
            init_primitive.push((0.125, 0.1, 0.0, 0.0, 0.5, -1.0));
        }
    }
    init_primitive
}

//fn cons_vec_from_prim(prims: Vec<(f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64)> {
//    let mut cons_vec = Vec::new();
//    for i in 0..(CELL_NUM as u8) {
//        let index: usize = (i).try_into().unwrap();
//        cons_vec.push(prim_to_cons(prims[index], a_index));
//    }
//    cons_vec
//}

//fn prim_vec_from_cons(cons: Vec<(f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64)> {
//    let mut prims_vec = Vec::new();
//    for i in 0..(CELL_NUM as u8) {
//        let index: usize = (i).try_into().unwrap();
//        prims_vec.push(cons_to_prim(cons[index], a_index));
//    }
//    prims_vec
//}

//fn read_config() {
//    let file_path = "/Users/toogan13/Desktop/MHD_Wind/magneto_config.txt".to_string();
//    let file = fs::read_to_string(&file_path);
//    let breaking = file.expect("REASON");
//    let file_rows: Vec<&str> = breaking.split('\n').collect();
//
//    println!("{:?}", file_rows);
//}

//fn write_checkpoint(prims: Vec<(f64, f64, f64)>, t: f64, check_count: i8) -> Result<(), Box<dyn std::error::Error>> {
//    let file_num = check_count.to_string();
//    let file_type = ".txt".to_string();
//    let file_name = format!("{}{}", file_num, file_type);
//    let file_path = "time_step_files/".to_string() + &file_name;
//    let output_file_path = Path::new(&file_path);
//
//    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));
//
//    let num_fill = format!("{} {}", "cell_num:".to_string(), &(CELL_NUM.to_string()+&" ".to_string()));
//
//    let mut p_string = "p: ".to_string();
//    let mut rho_string = "rho: ".to_string();
//    let mut v_string = "v: ".to_string();
//    for i in prims {
//        let p_fill = i.0.to_string() + &" ".to_string();
//        p_string.push_str(&p_fill);
//
//        let rho_fill = i.1.to_string() + &" ".to_string();
//        rho_string.push_str(&rho_fill);
//
//        let v_fill = i.2.to_string() + &" ".to_string();
//        v_string.push_str(&v_fill);
//    }
// 
//    let file = fs::OpenOptions::new()
//        .append(true)
//        .create(true)
//        .open(output_file_path)?;
//
//    let mut writer = BufWriter::new(file);
//
//    writeln!(writer, "{}", num_fill)?;
//    writeln!(writer, "{}", t_fill)?;
//    writeln!(writer, "{}", p_string)?;
//    writeln!(writer, "{}", rho_string)?;
//    writeln!(writer, "{}", v_string)?;
//    writer.flush()?;
//
//    Ok(())
//}

///////////////
// Simulation
///////////////
fn main() {
//    let mut t: f64 = 0.0;
//    let mut t_checkpoint = CHECK_INTERVAL;
//    let mut check_count: i8 = 0;
//
    let initial_primitives = init_prim();
    let test_val = math_functions::total_energy_density(initial_primitives[0], 1.4);
    println!("{:?}", test_val);
//    let mut conserved_vec = cons_vec_from_prim(initial_primitives.clone(), ADIABATIC);
//
//    while t < T_FINAL {
//        let primitives = prim_vec_from_cons(conserved_vec.clone(), ADIABATIC);
//        let conserve = cons_vec_from_prim(primitives.clone(), ADIABATIC);
//
//        let mut dt = 1.0;
//        for i in 0..((CELL_NUM - 1.0) as u8) {
//            let index_1: usize = (i).try_into().unwrap();
//            let index_2: usize = (i+1).try_into().unwrap();
//            let dt_check = compute_time_step(primitives[index_1], primitives[index_2], ADIABATIC);
//            if dt_check < dt {
//                dt = dt_check;
//            }
//        }
//        dt = CFL * dt;
//
//        conserved_vec = l_function(primitives.clone(), conserve, dt);
//
//        println!("First Cell: {:?}", primitives[0]);
//        println!("Middle Cell: {:?}", primitives[5]);
//        println!("Last Cell: {:?}", primitives[9]);
//        if t >= t_checkpoint {
//            write_checkpoint(primitives.clone(), t, check_count);
//            t_checkpoint += CHECK_INTERVAL;
//            check_count += 1;
//        }
//
//        t += dt;
//    }
} 