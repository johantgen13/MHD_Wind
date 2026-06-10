// This is a simple non-relativisetic 1D MHD code.
//
// Author: Brayden JoHantgen
// Last Update: 6/9/2026

use std::fs;
//use std::str;
use std::io::{BufWriter, Write};//, BufRead, BufReader, Read};
use std::path::Path;
use std::time::Instant;
use rand::Rng;
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
fn init_prims_1d(phys: &Physics, zones: usize, discon: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut init_primitive: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
    for i in 0..zones {
        if i < (((zones as f64) * discon) as usize) {
            init_primitive.push((phys.p.0, phys.rho.0, phys.vx.0, phys.vy.0, phys.vz.0, phys.bx.0, phys.by.0, phys.bz.0));
        } else {
            init_primitive.push((phys.p.1, phys.rho.1, phys.vx.1, phys.vy.1, phys.vz.1, phys.bx.1, phys.by.1, phys.bz.1));
        }
    }
    init_primitive
}

/// Input:
/// Output:
/// Description:
fn init_prims_2d(phys: &Physics, drive: &Driver) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut init_primitive: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    for i in 0..drive.num_zones_x {
        let mut prim_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        for j in 0..drive.num_zones_y {
            let mut rng = rand::thread_rng();
            let num1: f64 = rng.gen_range(0.0..0.01);
            let num2: f64 = rng.gen_range(0.0..0.01);
            if j < (((drive.num_zones_y as f64) * drive.discontinuity) as usize) {
                prim_fill.push((phys.p.1, phys.rho.1, phys.vx.1+num1, phys.vy.1+num2, phys.vz.1, phys.bx.1, phys.by.1, phys.bz.1));
            } else if j > (((drive.num_zones_y as f64) * (1.0 - drive.discontinuity)) as usize) {
                prim_fill.push((phys.p.1, phys.rho.1, phys.vx.1+num1, phys.vy.1+num2, phys.vz.1, phys.bx.1, phys.by.1, phys.bz.1));
            } else {
                prim_fill.push((phys.p.0, phys.rho.0, phys.vx.0+num1, phys.vy.0+num2, phys.vz.0, phys.bx.0, phys.by.0, phys.bz.0));
            }
        }
        init_primitive.push(prim_fill);
    }
    init_primitive
}

/// Input:
/// Output:
/// Description:
fn cons_vec_from_prim_1d(prims: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut cons_vec = Vec::new();
    for i in 0..(prims.len() as usize) {
        cons_vec.push(math_func::prim_to_cons(prims[i], a_index));
    }
    cons_vec
}

/// Input:
/// Output:
/// Description:
fn cons_vec_from_prim_2d(prims: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut cons_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    for i in 0..x_zone {
        let mut cons_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new(); 
        for j in 0..y_zone {
            cons_fill.push(math_func::prim_to_cons(prims[i][j], a_index));
        }
        cons_vec.push(cons_fill);
    }
    cons_vec
}

/// Input:
/// Output:
/// Description:
fn prim_vec_from_cons_1d(cons: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>, a_index: f64) -> Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> {
    let mut prims_vec = Vec::new();
    for i in 0..(cons.len() as usize) {
        prims_vec.push(math_func::cons_to_prim(cons[i], a_index));
    }
    prims_vec
}

/// Input:
/// Output:
/// Description:
fn prim_vec_from_cons_2d(cons: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut prims_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    for i in 0..x_zone {
        let mut prim_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        for j in 0..y_zone {
            prim_fill.push(math_func::cons_to_prim(cons[i][j], a_index));
        }
        prims_vec.push(prim_fill);
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
fn godonov_x_2d(prims: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut go_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    let mut go: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
    let mut count: usize = 0;
    while count < y_zone {
        go.push(hll_flux_x(prims[x_zone-1][count], prims[0][count], a_index));
        count += 1;
    }
    go_vec.push(go.clone());
    for i in 1..x_zone {
        let mut go_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        for j in 0..y_zone {
            go_fill.push(hll_flux_x(prims[i-1][j], prims[i][j], a_index));
        }
        go_vec.push(go_fill);
    }
    go_vec.push(go);
    go_vec
}

/// Input:
/// Output:
/// Description:
fn godonov_y_2d(prims: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut go_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    for i in 0..x_zone {
        let mut go_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        go_fill.push(hll_flux_y(prims[i][y_zone-1], prims[i][0], a_index));
        for j in 1..y_zone {
            go_fill.push(hll_flux_y(prims[i][j-1], prims[i][j], a_index));
        }
        go_fill.push(hll_flux_y(prims[i][y_zone-1], prims[i][0], a_index));
        go_vec.push(go_fill);
    }
    go_vec
}

/// Input:
/// Output:
/// Description:
fn l_function_2d(prims: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let mut new_cons_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    let f_vec = godonov_x_2d(prims.clone(), x_zone, y_zone, a_index);
    let g_vec = godonov_y_2d(prims.clone(), x_zone, y_zone, a_index);
    let dx: f64 = 1.0 / (x_zone as f64);
    let dy: f64 = 1.0 / (y_zone as f64);

    for i in 0..x_zone {
        let mut new_fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        for j in 0..y_zone {
            let new_0 = - (f_vec[i+1][j].0 - f_vec[i][j].0) / dx - (g_vec[i][j+1].0 - g_vec[i][j].0) / dy;
            let new_1 = - (f_vec[i+1][j].1 - f_vec[i][j].1) / dx - (g_vec[i][j+1].1 - g_vec[i][j].1) / dy;
            let new_2 = - (f_vec[i+1][j].2 - f_vec[i][j].2) / dx - (g_vec[i][j+1].2 - g_vec[i][j].2) / dy;
            let new_3 = - (f_vec[i+1][j].3 - f_vec[i][j].3) / dx - (g_vec[i][j+1].3 - g_vec[i][j].3) / dy;
            let new_4 = - (f_vec[i+1][j].4 - f_vec[i][j].4) / dx - (g_vec[i][j+1].4 - g_vec[i][j].4) / dy;
            let new_5 = - (f_vec[i+1][j].5 - f_vec[i][j].5) / dx - (g_vec[i][j+1].5 - g_vec[i][j].5) / dy;
            let new_6 = - (f_vec[i+1][j].6 - f_vec[i][j].6) / dx - (g_vec[i][j+1].6 - g_vec[i][j].6) / dy;
            let new_7 = - (f_vec[i+1][j].7 - f_vec[i][j].7) / dx - (g_vec[i][j+1].7 - g_vec[i][j].7) / dy;
            new_fill.push((new_0, new_1, new_2, new_3, new_4, new_5, new_6, new_7));
        }
        new_cons_vec.push(new_fill);
    }
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
fn euler_timestep_2d(prims_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, cons_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, a_index: f64, dt: f64) -> Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> {
    let l_cons = l_function_2d(prims_vec.clone(), x_zone, y_zone, a_index);
    let mut new_cons_vec: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>> = Vec::new();
    for i in 0..x_zone {
        let mut fill: Vec<(f64, f64, f64, f64, f64, f64, f64, f64)> = Vec::new();
        for j in 0..y_zone {
            let fill_0 = cons_vec[i][j].0 + dt * l_cons[i][j].0;
            let fill_1 = cons_vec[i][j].1 + dt * l_cons[i][j].1;
            let fill_2 = cons_vec[i][j].2 + dt * l_cons[i][j].2;
            let fill_3 = cons_vec[i][j].3 + dt * l_cons[i][j].3;
            let fill_4 = cons_vec[i][j].4 + dt * l_cons[i][j].4;
            let fill_5 = cons_vec[i][j].5 + dt * l_cons[i][j].5;
            let fill_6 = cons_vec[i][j].6 + dt * l_cons[i][j].6;
            let fill_7 = cons_vec[i][j].7 + dt * l_cons[i][j].7;
            fill.push((fill_0, fill_1, fill_2, fill_3, fill_4, fill_5, fill_6, fill_7));
        }
        new_cons_vec.push(fill);
    }
    new_cons_vec
}

/// Input:
/// Output:
/// Description:
fn write_checkpoint_2d(prims: Vec<Vec<(f64, f64, f64, f64, f64, f64, f64, f64)>>, x_zone: usize, y_zone: usize, t: f64, check_count: i8) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);

    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));

    let num_fill_x = format!("{} {}", "cell_num_x:".to_string(), &(x_zone.to_string()+&" ".to_string()));
    let num_fill_y = format!("{} {}", "cell_num_y:".to_string(), &(y_zone.to_string()+&" ".to_string()));

    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut vx_string = "vx: ".to_string();
    let mut vy_string = "vy: ".to_string();
    let mut vz_string = "vz: ".to_string();
    let mut bx_string = "Bx: ".to_string();
    let mut by_string = "By: ".to_string();
    let mut bz_string = "Bz: ".to_string();

    let mut prims_fill = Vec::new();
    for i in 0..x_zone {
        for j in 0..y_zone {
            prims_fill.push(prims[i][j]);
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

    writeln!(writer, "{}", num_fill_x)?;
    writeln!(writer, "{}", num_fill_y)?;
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
    let phys = Physics{adiabatic_index: 1.4, p: (2.5, 2.5), rho: (2.0, 1.0), vx: (0.5, -0.5), vy: (0.0, 0.0), vz: (0.0, 0.0), bx: (1.772, 1.772), by: (0.0, 0.0), bz: (0.0, 0.0)};
    let drive = Driver{cfl: 0.4, tfinal: 5.001, checkpoint: 0.0125, num_zones_x: 512, num_zones_y: 512, discontinuity: 0.25, dimensionality: "2D".to_string(), plm: false, grid_type: "Cartesian".to_string()};

    let before = Instant::now();

    //let path = "/Users/toogan13/Desktop/MHD_Wind/magneto_config.txt";
    //let config = read_config(path.to_string());
    //println!("{:?}", config);

    let mut t: f64 = 0.0;
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i8 = 0;

    let initial_primitives = init_prims_2d(&phys, &drive);
    let mut conserved_vec = cons_vec_from_prim_2d(initial_primitives.clone(), drive.num_zones_x, drive.num_zones_y, phys.adiabatic_index);

    while t < drive.tfinal {
        let primitives = prim_vec_from_cons_2d(conserved_vec.clone(), drive.num_zones_x, drive.num_zones_y, phys.adiabatic_index);
        let conserve = cons_vec_from_prim_2d(primitives.clone(), drive.num_zones_x, drive.num_zones_y, phys.adiabatic_index);
        
        let mut dt = 1.0;
        for i in 0..(drive.num_zones_x-1) {
            for j in 0..(drive.num_zones_y-1) {
                let dt_check = math_func::compute_time_step(primitives[i][j], primitives[i+1][j], phys.adiabatic_index, 1.0 / (drive.num_zones_x as f64));
                if dt_check < dt {
                    dt = dt_check; 
                    }
                }
            }
        dt = drive.cfl * dt;
        
        conserved_vec = euler_timestep_2d(primitives.clone(), conserve, drive.num_zones_x, drive.num_zones_y, phys.adiabatic_index, dt);

        if t >= t_checkpoint {
            let _ = write_checkpoint_2d(primitives.clone(), drive.num_zones_x, drive.num_zones_y, t, check_count);
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