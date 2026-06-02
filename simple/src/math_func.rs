// This file is full of functions to supplement the simple 1D mhd code.
//
// Author: Brayden JoHantgen
// Last Update: 6/2/2026

/// Input:
///     prim: the eight component array of the primitive variables
/// Output:
///     p: the total pressure
/// Description:
///     This function takes the primitive variable and solves for the 
///     total pressure using this equation: P = P_gas + 0.5 * B^2.
pub fn total_pressure(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let p = prim.0 + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    p
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     E: the total energy
/// Description:
///     This function takes the primitive variable and solves for the total
///     energy using this equation: E = 0.5 * rho * v^2 + P/(rho * (gamma - 1)) + 0.5 * B^2.
pub fn total_energy(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    e
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     con: the seven component array of the conserved variables
/// Description:
///     This function constructs the conserved variable from the primitive variable.
///     The conserved variable has the components rho, rho * vx, rho * vy, rho * vz,
///     By, Bz, and E. This function uses the total energy function.
pub fn prim_to_cons(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
    let e = total_energy(prim, a_index);
    let con = (prim.1, prim.1 * prim.2, prim.1 * prim.3, prim.1 * prim.4, prim.6, prim.7, e);
    con
}

/// Input:
///     con: the seven component array of the conserved variables
///     a index: the adiabatic index
///     bx: the x component of the magnetic field
/// Output:
///     prim: the eight component array of the primitive variables
/// Description:
///     This function reconstructs the primirive variable from the conserved variable.
///     This function solves for the components of velocity and the pressure. The x
///     component of the magnetic field is held constant.
pub fn cons_to_prim(con: (f64, f64, f64, f64, f64, f64, f64), a_index: f64, bx: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let vx = con.1 / con.0;
    let vy = con.2 / con.0;
    let vz = con.3 / con.0;
    let p = (a_index - 1.0) * (con.6 - 0.5 * (con.0 * (vx * vx + vy * vy + vz * vz) + (con.4 * con.4 + con.5 * con.5 + bx * bx)));
    let prim = (p, con.0, vx, vy, vz, bx, con.4, con.5);
    prim
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     flux: the flux
/// Description:
///     This function takes the primitive variable to construct the flux.
///     This function uses the total pressure and total energy functions.
pub fn flux(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
    let p = total_pressure(prim);
    let e = total_energy(prim, a_index);
    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + p - prim.5 * prim.5;
    let f2 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f3 = prim.1 * prim.2 * prim.4 - prim.5 * prim.7;
    let f4 = prim.6 * prim.2 - prim.5 * prim.3;
    let f5 = prim.7 * prim.2 - prim.5 * prim.4;
    let f6 = (e + p) * prim.2 - prim.5 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_arr = (f0, f1, f2, f3, f4, f5, f6);
    f_arr
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     a: the sound speed of the gas
/// Description:
///     This function uses the primitive variables and the adiabatic index
///     to calculate the sound speed of the gas with the following equation:
///     cs = sqrt(gamma * P / rho). 
pub fn sound_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let a = (a_index * prim.0 / prim.1).sqrt();
    a
}

/// Input:
///     prim: the eight component array of the primitive variables
/// Output:
///     ca: the Alfven speed of the gas
/// Description:
///     This function takes the primitive variable to calculate the Alfven
///     speed of the gas using this equation: ca = |B_x| / sqrt(rho).
pub fn alfven_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let ca = prim.5.abs() / (prim.1).sqrt();
    ca
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     cf: the fast magnetosonic speed
/// Description:
///     This function takes the primitive variable to calculate the fast 
///     magnetosonic speed using this equation: 
///     cf = sqrt(gamma * P + B^2 + sqrt((gamma * P + B^2)^2 - 4 * gamma * P * B_x^2) / (2 * rho))
pub fn fast_magsonic_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let cf = ((a_index * prim.0 + b_squared + ((a_index * prim.0 + b_squared) * (a_index * prim.0 + b_squared) - 4.0 * a_index * prim.0 * prim.5 * prim.5).sqrt()) / (2.0 * prim.1)).sqrt();
    cf
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     cs: the slow magnetosonic speed
/// Description:
///     This function takes the primitive variable to calculate the slow 
///     magnetosonic speed using this equation: 
///     cf = sqrt(gamma * P + B^2 - sqrt((gamma * P + B^2)^2 - 4 * gamma * P * B_x^2) / (2 * rho))
pub fn slow_magsonic_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let cs = ((a_index * prim.0 + b_squared - ((a_index * prim.0 + b_squared) * (a_index * prim.0 + b_squared) - 4.0 * a_index * prim.0 * prim.5 * prim.5).sqrt()) / (2.0 * prim.1)).sqrt();
    cs
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     max: this is the max eigen value
/// Description:
///     This function uses the fast magneto sonic speed function. It uses
///     this function and the primitive variables and adiabatic index. It 
///     uses the equation: max = vx + cf.
pub fn max_eigen(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let cf = fast_magsonic_speed(prim, a_index);
    let max = prim.2 + cf;
    max
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     min: this is the min eigen value
/// Description:
///     This function uses the fast magneto sonic speed function. It uses
///     this function and the primitive variables and adiabatic index. It 
///     uses the equation: min = vx - cf.
pub fn min_eigen(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let cf = fast_magsonic_speed(prim, a_index);
    let min = prim.2 - cf;
    min
}

/// Input: 
///     tup: A tuple that contains three elements of type f64.
/// Output: 
///     The maximum value of the elements in the tuple, a float.
/// Description: 
///     This function converts the tuple to an array and then iterates over the elements of the array. 
///     If an array/tuple element is larger than zero it is saved and compared to the rest of the 
///     elements. The largest element is returned. This function will fail if all elements are negative.
pub fn tuple_max(tup: (f64, f64, f64)) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut max_check: f64 = 0.0;
    for i in 0..3 {
        if arr[i] > max_check {
            max_check = arr[i];
        }
    }
    max_check
}

/// Input:
///     tup: A tuple that contains three elements of type f64.
/// Output:
///     The minimum value of the elements in the tuple, a float.
/// Description:
///     This function uses the tuple max function to determine the max element of the tuple. It then
///     converts the tuple to an array and iterates over the array to find the minimum value.
pub fn tuple_min(tup: (f64, f64, f64)) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut min_check = tuple_max(tup);
    for i in 0..3 {
        if arr[i] < min_check {
            min_check = arr[i]
        }
    }
    min_check
}

/// Input:
///     prim l: the primitive variables of the left cell
///     prim r: the primitive variables of the right cell
///     a index: the adiabatic index
///     dx: the size of a cell
/// Output:
///     dt: the timestep
/// Description:
pub fn compute_time_step(prim_l: (f64, f64, f64, f64, f64, f64, f64, f64), prim_r: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64, dx: f64) -> f64 {
        let plus_l = max_eigen(prim_l, a_index);
        let minus_l = min_eigen(prim_l, a_index);
    
        let plus_r = max_eigen(prim_r, a_index);
        let minus_r = min_eigen(prim_r, a_index);
    
        let a_plus = tuple_max((0.0, plus_l, plus_r));
        let a_minus = tuple_max((0.0, -minus_l, -minus_r));
    
        let mut dt: f64 = 0.0;
    
        if a_minus > a_plus {
            dt += dx / a_minus;
        } else {
            dt += dx / a_plus;
        }    
        dt
    }

/// Input:
/// Output:
/// Description:
pub fn sgn(num: f64) -> f64 {
    let sign: f64;
    if num > 0.0 {
        sign = 1.0;
    } else if num < 0.0 {
        sign = -1.0;
    }
    else {
        sign = 0.0;
    }
    sign
}

/// Input:
/// Output:
/// Description:
pub fn minmod(x: f64, y: f64, z: f64) -> f64 {
    let mm_1 = (sgn(x) + sgn(y)).abs();
    let mm_2 = sgn(x) + sgn(z);
    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
    let mm = 0.25 * mm_1 * mm_2 * mm_3;
    mm
}

/// Input:
/// Output:
/// Description:
pub fn left_reconstruction(c_min: f64, c_mid: f64, c_max: f64) -> f64 {
    let arg1 = 1.5 * (c_mid - c_min);
    let arg2 = 0.5 * (c_max - c_mid);
    let arg3 = 1.5 * (c_max - c_mid);
    let cl = c_mid + 0.5 * minmod(arg1, arg2, arg3);
    cl
}

/// Input:
/// Output:
/// Description:
pub fn right_reconstruction(c_min: f64, c_mid: f64, c_max: f64) -> f64 {
    let arg1 = 1.5 * (c_mid - c_min);
    let arg2 = 0.5 * (c_max - c_min);
    let arg3 = 1.5 * (c_max - c_mid);
    let cr = c_mid - 0.5 * minmod(arg1, arg2, arg3);
    cr
}

// Input:
// Output:
// Description:
//pub fn vector_index(vec: Vec<u8>, desired_val: u8) -> Vec<u8> {
//    let mut count: u8 = 0;
//    let mut index_vec = Vec::new();
//    for i in vec {
//        if i == desired_val{
//            index_vec.push(count);
//        }
//        count += 1;
//    }
//    index_vec
//}