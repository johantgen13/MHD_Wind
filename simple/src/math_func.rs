// This file is full of functions to supplement the simple 1D mhd code.
//
// Author: Brayden JoHantgen
// Last Update: 6/7/2026

/// Input:
///     prim: the eight component array of the primitive variables
/// Output:
///     p: the total pressure
/// Description:
///     This function takes the primitive variable and solves for the 
///     total pressure using this equation: P = P_gas + 0.5 * B^2.
#[inline(always)]
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
///     energy using this equation: E = 0.5 * rho * v^2 + P/(gamma - 1) + 0.5 * B^2.
#[inline(always)]
pub fn total_energy(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    e
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     con: the eight component array of the conserved variables
/// Description:
///     This function constructs the conserved variable from the primitive variable.
///     The conserved variable has the components rho, rho * vx, rho * vy, rho * vz,
///     By, Bz, and E. This function uses the total energy function.
pub fn prim_to_cons(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let e = total_energy(prim.clone(), a_index);
    let con = (prim.1, prim.1 * prim.2, prim.1 * prim.3, prim.1 * prim.4, prim.5, prim.6, prim.7, e);
    con
}

/// Input:
///     con: the eight component array of the conserved variables
///     a index: the adiabatic index
///     bx: the x component of the magnetic field
/// Output:
///     prim: the eight component array of the primitive variables
/// Description:
///     This function reconstructs the primirive variable from the conserved variable.
///     This function solves for the components of velocity and the pressure. The x
///     component of the magnetic field is held constant.
pub fn cons_to_prim(con: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let vx = con.1 / con.0;
    let vy = con.2 / con.0;
    let vz = con.3 / con.0;
    let p = (a_index - 1.0) * (con.7 - 0.5 * (con.0 * (vx * vx + vy * vy + vz * vz) + (con.4 * con.4 + con.5 * con.5 + con.6 * con.6)));
    let prim = (p, con.0, vx, vy, vz, con.4, con.5, con.6);
    prim
}

/// Input:
///     prim: the eight component array of the primitive variables
///     a index: the adiabatic index
/// Output:
///     flux: the flux in the x-direction
/// Description:
///     This function takes the primitive variable to construct the flux in the x-direction.
///     This function uses the total pressure and total energy functions.
pub fn flux_x(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let p = total_pressure(prim);
    let e = total_energy(prim, a_index);
    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + p - prim.5 * prim.5;
    let f2 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f3 = prim.1 * prim.2 * prim.4 - prim.5 * prim.7;
    let f4 = 0.0;
    let f5 = prim.6 * prim.2 - prim.5 * prim.3;
    let f6 = prim.7 * prim.2 - prim.5 * prim.4;
    let f7 = (e + p) * prim.2 - prim.5 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_x = (f0, f1, f2, f3, f4, f5, f6, f7);
    f_x
}

/// Input:
/// Output:
/// Description:
pub fn flux_y(prim:(f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let p = total_pressure(prim);
    let e = total_energy(prim, a_index);
    let f0 = prim.1 * prim.3;
    let f1 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f2 = prim.1 * prim.3 * prim.3 + p - prim.6 * prim.6;
    let f3 = prim.1 * prim.3 * prim.4 - prim.6 * prim.7;
    let f4 = prim.5 * prim.3 - prim.6 * prim.2;
    let f5 = 0.0;
    let f6 = prim.7 * prim.3 - prim.6 * prim.4;
    let f7 = (e + p) * prim.3 - prim.6 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_y = (f0, f1, f2, f3, f4, f5, f6, f7);
    f_y
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
#[inline(always)] 
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
#[inline(always)]
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
///     tup: A tuple that contains three elements of type f64.
/// Output: 
///     The maximum value of the elements in the tuple, a float.
/// Description: 
///     This function converts the tuple to an array and then iterates over the elements of the array. 
///     If an array/tuple element is larger than zero it is saved and compared to the rest of the 
///     elements. The largest element is returned. This function will fail if all elements are negative.
#[inline(always)]
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
#[inline(always)]
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
    let cf_l = fast_magsonic_speed(prim_l, a_index);    
    let plus_x_l = cf_l + prim_l.2;
    let plus_y_l = cf_l + prim_l.3; 
    let minus_x_l = prim_l.2 - cf_l;
    let minus_y_l = prim_l.3 - cf_l;

    let cf_r = fast_magsonic_speed(prim_r, a_index);    
    let plus_x_r = cf_r + prim_r.2;
    let plus_y_r = cf_r + prim_r.3; 
    let minus_x_r = prim_r.2 - cf_r;
    let minus_y_r = prim_r.3 - cf_r;
    
    let a_plus_l = tuple_max((0.0, plus_x_l, plus_y_l));
    let a_minus_l = tuple_max((0.0, -minus_x_l, -minus_y_l));
    
    let a_plus_r = tuple_max((0.0, plus_x_r, plus_y_r));
    let a_minus_r = tuple_max((0.0, -minus_x_r, -minus_y_r));

    let mut a_minus: f64 = 0.0;
    let mut a_plus: f64 = 0.0;

    if a_minus_l > a_minus_r {
        a_minus = a_minus_l;
    } else {
        a_minus = a_minus_r;
    }

    if a_plus_l > a_plus_r {
        a_plus = a_plus_l;
    } else {
        a_plus = a_plus_r;
    }

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
#[inline(always)]
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
#[inline(always)]
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
#[inline(always)]
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
#[inline(always)]
pub fn right_reconstruction(c_min: f64, c_mid: f64, c_max: f64) -> f64 {
    let arg1 = 1.5 * (c_mid - c_min);
    let arg2 = 0.5 * (c_max - c_min);
    let arg3 = 1.5 * (c_max - c_mid);
    let cr = c_mid - 0.5 * minmod(arg1, arg2, arg3);
    cr
}