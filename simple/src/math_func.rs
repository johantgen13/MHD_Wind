// This file is full of functions to supplement the simple mhd code.
//
// Author: Brayden JoHantgen
// Last Update: 5/22/2026

/// Input:
/// Output:
/// Description:
pub fn spec_int_energy_gas(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = prim.0 / (prim.1 * (a_index - 1.0));
    e
}

/// Input:
/// Output:
/// Description:
pub fn b_zero(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b0 = prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7;
    b0
}

/// Input:
/// Output:
/// Description:
pub fn b_x(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b0 = b_zero(prim.clone());
    let bx = prim.5 + b0 * prim.2;
    bx
}

/// Input:
/// Output:
/// Description:
pub fn b_y(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b0 = b_zero(prim.clone());
    let by = prim.6 + b0 * prim.3;
    by
}

/// Input:
/// Output:
/// Description:
pub fn b_z(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b0 = b_zero(prim.clone());
    let bz = prim.7 + b0 * prim.4;
    bz
}

/// Input:
/// Output:
/// Description:
pub fn b_squared(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b0 = b_zero(prim.clone());
    let bsq = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7 + b0 * b0;
    bsq
}

/// Input:
/// Output:
/// Description:
pub fn spec_int_energy_total(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = spec_int_energy_gas(prim.clone(), a_index);
    let bsq = b_squared(prim.clone());
    let e_star = e + bsq / (2.0 * prim.1);
    e_star
}

/// Input:
/// Output:
/// Description:
pub fn pressure_total(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let bsq = b_squared(prim.clone());
    let p_star = prim.0 + bsq / 2.0;
    p_star
}

/// Input:
/// Output:
/// Description:
pub fn spec_enthalpy_total(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e_star = spec_int_energy_total(prim.clone(), a_index);
    let p_star = pressure_total(prim.clone());
    let h_star = 1.0 + e_star + p_star / prim.1;
    h_star
}

/// Input:
/// Output:
/// Description:
pub fn mom_den_x(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b0 = b_zero(prim.clone());
    let bx = b_x(prim.clone());
    let h_star = spec_enthalpy_total(prim.clone(), a_index);
    let sx = prim.1 * h_star * prim.2 - b0 * bx;
    sx
}

/// Input:
/// Output:
/// Description:
pub fn mom_den_y(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b0 = b_zero(prim.clone());
    let by = b_y(prim.clone());
    let h_star = spec_enthalpy_total(prim.clone(), a_index);
    let sy = prim.1 * h_star * prim.3 - b0 * by;
    sy
}

/// Input:
/// Output:
/// Description:
pub fn mom_den_z(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b0 = b_zero(prim.clone());
    let bz = b_z(prim.clone());
    let h_star = spec_enthalpy_total(prim.clone(), a_index);
    let sz = prim.1 * h_star * prim.4 - b0 * bz;
    sz
}

/// Input:
/// Output:
/// Description:
pub fn energy_den_total(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b0 = b_zero(prim.clone());
    let h_star = spec_enthalpy_total(prim.clone(), a_index);
    let p_star = pressure_total(prim.clone());
    let tau = prim.1 * h_star - p_star - b0 * b0 - prim.1;
    tau
}

/// Input:
/// Output:
/// Description:
pub fn prim_to_cons(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let sx = mom_den_x(prim.clone(), a_index);
    let sy = mom_den_y(prim.clone(), a_index);
    let sz = mom_den_z(prim.clone(), a_index);
    let tau = energy_den_total(prim.clone(), a_index);
    let con = (prim.1, sx, sy, sz, tau, prim.5, prim.6, prim.7);
    con
}

/// Input:
/// Output:
/// Description:
pub fn cubic_solver(a: f64, b: f64, c:f64, d:f64) -> f64 {
    let delta_0 = b * b - 3.0 * a * c;
    let delta_1 = 2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d;
    let delta_2 = (0.5 * (delta_1 - (delta_1 * delta_1 - 4.0 * delta_0 * delta_0 * delta_0).sqrt())).cbrt();
    let x = - 0.3333 * (b + delta_2 + delta_0 / delta_2) / a;
    x
}

/// Input:
/// Output:
/// Description:
pub fn cons_to_prim(con: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b2 = con.5 * con.5 + con.6 * con.6 + con.7 * con.7;
    let a = 2.0 / (a_index - 1.0);
    let b = 4.0 * con.0 / (a_index - 1.0) + b2 - 2.0 * con.4 - 1.0 / ((a_index - 1.0) * (a_index - 1.0));
    let c = 2.0 * con.0 * con.0 / (a_index - 1.0) + 2.0 * con.0 * b2 - 4.0 * con.0 * con.4 + 2.0 * (con.1 * con.5 + con.2 * con.6 + con.3 * con.7) / (a_index - 1.0);
    let d = b2 * con.0 * con.0 - 2.0 * con.0 * con.0 * con.4 - (con.1 * con.5 + con.2 * con.6 + con.3 * con.7) * (con.1 * con.5 + con.2 * con.6 + con.3 * con.7);
    let p = cubic_solver(a, b, c, d);
    p
}