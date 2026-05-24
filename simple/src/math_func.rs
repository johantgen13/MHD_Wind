// This file is full of functions to supplement the simple mhd code.
//
// Author: Brayden JoHantgen
// Last Update: 5/24/2026

/// Input:
/// Output:
/// Description:
pub fn total_pressure(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let p = prim.0 + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    p
}

/// Input:
/// Output:
/// Description:
pub fn total_energy(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    e
}

/// Input:
/// Output:
/// Description:
pub fn prim_to_cons(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let e = total_energy(prim.clone(), a_index);
    let con = (prim.1, prim.1 * prim.2, prim.1 * prim.3, prim.1 * prim.4, prim.5, prim.6, prim.7, e);
    con
}

/// Input:
/// Output:
/// Description:
pub fn cons_to_prim(con: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let vx = con.1 / con.0;
    let vy = con.2 / con.0;
    let vz = con.3 / con.0;
    let p = (a_index - 1.0) * (con.7 - 0.5 * (con.0 * (vx * vx + vy * vy + vz * vz) + (con.4 * con.4 + con.5 * con.5 + con.6 * con.6)));
    let prim = (p, con.0, vx, vy, vz, con.4, con.5, con.6);
    prim
}

/// Input:
/// Output:
/// Description:
pub fn flux(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let p = total_pressure(prim.clone());
    let e = total_energy(prim.clone(), a_index);
    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + p;
    let f2 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f3 = prim.1 * prim.2 * prim.4 - prim.5 * prim.7;
    let f5 = prim.6 * prim.2 - prim.5 * prim.3;
    let f6 = prim.7 * prim.2 - prim.5 * prim.4;
    let f7 = (e + p) * prim.2 - prim.5 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_arr = (f0, f1, f2, f3, 0.0, f5, f6, f7);
    f_arr
}