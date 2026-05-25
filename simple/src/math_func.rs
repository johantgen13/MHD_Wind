// This file is full of functions to supplement the simple mhd code.
//
// Author: Brayden JoHantgen
// Last Update: 5/25/2026

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

/// Input:
/// Output:
/// Description:
pub fn sound_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let a = (a_index * prim.0 / prim.1).sqrt();
    a
}

/// Input:
/// Output:
/// Description:
pub fn alfven_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let ca = prim.5.abs() / (prim.1).sqrt();
    ca
}

/// Input:
/// Output:
/// Description:
pub fn fast_magsonic_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let cf = ((a_index * prim.0 + b_squared + ((a_index * prim.0 + b_squared) * (a_index * prim.0 + b_squared) - 4.0 * a_index * prim.0 * prim.5 * prim.5).sqrt()) / (2.0 * prim.1)).sqrt();
    cf
}

/// Input:
/// Output:
/// Description:
pub fn slow_magsonic_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let cs = ((a_index * prim.0 + b_squared - ((a_index * prim.0 + b_squared) * (a_index * prim.0 + b_squared) - 4.0 * a_index * prim.0 * prim.5 * prim.5).sqrt()) / (2.0 * prim.1)).sqrt();
    cs
}

/// Input:
/// Output:
/// Description:
pub fn max_eigen(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let cf = fast_magsonic_speed(prim.clone(), a_index);
    let max = prim.2 + cf;
    max
}

/// Input:
/// Output:
/// Description:
pub fn min_eigen(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let cf = fast_magsonic_speed(prim.clone(), a_index);
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
/// Output:
/// Description:
pub fn compute_time_step(prim_l: (f64, f64, f64, f64, f64, f64, f64, f64), prim_r: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64, dx: f64) -> f64 {
        let plus_l = max_eigen(prim_l.clone(), a_index);
        let minus_l = min_eigen(prim_l.clone(), a_index);
    
        let plus_r = max_eigen(prim_r.clone(), a_index);
        let minus_r = min_eigen(prim_r.clone(), a_index);
    
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