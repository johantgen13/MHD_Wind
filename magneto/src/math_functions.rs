// This is the rust file containing the math functions for the simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/18/2026

/// Input:
///     vx: The x-component of the fluid velocity.
///     vy: The y-component of the fluid velocity.
/// Output:
///     w: The Lorentz factor
/// Description:
///     This function produces the Lorentz factor of the fluid.
///     It uses the equation: (1 - v^2)^(-1/2).
pub fn lorentz_factor(vx: f64, vy: f64, vz: f64) -> f64 {
    let w = 1.0 / (1.0 - (vx * vx + vy * vy + vz * vz)).sqrt();
    w
}

/// Input:
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     d: The rest_mass_density
/// Description:
///     This function uses the lorentz factor function to solve for the rest
///     mass density by using the following equation: D = rho * w.
pub fn rest_mass_density(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let d = prim.1 * w;
    d
}

/// Input: 
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     b_0: The zeroth-component of the magnetic four-vector.
/// Description:
///     This function uses the lorentz factor function. The zeroth-component of the
///     magnetic four-vector is given by: b_0 = w * B^i * v_i.
pub fn b_zero(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = w * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    b_0
}

/// Input:
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     b_x: The x-component of the magnetic four-vector.
/// Description:
///     This function uses the lorentz factor function. The x-component of the 
///     magnetic four-vector is given by : B_x / w + b_0 * v_x.
pub fn b_x(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let x_component = (prim.5 / w) + b_0 * prim.2;
    x_component
}

/// Input:
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     b_y: The y-component of the magnetic four-vector.
/// Description:
///     This function uses the lorentz factor function. The y-component of the 
///     magnetic four-vector is given by : B_y / w + b_0 * v_y.
pub fn b_y(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let y_component = (prim.6 / w) + b_0 * prim.3;
    y_component
}

/// Input:
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     b_z: The z-component of the magnetic four-vector.
/// Description:
///     This function uses the lorentz factor function. The z-component of the 
///     magnetic four-vector is given by : B_z / w + b_0 * v_z.
pub fn b_z(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let z_component = (prim.7 / w) + b_0 * prim.4;
    z_component
}

/// Input:
///     prim: The eight component primitive variable tuple for an individual cell.
/// Output:
///     b^2: The square of the magnetic four-vector
/// Description:
///     This function uses the lorentz factor function, the b_zero function, the 
///     b_x function, the b_y function, and the b_z function. This function uses
///     the equation: (b^mu * b_mu) / (4 * pi).
pub fn b_square(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let b_0 = b_zero(prim.clone());
    let b_1 = b_x(prim.clone());
    let b_2 = b_y(prim.clone());
    let b_3 = b_z(prim.clone());
    let bsq = (b_0 * b_0 + b_1 * b_1 + b_2 * b_2 + b_3 * b_3) / (4.0 * 3.14159);
    bsq
}

/// Input:
/// Output:
/// Description:
pub fn specific_energy_gas(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 { 
    let e = prim.0 / ((a_index - 1.0) * prim.1);
    e
}

/// Input:
/// Output:
/// Description:
pub fn sound_speed(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let cs = ((a_index * prim.0) / prim.1).sqrt();
    cs
}

/// Input:
/// Output:
/// Description:
pub fn specific_energy_total(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e = specific_energy_gas(prim.clone(), a_index);
    let bsq = b_square(prim.clone());
    let e_star = e + (bsq / (2.0 * prim.1));
    e_star
}

/// Input:
/// Output:
/// Description:
pub fn total_pressure(prim: (f64, f64, f64, f64, f64, f64, f64, f64)) -> f64 {
    let bsq = b_square(prim.clone());
    let p_star = prim.0 + (bsq / 2.0);
    p_star
}

/// Input:
/// Output:
/// Description:
pub fn specific_enthalpy(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let e_star = specific_energy_total(prim.clone(), a_index);
    let p_star = total_pressure(prim.clone());
    let h_star = 1.0 + e_star + (p_star / prim.1);
    h_star
}

/// Input:
/// Output:
/// Description:
pub fn momentum_density_x(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let h_star = specific_enthalpy(prim.clone(), a_index);
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let b_xcomp = b_x(prim.clone());
    let s_x = (prim.1 * h_star * w * w * prim.2) - ((b_0 * b_xcomp) / (4.0 * 3.14159));
    s_x
}

/// Input:
/// Output:
/// Description:
pub fn momentum_density_y(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let h_star = specific_enthalpy(prim.clone(), a_index);
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let b_ycomp = b_y(prim.clone());
    let s_y = (prim.1 * h_star * w * w * prim.3) - ((b_0 * b_ycomp) / (4.0 * 3.14159));
    s_y
}

/// Input:
/// Output:
/// Description:
pub fn momentum_density_z(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let h_star = specific_enthalpy(prim.clone(), a_index);
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let b_zcomp = b_z(prim.clone());
    let s_z = (prim.1 * h_star * w * w * prim.4) - ((b_0 * b_zcomp) / (4.0 * 3.14159));
    s_z
}

/// Input:
/// Output:
/// Description:
pub fn total_energy_density(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> f64 {
    let h_star = specific_enthalpy(prim.clone(), a_index);
    let w = lorentz_factor(prim.2, prim.3, prim.4);
    let b_0 = b_zero(prim.clone());
    let p_star = total_pressure(prim.clone());
    let en_den = (prim.1 * h_star * w * w) - p_star - ((b_0 * b_0) / (4.0 * 3.14159)) - (prim.1 * w);
    en_den
}

/// Input:
/// Output:
/// Description:
pub fn prim_to_cons(prim: (f64, f64, f64, f64, f64, f64, f64, f64), a_index: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let d = rest_mass_density(prim.clone());
    let s_x = momentum_density_x(prim.clone(), a_index);
    let s_y = momentum_density_y(prim.clone(), a_index);
    let s_z = momentum_density_z(prim.clone(), a_index);
    let tau = total_energy_density(prim.clone(), a_index);
    let cons = (d, s_x, s_y, s_z, tau, prim.5, prim.6, prim.7);
    cons
}