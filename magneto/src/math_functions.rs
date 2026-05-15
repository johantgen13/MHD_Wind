// This is the rust file containing the math functions for the simulation.
// 
// Author: Brayden JoHantgen
// Last Update: 5/15/2026

/// Input:
///     vx: The x-component of the fluid velocity.
///     vy: The y-component of the fluid velocity.
/// Output:
///     w: The Lorentz factor
/// Description:
///     This function produces the Lorentz factor of the fluid.
///     It uses the equation: (1 - v^2)^(-1/2)
fn lorentz_factor(vx: f64, vy: f64) -> f64 {
    let w = 1.0 / (1 - (vx * vx + vy * vy)).sqrt();
    w
}