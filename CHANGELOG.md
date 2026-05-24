 ## May 22 2026
 - Began writing a simple nonrelativistic mhd code for example purposes.
 - Two files were written for this code: main.rs and math_func.rs
 - The following functions were written for this code:
   - total pressure
   - total energy
   - prim to cons
   - cons to prim
   - flux
 - Updated the CHANGELOG file.
 
 ## May 21 2026
 - Wrote or updated the following functions for MHD:
   - round float (location: magneto/math_functions.rs)
   - adjugate (location: magneto/math_functions.rs)
   - Newton root finder (location: magneto/math_functions.rs)
   - cons to prim (location: magneto/math_functions.rs)
 
 ## May 19 2026
 - Wrote or updated the following functions for MHD:
   - cons to prim (location: magneto/math_functions.rs)
   - two_determinant (location: magneto/math_functions.rs)
   - three_determinant (location: magneto/math_functions.rs)
 
 ## May 18 2026
 - Wrote or updated the following functions for MHD:
   - prim to cons (location: magneto/math_functions.rs)
   - initial params (location: magneto/main.rs)
   - lorentz factor (location: magneto/math_functions.rs)
   - rest mass density (location: magneto/math_functions.rs)
   - b_zero (location: magneto/math_functions.rs)
   - b_x (location: magneto/math_functions.rs)
   - b_y (location: magneto/math_functions.rs)
   - b_z (location: magneto/math_functions.rs)
   - b_square (location: magneto/math_functions.rs)
   - specific_energy_gas (location: magneto/math_functions.rs)
   - sound_speed (location: magneto/math_functions.rs)
   - specific_energy_total (location: magneto/math_functions.rs)
   - total_pressure (location: magneto/math_functions.rs)
   - specific_enthalpy (location: magneto/math_functions.rs)
   - momentum_density_x (location: magneto/math_functions.rs)
   - momentum_density_y (location: magneto/math_functions.rs)
   - momentum_density_z (location: magneto/math_functions.rs)
   - total_energy_density (location: magneto/math_functions.rs)
 - Wrote or updated the documentation for the following functions.
   - b_x (location: magneto/math_functions.rs)
   - b_y (location: magneto/math_functions.rs)
   - b_z (location: magneto/math_functions.rs)
   - b_square (location: magneto/math_functions.rs)

 ## May 15 2026
 - Paused implementation of the higher order solution.
 - Began implementation of MHD. The following functions were written or modified:
   - initial params (location: magneto/main.rs)
   - lorentz factor (location: magneto/math_functions.rs)
   - b_zero (location: magneto/math_functions.rs)
   - b_x (location: magneto/math_functions.rs)
   - b_y (location: magneto/math_functions.rs)
   - b_square (location: magneto/math_functions.rs)
   - specific_energy_gas (location: magneto/math_functions.rs)
   - sound_speed (location: magneto/math_functions.rs)
   - specific_energy_total (location: magneto/math_functions.rs)
   - total_pressure (location: magneto/math_functions.rs)
   - specific_enthalpy (location: magneto/math_functions.rs)
   - momentum_density_x (location: magneto/math_functions.rs)
   - momentum_density_y (location: magneto/math_functions.rs)
   - total_energy_density (location: magneto/math_functions.rs)
 - Documentation has been written for the following functions:
   - initial params (location: magneto/main.rs)
   - lorentz factor (location: magneto/math_functions.rs)
   - b_zero (location: magneto/math_functions.rs)
 - Seperation of the simulation into multiple files.
 - Created the math_functions.rs file.
 - Updated the CHANGELOG file.
 - Updated the README file.
 
 ## May 7 2026
 - Began implementation of the higher order solution.
 - Wrote the following functions for the rust simulation.
    - tuple min
    - sign
    - rk4 step
    - minmod
    - high order hll flux
 - Finished the write checkpoint function.
 - Created a folder in the magneto folder to hold the checkpoint txt files.
 - Wrote the animation file in python to animate the data from the simulation. This includes the following functions.
    - plot params
    - line num 
    - read txt files
    - animation function to animate each individual frame
 - Updated the CHANGELOG file.
 - Updated the README file.
 
 ## May 6 2026
 - Wrote the following functions.
    - Conserved Vector Build
    - Primitive Vector build from the conserved vector
    - main containing the simulation
    - write checkpoint (incomplete)
 - Finished the L function.
 - Fixed a bug in the godonov flux function.
 - Updated the CHANGELOG file.
 
 ## May 5 2026
 - Continuing to write one dimensional hydrocode functions.
    - tuple max
    - inverse energy density
    - prim to cons
    - cons to prim
    - flux
    - positive eigenvalue
    - negative eigenvalue
    - hll flux
    - godonov flux
    - compute time step
    - L function (incomplete)
- Wrote a Doc String for the following function.
    - tuple max
- Updated the CHANGELOG file.
- Updated the README file.

 ## May 4 2026
 - Began writing functions for a one dimensional hydrocode.
    - specific energy 
    - sound speed 
    - total energy density
    - intial primitive 
- Created the CHANGELOG file.
- Wrote the README description.

 ## May 1 2026
 - Created Repository with README, config, and rust files.