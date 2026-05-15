 ## May 15 2026
 - Paused implementation of the higher order solution.
 - Began implementation of MHD. The following functions were written or modified:
   - initial params (location: main.rs)
   - lorentz factor (location: math_functions.rs)
 - Documentation has been written for the following functions:
   - initial params (location: main.rs)
   - lorentz factor (location: math_functions.rs)
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