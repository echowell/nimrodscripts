# biotSavart
These scripts calculate the edge magnetic field from external rmp coils
by integrating the Biot Savart law 

## Inputs:
  - Coil specification
  - Coil currents
  - NIMROD boundary node locations

## Outputs:
  - Edge magnetic field decomposed by toroidal mode

## Classes:
  - A generic coil class
  - Specific ways to specify C-coil and I-coil like perturbations
  - Biot-Savart Integrator

## Key Steps:
  - Read in node locations and generate a list of r' points
  - Generate coils
  - Integrate Biot-Savart Equation
  - FFT B-Field to get Fourier decomposition
  - Write bmn file

## Status: 
  - Can setup simple circular coils
  - ODE integrator up and running
  - Tested with circular coil and point at the center of coil
  - Simple test of C coil with n=0 and n=1 pert at 0,0,0
  - Code runs, n=1 c coil fields look reasonable
  - Much of the functionality needs to be moved to functions, to clean up

## Todo:
  - [x] Read node locations at generate r' points
  - [x] Write coil class
  - [x] Write initialization for circular coils
  - [x] Write initialization for C-coils  
  - [ ] Write initialization for I-coils
  - [x] Write Biot-Savart Integrator
  - [x] Test integrator for planar coils
  - [x] Write FFT and write functionality (can copy lots from trip2NIM)
  - [ ] Plot coils