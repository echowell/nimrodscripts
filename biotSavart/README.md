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
  - Just starting out

## Todo:
  - [ ] Read node locations at generate r' points
  - [ ] Write coil class
  - [x] Write initialization for circular coils
  - [ ] Write initialization for C-coils  
  - [ ] Write initialization for I-coils
  - [ ] Write Biot-Savart Integrator
  - [ ] Test integrator for planar coils
  - [ ] Write FFT and write functionality (can copy lots from trip2NIM)
