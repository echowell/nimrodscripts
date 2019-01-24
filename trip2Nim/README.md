# TripToNim
TripToNim is a collection of scripts for processing trip3D external magnetic 
field data and converting it into a file formate that NIMROD can read.

## Inputs:
  - trip3D probe_gb.out
## Outputs:
  - nimrod brmpn##.dat

## Key Steps:
  - Read Trip3D probe file files 
  - Fourier Transform data from physical space to configuration space
  - Write NIMROD brmp files 

## Status: 
  - Code runs but needs to be tested.

## Todo:
  - [ ] Test code
  - [ ] Verify sign of B toroidal
  - [ ] Verify sign of phi
  - [ ] Verify that I'm using the correct mode number (n or -n)
  - [ ] Verify FFT normalization (where do I divide by N)

## Future Work: 
  - Possibly use vector potential instead of B Field