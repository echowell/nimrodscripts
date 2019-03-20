# TripToNim
A collection of scripts for processing trip3D external magnetic field data and converting it into a format that NIMROD reads.

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
  - [x] Test code
  - [x] Verify sign of B toroidal
  - [x] Verify sign of phi
  - [x] Verify that I'm using the correct mode number (n or -n)
  - [x] Verify FFT normalization (where do I divide by N)

## Possible Future Work: 
  - Use vector potential instead of B Field