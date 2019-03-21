# combineDump
These scripts combine dumpfiles with multiple Fourier modes, into one larger dump file.
## Inputs:
  - nimrod h5 dump with each perturbation
## Outputs:
  - dumpgll.00000.h5 with all perturbations

## Key Steps:
  - Read dump files
  - Combine Fourier modes
  - Write new h5 dump files

## Status: 
  - Code runs and works for 2 dumps files with 1 mode each
  - Needs more testing
  - Slow

## Todo:
  - [x] Make sure I have h5 py
  - [x] Write dumpTime
  - [x] Combine and write keff
  - [x] Write seams
  - [x] Write file attributes
  - [x] Combine and write rblocks
  - [x] Write code to combine two dumpfiles
  - [x] Test the output with one fourier mode in each dump file 
  - [ ] Test the output with multiple fourier mode in a dump file 
  - [ ] Speed up copy (can I avoid the looping)