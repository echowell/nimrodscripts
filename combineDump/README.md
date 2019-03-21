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
  - Code is in early development
  - Copied the script pert merge from nimdevel scripts and rename

## Todo:
  - [ ] Make sure I have h5 py
  - [ ] Write code to combine two dumpfiles
  - [ ] Test the output 
  - [ ] Expand to combing multiple dumpfiles
