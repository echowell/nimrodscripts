#!/usr/local/bin/python3
''' This script reads the external magnetic fields from trip3D output files
    and converts that data to a formate readible by NIMROD '''

import sys
sys.path.insert(0, "./")
import tripClass as tc

readDirectory = "/Users/ehowell/SCRATCH/nimruns/echowell_runs/heatwidthscaling/166439/03300/EF_GRID_18121801/"
writeDirectory = "/Users/ehowell/SCRATCH/testingjunk/"
shotNumber = "166439"
timeSlice = "03300"

rzFileSuffix = "probe.points.rz.in"
aFileSuffix = "probe_ga.out"
bFileSuffix = "probe_gb.out"

nimrodRmpFile = "brmpn"
nimrodRmpSuffix = ".dat"

rzProbeFile = readDirectory + shotNumber + "." + timeSlice + "." + rzFileSuffix
aProbeFile = readDirectory + shotNumber + "." + timeSlice + "." + aFileSuffix
bProbeFile = readDirectory + shotNumber + "." + timeSlice + "." + bFileSuffix

tripData = tc.TripClass(rzProbeFile,aProbeFile,bProbeFile)
tripData.processBFile()
tripData.writeNimrodBext(writeDirectory,nimrodRmpFile,nimrodRmpSuffix)