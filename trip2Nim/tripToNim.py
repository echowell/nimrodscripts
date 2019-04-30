#!/usr/local/bin/python3
''' This script reads the external magnetic fields from trip3D output files
    and converts that data to a formate readible by NIMROD '''

import sys
sys.path.insert(0, "./")
import os
import tripClass as tc

homeDir = os.environ['HOME']

readDirectory = homeDir + "/SCRATCH/nimruns/echowell_runs/heatwidthscaling/166439/03300/EF_GRID_18121801/"

#writeDirectory = homeDir + "/SCRATCH/testingjunk/"
writeDirectory = homeDir + "/SCRATCH/166439/03300_vac_eq/normal_rmp/"
writeDirectory = homeDir + "/SCRATCH/166439/03300_vac_eq/complexconj_rmp/"

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

indexShift=5
complexCon = True

tripData = tc.TripClass(rzProbeFile,aProbeFile,bProbeFile,indexShift,complexCon)
tripData.processBFile()
print(tripData.brPhase[0,:])
tripData.writeNimrodBext(writeDirectory,nimrodRmpFile,nimrodRmpSuffix)