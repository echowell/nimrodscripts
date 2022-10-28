#!/usr/local/bin/python3
''' This script reads the external magnetic fields from trip3D output files
    and converts that data to a formate readible by NIMROD '''

import sys
sys.path.insert(0, "./")
import os
import tripClass as tc

homeDir = os.environ['HOME']

readDirectory = homeDir + "/SCRATCH/166439/03300_2_equilbria/19091201_probeg/"

#writeDirectory = homeDir + "/SCRATCH/testingjunk/"
#writeDirectory = homeDir + "/SCRATCH/166439/03300_vac_eq/normal_rmp/"
writeDirectory = homeDir + "/SCRATCH/166439/03300_2_equilbria/19091201_probeg/"
readDirectory = homeDir + "/SCRATCH/166439/03300_2_equilbria/19100401_probe_gb/"
writeDirectory = homeDir + "/SCRATCH/166439/03300_2_equilbria/19100401_probe_gb/"

readDirectory = homeDir + "/SCRATCH/KSTAR/19118_2800ms/22032403_probeg/"
writeDirectory = homeDir + "/SCRATCH/KSTAR/19118_2800ms/22032403_rmp_v1/"

readDirectory = homeDir + "/SCRATCH/KSTAR/19118_2950_C1wall/22052401_brmp/"
writeDirectory = homeDir + "/SCRATCH/KSTAR/19118_2950_C1wall/22052401_brmp/"

readDirectory = homeDir + "/SCRATCH/KSTAR/19118_2950_C1wall/22062201_probeg/"
writeDirectory = homeDir + "/SCRATCH/KSTAR/19118_2950_C1wall/22062201_brmp/"

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

rzProbeFile = readDirectory + rzFileSuffix
aProbeFile = readDirectory + aFileSuffix
bProbeFile = readDirectory + bFileSuffix

indexShift=0 #was 5
complexCon = True

tripData = tc.TripClass(rzProbeFile,aProbeFile,bProbeFile,indexShift,complexCon)
tripData.processBFile()
print(tripData.brPhase[0,:])
tripData.writeNimrodBext(writeDirectory,nimrodRmpFile,nimrodRmpSuffix)