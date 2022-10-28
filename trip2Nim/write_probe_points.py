#!/usr/local/bin/python3



import os
import numpy as np
import sys
import matplotlib.pyplot as plt



homeDir = os.environ['HOME']
scratchDir = homeDir + '/SCRATCH/'
wrkDir = 'KSTAR/19118_2950_C1wall/22062201_probeg/'

nimBdry = 'nimrod_bdry_rz.txt'
rzFile = 'probe.points.rz.in'

phiPlanes = 256



############## Start of code ###################
fullNimBdry = scratchDir+wrkDir+nimBdry
fullRzFile = scratchDir+wrkDir+rzFile
nodeRZ = np.loadtxt(fullNimBdry,delimiter=',',skiprows=1)
print(nodeRZ.shape)

rzpts = nodeRZ.shape[0]
probeFile = open(fullRzFile,'w')
probeFile.write("&startptsrz\n")
probeFile.write("\n")
probeFile.write("  nptsrz = " + str(rzpts) + "\n")
probeFile.write("  phinrz = 0 \n")
probeFile.write("  nphinrz = " + str(phiPlanes)+ "\n")
probeFile.write("\n")
probeFile.write("  ptsrz(1:2, 1:" + str(rzpts) + ") = \n")
for ii in range(rzpts):
  thisStr = "   " + str(nodeRZ[ii,0]) + "   " + str(nodeRZ[ii,1]) 
  probeFile.write(thisStr + "\n")
probeFile.write("\n")
probeFile.write("&END")
probeFile.close()
