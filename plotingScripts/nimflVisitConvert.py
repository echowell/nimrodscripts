#!/usr/local/bin/python3
''' This script reads a reads a nimfl.dat file and converts it to a plot 3D'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
import os

homeDir = os.environ['HOME']
relDir = "/SCRATCH/166439/footpoint_03300_q104/lphi5/S7Pr1e2_nimfl2/"
fileName = "nimfl0000000.dat"
plot3Dfile = "nimfl000000.3D"

fullFile = homeDir+relDir+fileName
fullOutFile = homeDir+relDir+plot3Dfile

pssData = np.loadtxt(fullFile)
print(pssData.shape)

f = open(fullOutFile, "wt")
f.write("x y z value\n")
for ii in range(pssData.shape[0]):
  f.write("%g %g %g %g\n" % (pssData[ii,0],pssData[ii,1],0.0,pssData[ii,2]))
f.close()