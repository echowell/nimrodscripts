#!/usr/local/bin/python3
''' This script reads a reads a nimfl.dat file and converts it to a plot 3D'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
import os

homeDir = os.environ['HOME']
relDir = "/SCRATCH/166439/03300_2_equilbria/19091201_vac_lphi5_fp_deg50/"
#relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_rmp_cfl_b/200000_50deg/"
#relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_nolinear_restart/58000_50deg/"
relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_nolinear_fresh/26000/"
fileName = "nimfl0026000.dat"
plot3Dfile = "nimfl0026000.3D"

fullFile = homeDir+relDir+fileName
fullOutFile = homeDir+relDir+plot3Dfile

pssData = np.loadtxt(fullFile)
print(pssData.shape)

f = open(fullOutFile, "wt")
f.write("x y z value\n")
for ii in range(pssData.shape[0]):
  f.write("%g %g %g %g\n" % (pssData[ii,0],pssData[ii,1],0.0,pssData[ii,2]))
f.close()