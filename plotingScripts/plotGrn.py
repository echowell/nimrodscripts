#!/usr/local/bin/python3
''' This script reads a .grn file and plots the seperatrix and wall'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
import os
import sys as sys

homeDir = os.environ['HOME']
relDir = "/SCRATCH/166439/03300_vac_eq/complexconj_rmp/"
fileName = "sol.grn"




# Begin actual code
fullFilePath = homeDir + relDir + fileName
xyCon = [] # list of contours
it = 0
with open(fullFilePath, 'r') as grnFile:
  while(True):
    thisLine = grnFile.readline()
    thisWords = thisLine.split()
    if(len(thisWords)==0): break
    if (thisWords[0] != 'x'+str(it)):
      sys.exit("Expecting x" +str(it) + " read " + thisWords[0])
    nSep = int(thisWords[1])
    thisCon = np.zeros([nSep,2])
    thisLine = grnFile.readline()
    thisWords = thisLine.split()
    for ix, xx in enumerate(thisWords):
      thisCon[ix,0]= float(xx)
    thisLine = grnFile.readline()
    thisWords = thisLine.split()
    if (thisWords[0] != 'y'+str(it)):
      sys.exit("Expecting y" +str(it) + " read " + thisWords[0])
    if (int(thisWords[1])!=nSep): sys.exit("nSep x != nSep y")
    thisLine = grnFile.readline()
    thisWords = thisLine.split()
    for iy, yy in enumerate(thisWords):
      thisCon[iy,1]= float(yy)
    xyCon.append(thisCon)
    it+=1
    if (it>5):
      break

for ii, iCon in enumerate(xyCon):
  plt.scatter(iCon[:,0],iCon[:,1],s=1, label = "Contour :" + str(ii))
  plt.legend()
plt.show()

plt.scatter(xyCon[1][:,0],xyCon[1][:,1])
plt.show()