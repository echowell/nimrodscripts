#!/usr/local/bin/python3
''' This script randomly generates a collection of seed locations basd on the sol file '''

################################################################################
#  Set up envirment and import modules
################################################################################
import sys
sys.path.insert(0, "./")
import os
pwd = os.getcwd()
homeDir = os.environ['HOME']
import numpy as np
import random as ra
import matplotlib.pyplot as plt
################################################################################
#  User defined input
################################################################################
solDir = "/SCRATCH/166439/03300_vac_eq/complexconj_rmp/"
solFileName = "sol.grn"
solCon = 1 # This is the seperatrix

writeFileName = "start_positions.dat_sep1"
writeDir = "/SCRATCH/166439/03300_vac_eq/start_pos_files/"

nPoints = 3000
sigmaDist = 0.05
phiZero = 0.0
randomPhi = True # Ture

################################################################################
#  Start of code
################################################################################

fullSolFile = homeDir+solDir+solFileName
fullWriteFile = homeDir + writeDir + writeFileName
xyCon = [] # list of contours
it = 0
with open(fullSolFile, 'r') as grnFile:
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

thisCon = xyCon[solCon]
nCon = thisCon.shape[0]
thisRZP = np.zeros([nPoints,3])
for ii in range(nPoints):
  thisIn = ra.randint(0,nCon-1)
  thisR = ra.gauss(0,sigmaDist)
  thisT = ra.uniform(0,np.pi*2.0)
  if randomPhi:
    thisPhi = ra.uniform(0,np.pi*2.0)
  else:
    thisPhi = phi0
  thisRZP[ii,0]=thisCon[thisIn,0]+thisR * np.cos(thisT)
  thisRZP[ii,1]=thisCon[thisIn,1]+thisR * np.sin(thisT)
  thisRZP[ii,2]=thisPhi

for ii, iCon in enumerate(xyCon):
  plt.scatter(iCon[:,0],iCon[:,1],s=1, label = "Contour :" + str(ii))
  plt.legend()
plt.scatter(thisRZP[:,0], thisRZP[:,1],s=1)
plt.show()

with open(fullWriteFile,'w') as thisFile:
  thisFile.write(str(nPoints)+"\n")
  for jj in range(nPoints):
    thisLine = '{: 16.16e}'.format(thisRZP[jj,0]) + ", " 
    thisLine+= '{: 16.16e}'.format(thisRZP[jj,1]) + ", " 
    thisLine+= '{: 16.16e}'.format(thisRZP[jj,2]) + "\n" 
    thisFile.write(thisLine)