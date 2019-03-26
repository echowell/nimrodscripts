#!/usr/local/bin/python3
#
# Input files:
# Ouput file:


import os
import numpy as np
import sys
import coilClass as cc
import integrateBS as bs

homeDir = os.environ['HOME']

segmentsPerCoil = 1000
baseCurrent = 1.0
nPert = 1

nodeRZ=np.zeros([2,1])
bvec = np.zeros(3)
coilList=[]
#        print(ii, bxyzVec)
#        sys.stdout.flush()

# set up c coils
# in the future I can move this to an function
for ii in range(6):
    thisCurrent = baseCurrent * np.cos(np.pi*(ii)*nPert/3.0)
    thisCoil = cc.coil(thisCurrent,segmentsPerCoil)
    coilList.append(thisCoil)
    coilList[ii].cCoil(ii+1)

for iCoil in coilList:
    bvec+=bs.intCoil(iCoil,np.asarray([0.0,0.0,0.0]))

print(bvec)
print(coilList[0].segments)
print(np.pi *2)
sys.exit("debug")