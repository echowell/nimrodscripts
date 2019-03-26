#!/usr/local/bin/python3
#
# Input files:
# Ouput file:


import os
import numpy as np
import sys
import coilClass as cc
homeDir = os.environ['HOME']


nodeRZ=np.zeros([2,1])
coilList=[]

coilOne = cc.coil('planar',1e3,100)
coilList.append(coilOne)
coilList[0].planarCoil(1,0,0,1,np.pi/2.,np.pi/2.)
print(coilList[0].segments)
sys.exit("debug")