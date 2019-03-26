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
filePath = homeDir + '/SCRATCH/174446/fgnimeq_iter5_3_brmp/'
rzfile = filePath + 'nimrod_bdry_rz.txt'
baseFileName = "brmpn"
fileExt = ".dat"


phiPlanes = 4
segmentsPerCoil = 1000
baseCurrent = 1.0
nPert = 1


########### begin code #########
#nimrod RZ node locations
nodeRZ = np.loadtxt(rzfile,comments='%',delimiter=',', skiprows=1)
# allocate xyz and b vectors
print(nodeRZ.shape)

nodeXYZ = np.zeros([3,nodeRZ.shape[0],phiPlanes])
bRZPhi = np.zeros([3,nodeRZ.shape[0],phiPlanes])
bRPhase = np.zeros(1,dtype=np.complex_)
bZPhase = np.zeros(1,dtype=np.complex_)
bPhiPhase = np.zeros(1,dtype=np.complex_)
#convert node locations to xyz coordinates at multiple phi planes
for iPhi in range(phiPlanes):
    sinPhi = np.sin(iPhi*2.0*np.pi/phiPlanes)
    cosPhi = np.cos(iPhi*2.0*np.pi/phiPlanes)
    nodeXYZ[0,:,iPhi]=nodeRZ[:,0]*cosPhi
    nodeXYZ[1,:,iPhi]=nodeRZ[:,0]*sinPhi
    nodeXYZ[2,:,iPhi]=nodeRZ[:,1]

# set up c coils
# in the future I can move this to an function
coilList=[]
for ii in range(6):
    thisCurrent = baseCurrent * np.cos(np.pi*(ii)*nPert/3.0)
    thisCoil = cc.coil(thisCurrent,segmentsPerCoil)
    coilList.append(thisCoil)
    coilList[ii].cCoil(ii+1)

for iNode in range(nodeXYZ.shape[1]):
    print("Calculating node: " + str(iNode))
    for iPhi in range(nodeXYZ.shape[2]):
        print("Calculating plane: " + str(iPhi))
        sys.stdout.flush()
        bXYZ=np.zeros(3)
        for iCoil in coilList:
            bXYZ[:]+=bs.intCoil(iCoil,nodeXYZ[:,iNode,iPhi])
        phi = nodeXYZ[2,iNode,iPhi]
### transform to bRZPhi
# bRZPhi accounts for the negative in Bphi due to rzphi coordinates
        bRZPhi[0,iNode,iPhi] = bXYZ[0]*np.cos(phi)+bXYZ[1]*np.sin(phi)
        bRZPhi[1,iNode,iPhi] = bXYZ[2]
        bRZPhi[2,iNode,iPhi] = bXYZ[0]*np.sin(phi)-bXYZ[1]*np.cos(phi)

bRPhase=np.fft.fft(bRZPhi[0,:,:],axis=1)/(float(phiPlanes))
bZPhase=np.fft.fft(bRZPhi[1,:,:],axis=1)/(float(phiPlanes))
bPhiPhase=np.fft.fft(bRZPhi[2,:,:],axis=1)/(float(phiPlanes))

### write brmp files
if (phiPlanse % 2 == 0): #even
    maxnphi = int(phiPlanes/2)
else: #odd
    maxnphi = int((phiPlanes+1)/2)
for ii in range (maxnphi +1):
    if ii==maxnphi:
        fac=0.5
    else:
        fac=1.0
    print(ii, maxnphi, fac)
    tempFileName = filePath + baseFileName +"{0:0=2d}".format(ii)  + fileExt
    thisFile = open(tempFileName,'w')
    for jj in range(bRPhase.shape[0]):
        thisLine ='{: 16.16e}'.format(fac*bRPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bRPhase[jj,ii].imag) + ", "
        thisLine+='{: 16.16e}'.format(fac*bZPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bZPhase[jj,ii].imag) + ", "
        thisLine+='{: 16.16e}'.format(fac*bPhiPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bPhiPhase[jj,ii].imag) + "\n"
        thisFile.write(thisLine)
    thisFile.close()

#for iCoil in coilList:
#    bvec+=bs.intCoil(iCoil,np.asarray([0.0,0.0,0.0]))

print(coilList[0].segments)
print(np.pi *2)
sys.exit("debug")