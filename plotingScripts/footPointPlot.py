#!/usr/local/bin/python3
''' This script reads a surfcross.txt file and plots the magnetic footpoint'''

import numpy as np
import matplotlib.pyplot as plt
import math as m
import os

def rzToS(r,z):
    tol=1e-4
    r1=1.016
    z1=-1.223
    r2=1.153
    z2=-1.363
    r3=1.372
    z3=-1.363
    r4=1.372
    r5=1.682
    z4=-1.25
    m2=(z2-z1)/(r2-r1)
    b2=z1-m2*r1
    l2 =m.sqrt((z2-z1)**2 +(r2-r1)**2)
    if ( (abs(r-r1)<tol) and (z>=z1)):
        '''on the first segment'''
        s=0-z
    elif (abs(z-m2*r-b2)<tol):
#        print ('slant')
        s=-z1 + m.sqrt((z-z1)**2+(r-r1)**2)
    elif ((abs(z-z3)<tol) and (r>r2-tol) and (r<r3+tol)):
    #    print('bottom')
        s=-z1+l2+(r-r2)
    elif((abs(z-z4)<tol) and (r>r4-tol) and (r<r5+tol)):
        s=-z1+l2+(r3-r2)+(r-r4)
    else:
        print("bad point")
        s=100
    #print(r,z,s)
    return s

def sortData(rawData,pssData):
    lastLine = -1
    numLines=0
    maxHits=1
    thisHits=0
    for ii in range(rawData.shape[0]):
        if rawData[ii,3]==lastLine:
            thisHits+=1
            if thisHits>maxHits:
                maxHits=thisHits
        else:
            lastLine=rawData[ii,3]
            numLines+=1
            thisHits=1
#        msg = "This line " + str(rawData[ii,3]) +" thisHits: " + str(thisHits)
#        print(msg)
#    print(numLines)
#    print(maxHits)
    prosData = np.zeros([numLines,maxHits,3])
    pssDict = {}
#    tempPss=np.zeros([int(np.amax(pssData[:,3])),2])
    
    lastLine = -1
    for ii in range(pssData.shape[0]):
        if int(pssData[ii,3])!=lastLine:
            lastLine=int(pssData[ii,3])
            pssDict[lastLine] = pssData[ii,2]

# for physical data phi in [0,2pi] use a large phi to hide data from plot
    prosData[:,:,1]=1000.
    lastLine=-1
    lineIndex=-1
    hitIndex=-1
    for ii in range(rawData.shape[0]):
        if rawData[ii,3]==lastLine:
            hitIndex+=1
        else:
            lineIndex+=1
            hitIndex=0
            lastLine=rawData[ii,3]
        prosData[lineIndex,hitIndex,0]=rzToS(rawData[ii,0],rawData[ii,1])
        prosData[lineIndex,hitIndex,1]=rawData[ii,2]*360/(2*np.pi)
        prosData[lineIndex,hitIndex,2]=pssDict[lastLine]
    return prosData
    
homeDir = os.environ['HOME']
relDir = "/SCRATCH/166439/footpoint_03300_q104/lphi4/S7Pr1e2/"
fileName = "surfcross0000000.txt"
pssFileName = "nimfl0000000.dat"
fullFileName = homeDir+relDir+fileName
pssFullFileName = homeDir+relDir+pssFileName
plotTitle = "Vacuum Magnetic Footpoint with n=0-21 no sol"

pltt0=0.
plttf=360.
plts0=1.15
pltsf=1.3
plts0=0
pltsf=3
minLength=10
vMax=1e4
rawData = np.loadtxt(fullFileName)
pssData = np.loadtxt(pssFullFileName)
prosData = sortData(rawData,pssData)


for ii in range(prosData.shape[0]):
    if prosData[ii,0,2]<minLength: continue
    plt.scatter(prosData[ii,:,1],prosData[ii,:,0],c=np.log10(prosData[ii,:,2]),cmap='tab20c',vmin=np.log10(minLength),vmax=np.log10(vMax),s=1)

plt.axis([pltt0,plttf,plts0,pltsf])
plt.xlabel('Toroidal Angle (deg)')
plt.ylabel('Distance along wall (m)')
plt.title(plotTitle)
plt.show()


#print(prosData)
#print(fullFileName)