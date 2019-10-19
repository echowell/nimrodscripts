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
#        prosData[lineIndex,hitIndex,1]=rawData[ii,2]*360/(2*np.pi)+180
#        prosData[lineIndex,hitIndex,1]=rawData[ii,2]*360/(2*np.pi)
# Fix add 30 to be consistent
        prosData[lineIndex,hitIndex,1]=rawData[ii,2]*360/(2*np.pi)+30
        if prosData[lineIndex,hitIndex,1]>360:
          prosData[lineIndex,hitIndex,1]=prosData[lineIndex,hitIndex,1]-360
        prosData[lineIndex,hitIndex,2]=pssDict[lastLine]
    return prosData

plasma_model = "nonlinear"
  
homeDir = os.environ['HOME']
relDir = "/SCRATCH/166439/03300_vac_eq/normal_rmp_vac5_fpsep2/"
relDir = "/SCRATCH/166439/03300_2_equilbria/19091201_vac_lphi5_fp/"
relDir = "/SCRATCH/166439/03300_2_fl//19091702/lphi5_nolinear_restart/58000/"
#relDir = "/SCRATCH/166439/03300_2_equilbria/19091702_fl/linear/lphi5_2/50000/"
#relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_rmp_cfl_b/" +dump_num +"/"
#relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_nolinear_restart/" +dump_num +"/"
#relDir = "/SCRATCH/166439/03300/normal_rmp_vac5_fpsep2/"
#relDir = "/SCRATCH/166439/03300_vac_eq/complexconj_rmp_vac_fpsep2/"
fileName = "surfcross0058000.txt"
pssFileName = "nimfl0058000.dat"
#fileName = "surfcross0110000.txt"
#pssFileName = "nimfl0110000.dat"


if plasma_model == "vacuum":
  plotTitle = "Vacuum Response Footprint"
  relDir = "/SCRATCH/166439/03300_2_equilbria/19091201_vac_lphi5_fp_deg50/"
  dump_num = "00100"
elif plasma_model == "linear":
  plotTitle = "Linear Response Footprint"
  relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_rmp_cfl_b/200000_50deg/"
  dump_num="200000"
elif plasma_model =="nonlinear":
  plotTitle = "Nonlinear Response Footprint"
  relDir = "/SCRATCH/166439/03300_2_fl/19091702/lphi5_nolinear_fresh/32000/"
  dump_num = "32000"
else:
  plotTitle = "Vacuum Response Footprint"
  dump_num="58000"  

if len(dump_num)==6:
  fileName = "surfcross0"+dump_num+".txt"
  pssFileName = "nimfl0"+dump_num+".dat"
elif len(dump_num)==5:
  fileName = "surfcross00"+dump_num+".txt"
  pssFileName = "nimfl00"+dump_num+".dat"

#fileName = "surfcross0000100.txt"
#pssFileName = "nimfl0000100.dat"
fullFileName = homeDir+relDir+fileName
pssFullFileName = homeDir+relDir+pssFileName

pltt0=0.
plttf=360.
plts0=1.15
pltsf=1.3
minLength=70.
vMax=1e5
rawData = np.loadtxt(fullFileName)
pssData = np.loadtxt(pssFullFileName)
prosData = sortData(rawData,pssData)


for ii in range(prosData.shape[0]):
    if prosData[ii,0,2]<minLength: continue
    plt.scatter(prosData[ii,:,1],prosData[ii,:,0],c=np.log10(prosData[ii,:,2]),cmap='prism',vmin=np.log10(minLength),vmax=np.log10(vMax),s=1)
plt.vlines(80,plts0,pltsf,linestyles='dotted',linewidths=1)
plt.hlines(1.285,0,360,linestyles='-.',linewidths=1)
plt.text(150,1.2875, "Inner strike point",fontsize=12)
plt.axis([pltt0,plttf,plts0,pltsf])
plt.xlabel('Toroidal Angle (deg)')
plt.ylabel('Distance along wall (m)')
plt.title(plotTitle)
plt.show()


#print(prosData)
#print(fullFileName)
