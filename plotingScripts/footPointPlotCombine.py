#!/usr/bin/env python3

''' This script reads a surfcross.txt file and plots the magnetic footpoint'''

import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob
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
        s=-100
    #print(r,z,s)
    return s

def shiftphi(phiin, shift=30):
    # Fix add 30 to be consistent
    phi=phiin*360/(2*np.pi)+shift
    if phi >360:
        phi = phi-360
    return phi

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
def test(value, min, max):
    return  np.logical_and(value > min, value < max)

def main(step,**kwargs):
    DIRNAME = "temp_dir"
    SURFPRE = "surfcross"
    SURFPOST = ".txt"

    surffile = SURFPRE + str(step).zfill(7)+SURFPOST
    print(surffile)
    surffile = "surfcross0000100.txt"
    print(surffile)
    data_list = []
    for item in glob.glob(DIRNAME+"*"):
        os.chdir(item)
        print(os.getcwd())
        file = surffile

        try:
            raw = np.loadtxt(file)
            pros = np.array( [[rzToS(raw[i,0],raw[i,1]),
                              shiftphi(raw[i,2]),
                              raw[i,3]]
                              for i in range(1,raw.shape[0]) ] )
            data_list.append(pros)
        except:
            print(f"file {file} not found")
            raise IOError

        os.chdir('../')

    plasma_model = "nonlinear"

    pltt0=0.
    plttf=360.
    plts0=1.15
    plts0=1.0
    pltsf=1.3
    pltsf=1.3
    minLength=70.
    vMax=1e5

    min=40
    max=100
    soff=1.285
    for data in data_list:
        mask = test(data[:,2],min,max)
        #for ii in range(data.shape[0]):
    #    if prosData[ii,0,2]<minLength: continue
    #cmap='prism'
        #for ii in range(data.shape[0]):
    #        if data[ii,2] > min and data[ii,2]<max:
        plt.scatter(data[mask,1],data[mask,0]-soff,c=np.power(data[mask,2],1),vmin=np.power(min,1),vmax=np.power(max,1),s=1)
    plt.vlines(80,plts0-soff,pltsf-soff,linestyles='-',linewidths=5)
    plt.hlines(1.285-soff,0,360,linestyles='-.',linewidths=1)
    plt.text(150,1.2875, "Inner strike point",fontsize=12)
    plt.axis([pltt0,plttf,plts0-soff,pltsf-soff])
    plt.xlabel('Toroidal Angle (deg)')
    plt.ylabel('Approx Distance from ISP (m)')
    plt.colorbar()
    #plt.title(plotTitle)
    plt.show()

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='NIMFL divertor plate runner')
        parser.add_argument('dumpstep',
                            action='store',
                            type=int,
                            help="dumpstep" )

        args = vars(parser.parse_args())

        main(args['dumpstep'], **args)
