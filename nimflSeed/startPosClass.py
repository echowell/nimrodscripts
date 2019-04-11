#!/usr/local/bin/python3
''' Base class for generating start position data '''
import numpy as np
import random as ran

class startPosClass:
    numPoints=-1
    randomPhi = False
    phiPlane = 0.0
    rzp = np.zeros([1,3])
    geom = "none"
    def __init__ (self,numPoints,geom,randomPhi,phiPlane):
        ''' Initialize start position class '''
        self.numPoints = numPoints
        self.geom = geom
        self.randomPhi = randomPhi
        self.phiPlane = phiPlane
        self.rzp = np.zeros([self.numPoints,3])
        if (not self.randomPhi):
            self.rzp[:,2]=self.phiPlane
    def d3dlowerRZPhi(self):
        ''' generate a collection of rzp points for the lower d3d divertor '''
        for ii in range(self.numPoints):
            top=-0.7
            bottom= -1.2
            left=1.2
            right=1.5
            self.rzp[ii,0]=ran.uniform(left,right)
            self.rzp[ii,1]=ran.uniform(bottom,top)
            if self.randomPhi:
                self.rzp[ii,2] = ran.uniform(0,2*np.pi)
    def calculateRZPhi(self):
        if self.geom=='d3dlower':
            self.d3dlowerRZPhi()
        elif self.geom=='d3d':
            print("d3d geom")
        else: # do nothing
            print("Geom: ", geom, " not reconized")
    def writeStartPos(self,fileName):
        thisFile = open(fileName,'w')
        thisFile.write(str(self.numPoints)+"\n")
        for jj in range(self.numPoints):
            thisLine = '{: 16.16e}'.format(self.rzp[jj,0]) + ", " 
            thisLine+= '{: 16.16e}'.format(self.rzp[jj,1]) + ", " 
            thisLine+= '{: 16.16e}'.format(self.rzp[jj,2]) + "\n" 
            thisFile.write(thisLine)
        thisFile.close()