#!/usr/local/bin/python3
# coil class
# 

import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class coil:
    coilType=''
    coilId=0
    current=0.0
    segments=0
    xyz=np.zeros([0])
    def __init__(self,current,segment):
        self.current=current
        self.segments=segment
        self.xyz=np.zeros([3,self.segments+1])
    def planarCoil(self,x0,y0,z0,r,tx,ty,tz):
    #Initalize an coil around r0, with radius r
    #    t1 and t2 are the rotation angles are the x0 and y0 axis

        self.coilType='planar'
        costx=np.cos(tx)
        costy=np.cos(ty)
        costz=np.cos(tz)
        sintx=np.sin(tx)
        sinty=np.sin(ty)
        sintz=np.sin(tz)
        rotxMat = np.zeros([3,3])
        rotyMat = np.zeros([3,3])
        rotzMat = np.zeros([3,3])
        
        rotxMat[0,0] = 1.0
        rotxMat[1,1] = costx
        rotxMat[1,2] = -sintx
        rotxMat[2,1] = sintx
        rotxMat[2,2] = costx

        rotyMat[0,0] = costy
        rotyMat[0,2] = sinty
        rotyMat[1,1] = 1.0
        rotyMat[2,0] = -sinty
        rotyMat[2,2] = costy

        rotzMat[0,0] = costz
        rotzMat[0,1] = -sintz
        rotzMat[1,0] = sintz
        rotzMat[1,1] = costz
        rotzMat[2,2] = 1.0

        rotMat = np.matmul(rotzMat,np.matmul(rotyMat,rotxMat))
        thisXYZ=np.zeros([3])
        theta=np.linspace(0,2*np.pi,num=self.segments+1)
        for iT, thisT in enumerate(theta):
            thisXYZ[0]= r * np.cos(thisT)
            thisXYZ[1]= r * np.sin(thisT)
            thisXYZ[2]= 0.0
            self.xyz[:,iT]=np.matmul(rotMat,thisXYZ)
            self.xyz[0,iT]+=x0
            self.xyz[1,iT]+=y0
            self.xyz[2,iT]+=z0
# make sure the coil is periodic
        self.xyz[:,self.segments]=self.xyz[:,0]
    def cCoil(self,cCoilId):
        self.coilType='cCoil'
        self.coilId = cCoilId
        r0 = 3.2 #radius of C Coil
        zTop = 0.8 #C coils have a height of 1.6m
        zBottom = -0.8

# I assume that the coils are composed of 4 segments
# each has equal number of line segments. Test to see if the 
# number of segments is divisble by 4
        if(self.segments%4 !=0):
            sys.exit("c Coil segments must be divisble by 4")
        segsPerSec = int(self.segments/4)

        phiStart = (-1.0 + (cCoilId-1) * 2.0)/(6.) *np.pi
        phiEnd =   ( 1.0 + (cCoilId-1) * 2.0)/(6.) *np.pi
        phi=np.linspace(phiStart,phiEnd,num=segsPerSec+1)
        z  =np.linspace(zBottom, zTop,  num=segsPerSec+1)
        for ii in range(segsPerSec):
#top
            self.xyz[0,ii] = r0 * np.cos(phi[ii])
            self.xyz[1,ii] = r0 * np.sin(phi[ii])
            self.xyz[2,ii] = zTop
#right
            self.xyz[0,ii+segsPerSec]=r0 * np.cos(phiEnd)
            self.xyz[1,ii+segsPerSec]=r0 * np.sin(phiEnd)
            self.xyz[2,ii+segsPerSec]=z[segsPerSec-ii]
#bottom
            self.xyz[0,ii+2*segsPerSec]=r0 * np.cos(phi[segsPerSec-ii])
            self.xyz[1,ii+2*segsPerSec]=r0 * np.sin(phi[segsPerSec-ii])
            self.xyz[2,ii+2*segsPerSec]=zBottom
#left
            self.xyz[0,ii+3*segsPerSec]=r0 * np.cos(phiStart)
            self.xyz[1,ii+3*segsPerSec]=r0 * np.sin(phiStart)
            self.xyz[2,ii+3*segsPerSec]=z[ii]
        self.xyz[:,self.segments]=self.xyz[:,0]
    def plot_coil(self,axis):
      axis.plot(self.xyz[0,:],self.xyz[1,:],self.xyz[2,:])

    def plot_2D_coil(self,axis,coilcolor):
      axis.plot(self.xyz[0,:],self.xyz[2,:],color=coilcolor)
