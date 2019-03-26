#!/usr/local/bin/python3
# coil class
# 

import numpy as np

class coil:
    coilType=''
    current=0.0
    segments=0
    xyz=np.zeros([0])
    def __init__(self,type,current,segment):
        self.coilType=type
        self.current=current
        self.segments=segment
        self.xyz=np.zeros([3,self.segments+1])
    def planarCoil(self,x0,y0,z0,r,t1,t2):
    #Initalize an coil around r0, with radius r
    #    t1 and t2 are the rotation angles are the x0 and y0 axis

        cost1=np.cos(t1)
        cost2=np.cos(t2)
        sint1=np.sin(t1)
        sint2=np.sin(t2)
        rotMat = np.zeros([3,3])
        rotMat[0,0] = cost2
        rotMat[0,1] = sint1 * sint2
        rotMat[0,2] = cost1 * sint2
        rotMat[1,0] = 0.0
        rotMat[1,1] = cost1
        rotMat[1,2] = -sint1
        rotMat[2,0] = -sint2
        rotMat[2,1] = sint1 * cost2
        rotMat[2,2] = cost1 * cost2

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
