#!/usr/local/bin/python3
# coil class
# 

# the line integral uses a solution vector
# of length 6
# the first 3 index are the x,y,z locations at (l)
# the next 3 indexs are the values of bx,by,bz

import sys
import numpy as np
import coilClass as cc
from scipy.integrate import odeint
def dBFunc(solVec,l,dxyz,xyzp):
# integrand of Biot Savart for a line segment
    k = 1.0e-7 #mu_0/4pi
    x,y,z,bx,by,bz = solVec
    dx,dy,dz=dxyz
    xp,yp,zp=xyzp
    xVec = np.asarray([xp-x,yp-y,zp-z]) 
    xNorm =np.linalg.norm(xVec)
    dlVec =np.asarray([dx,dy,dz])
    dBvec = np.cross(dlVec,xVec) * k /(xNorm**3)  
    dS = np.asarray([dx,dy,dz,dBvec[0],dBvec[1],dBvec[2]])
    return dS

def intCoil(thisCoil,xyzp):
    bxyzVec = np.zeros(3)
    sol0 = np.zeros(6)
    for ii in range(thisCoil.segments):
        dx = thisCoil.xyz[0,ii+1]-thisCoil.xyz[0,ii]
        dy = thisCoil.xyz[1,ii+1]-thisCoil.xyz[1,ii]
        dz = thisCoil.xyz[2,ii+1]-thisCoil.xyz[2,ii]
        dxyz=np.asarray([dx,dy,dz])
        sol0[0]=thisCoil.xyz[0,ii]
        sol0[1]=thisCoil.xyz[1,ii]
        sol0[2]=thisCoil.xyz[2,ii]
        sol0[3]=bxyzVec[0]
        sol0[4]=bxyzVec[1]
        sol0[5]=bxyzVec[2]
        sol=odeint(dBFunc,sol0,[0.0,1.0],args=(dxyz,xyzp))
        bxyzVec[0]=sol[1,3]
        bxyzVec[1]=sol[1,4]
        bxyzVec[2]=sol[1,5]
    bxyzVec=bxyzVec * thisCoil.current
    return bxyzVec