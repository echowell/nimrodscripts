#!/usr/bin/env python3
#
# Input files:
# Ouput file:


import os
import numpy as np
import sys
import coilClass as cc
import integrateBS as bs
import h5py
from scipy.interpolate import interp1d
import scipy.optimize as opt

class coil_opt:
  rmag_axis=1.7682
  zmag_axis=0.0
  coil_r = 0.6
  distance_coil = 1.2
  delta_theta=0.20
  delta_tilt=.15
  coil_theta = [0, delta_theta, .5, 1-delta_theta]
  coil_tilt = [0.25, .5+delta_tilt, .75, 1-delta_tilt]
  phiPlanes = 4 
  segmentsPerCoil = 50
  baseCurrent = 1.0
  nPert = 1
  npCoil = len(coil_theta)
  nodeXYZ = np.zeros(0)
  coilBRZPhi = np.zeros(0)
  coilList = []
  filePath = ""
  baseFileName = "brmpn"
  fileExt = ".dat"

  def __init__(self,filePath):
    self.filePath =filePath
    rzfile = filePath + 'nimrod_bdry_rz.txt'
########### begin code #########
    nodeRZ = np.loadtxt(rzfile,comments='%',delimiter=',', skiprows=1)
    self.nodeXYZ = np.zeros([3,nodeRZ.shape[0],self.phiPlanes])
    self.coilBRZPhi = np.zeros([3,nodeRZ.shape[0],self.phiPlanes,self.npCoil])
    for iPhi in range(self.phiPlanes):
        sinPhi = np.sin(iPhi*2.0*np.pi/self.phiPlanes)
        cosPhi = np.cos(iPhi*2.0*np.pi/self.phiPlanes)
        self.nodeXYZ[0,:,iPhi]=nodeRZ[:,0]*cosPhi
        self.nodeXYZ[1,:,iPhi]=nodeRZ[:,0]*sinPhi
        self.nodeXYZ[2,:,iPhi]=nodeRZ[:,1]

    for ipcoil in range(self.npCoil):
      self.coilList.append([])

  #convert node locations to xyz coordinates at multiple phi planes

      for ii in range(6):
        phi = np.pi * ii/3.0
        thisCurrent = self.baseCurrent * np.cos(phi*self.nPert)
        theta = self.coil_theta[ipcoil] * 2.0 * np.pi
        
        r0 = self.rmag_axis + self.distance_coil * np.cos(theta)
        z0 = self.zmag_axis + self.distance_coil * np.sin(theta)
        x0 = r0 * np.cos(phi)
        y0 = r0 * np.sin(phi)
        tx = 0.0 
        ty = self.coil_tilt[ipcoil] * 2.0 * np.pi
        tz = phi

        thisCoil = cc.coil(thisCurrent,self.segmentsPerCoil)
        thisCoil.planarCoil(x0,y0,z0,self.coil_r,tx,ty,tz)
        self.coilList[ipcoil].append(thisCoil)

    for iNode in range(self.nodeXYZ.shape[1]):
      print("Calculating node: " + str(iNode))
      for iPhi in range(self.nodeXYZ.shape[2]):
        print("Calculating plane: " + str(iPhi))
        sys.stdout.flush()
        bXYZ=np.zeros(3)
        for ipcoil in range(self.npCoil):
          for iCoil in self.coilList[ipcoil]:
              bXYZ[:]+=bs.intCoil(iCoil,self.nodeXYZ[:,iNode,iPhi])
          phi = 2.0*np.pi*iPhi/self.phiPlanes
  ### transform to bRZPhi
  # bRZPhi accounts for the negative in Bphi due to rzphi coordinates
          self.coilBRZPhi[0,iNode,iPhi,ipcoil] = bXYZ[0]*np.cos(phi)+bXYZ[1]*np.sin(phi)
          self.coilBRZPhi[1,iNode,iPhi,ipcoil] = bXYZ[2]
          self.coilBRZPhi[2,iNode,iPhi,ipcoil] = bXYZ[0]*np.sin(phi)-bXYZ[1]*np.cos(phi)



  def coil_calc(self,coil_currents):
    bRZPhi=np.zeros([self.coilBRZPhi.shape[0],self.coilBRZPhi.shape[1],self.coilBRZPhi.shape[2]])
    bRPhase = np.zeros(1,dtype=np.complex_)
    bZPhase = np.zeros(1,dtype=np.complex_)
    bPhiPhase = np.zeros(1,dtype=np.complex_)
    for ipcoil, icur in enumerate(coil_currents):
      bRZPhi[:,:,:] += self.coilBRZPhi[:,:,:,ipcoil]*icur
    bRPhase=np.fft.fft(bRZPhi[0,:,:],axis=1)/(float(self.phiPlanes))
    bZPhase=np.fft.fft(bRZPhi[1,:,:],axis=1)/(float(self.phiPlanes))
    bPhiPhase=np.fft.fft(bRZPhi[2,:,:],axis=1)/(float(self.phiPlanes))

### write brmp files
    if (self.phiPlanes % 2 == 0): #even
      maxnphi = int(self.phiPlanes/2)
    else: #odd
      maxnphi = int((self.phiPlanes+1)/2)
    for ii in range (maxnphi +1):
      if ii==maxnphi:
        fac=0.5
      else:
        fac=1.0
      tempFileName = self.filePath + self.baseFileName +"{0:0=2d}".format(ii)  + self.fileExt
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

def surfmn_eval(fileName):
  profNames = ['Vprime','q']
  with h5py.File(fileName,'r') as fc:
    for aname, avalue in fc.attrs.items():
      print(aname,avalue)
    mrGrid = fc['surfmnGrid'][:]
    bmn = fc['Bmn001'][:]
    rho = fc['rho'][:]
    profs = fc['prof'][:]
    print(fc.keys())
    
#prof 1 is q
    fq = interp1d(profs[1,:],rho)
    rho_2 = fq(-2.)
    rho_3 = fq(-3)
    rho_4 = fq(-4)
    rho_5 = fq(-4)

#    print(rho_2,rho_3,rho_4)
    bmn1=interp1d(mrGrid[1,:,1],bmn[:,9])
    bmn2=interp1d(mrGrid[1,:,1],bmn[:,8])
    bmn3=interp1d(mrGrid[1,:,1],bmn[:,7])
    bmn4=interp1d(mrGrid[1,:,1],bmn[:,6])
    bmn5=interp1d(mrGrid[1,:,1],bmn[:,5])
    bres1=bmn1(rho_2) #evaluamte m=1 at rho=2 surface
    bres2=bmn2(rho_2)
    bres3=bmn3(rho_3) 
    bres4=bmn4(rho_4)
    bres5=bmn5(rho_5)

#    print(bres1,bres2,bres3,bres4)
    small2=(1.0e-10**2)
    value=(bres1**2 + bres3**2 + bres4**2 + bres5**2)/(bres2**2+small2)
    return value