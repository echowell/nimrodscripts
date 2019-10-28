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

def coil_calc(filePath,poloidal_coils,coil_theta,coil_tilt,distance_coil,coil_r):

  rzfile = filePath + 'nimrod_bdry_rz.txt'
  baseFileName = "brmpn"
  fileExt = ".dat"

  coilList=[]
  rmag_axis=1.7682
  zmag_axis=0.0
  phiPlanes = 4 #n=0->2
  segmentsPerCoil = 50
  segmentsPerCoil = 25
  baseCurrent = 1.0
  nPert = 1

########### begin code #########
  nodeRZ = np.loadtxt(rzfile,comments='%',delimiter=',', skiprows=1)

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


  for ii in range(6):
    phi = np.pi * ii/3.0
    thisCurrent = baseCurrent * np.cos(phi*nPert)
    for jj in range(4): 
      theta = coil_theta[jj] * 2.0 * np.pi
  #    theta = np.pi * jj/2.0
      if (jj%2 == 0):
        this_m_current = -thisCurrent 
      else:
        this_m_current= thisCurrent

      r0 = rmag_axis + distance_coil * np.cos(theta)
      z0 = zmag_axis + distance_coil * np.sin(theta)
      x0 = r0 * np.cos(phi)
      y0 = r0 * np.sin(phi)
      tx = 0.0 
      ty = coil_tilt[jj] * 2.0 * np.pi
      tz = phi

      thisCoil = cc.coil(this_m_current,segmentsPerCoil)
      thisCoil.planarCoil(x0,y0,z0,coil_r,tx,ty,tz)
      coilList.append(thisCoil)

  for iNode in range(nodeXYZ.shape[1]):
      print("Calculating node: " + str(iNode))
      for iPhi in range(nodeXYZ.shape[2]):
          print("Calculating plane: " + str(iPhi))
          sys.stdout.flush()
          bXYZ=np.zeros(3)
          for iCoil in coilList:
              bXYZ[:]+=bs.intCoil(iCoil,nodeXYZ[:,iNode,iPhi])
          phi = 2.0*np.pi*iPhi/phiPlanes
  ### transform to bRZPhi
  # bRZPhi accounts for the negative in Bphi due to rzphi coordinates
          bRZPhi[0,iNode,iPhi] = bXYZ[0]*np.cos(phi)+bXYZ[1]*np.sin(phi)
          bRZPhi[1,iNode,iPhi] = bXYZ[2]
          bRZPhi[2,iNode,iPhi] = bXYZ[0]*np.sin(phi)-bXYZ[1]*np.cos(phi)

  bRPhase=np.fft.fft(bRZPhi[0,:,:],axis=1)/(float(phiPlanes))
  bZPhase=np.fft.fft(bRZPhi[1,:,:],axis=1)/(float(phiPlanes))
  bPhiPhase=np.fft.fft(bRZPhi[2,:,:],axis=1)/(float(phiPlanes))

### write brmp files
  if (phiPlanes % 2 == 0): #even
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

#    print(rho_2,rho_3,rho_4)
    bmn1=interp1d(mrGrid[1,:,1],bmn[:,9])
    bmn2=interp1d(mrGrid[1,:,1],bmn[:,8])
    bmn3=interp1d(mrGrid[1,:,1],bmn[:,7])
    bmn4=interp1d(mrGrid[1,:,1],bmn[:,6])
    bres1=bmn1(rho_2) #evaluamte m=1 at rho=2 surface
    bres2=bmn2(rho_2)
    bres3=bmn3(rho_3) 
    bres4=bmn4(rho_4)
#    print(bres1,bres2,bres3,bres4)
    small2=(1.0e-10**2)
    value=(bres3**2 + bres4**2)/(bres2**2+small2)
    return value