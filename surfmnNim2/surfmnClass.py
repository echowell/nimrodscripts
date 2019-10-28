#!/usr/local/bin/python3

import sys, getopt
import os
import xySliceClass as xy
import numpy as np
import math
import matplotlib.pyplot as plt

class SurfmnClass:
  ''' Base class for surfmn calculations '''
  def __init__(self,argv):
    ''' initialize an surfmn instance with default inputs'''
    self.pd=-1
    self.mx=-1
    self.my=-1
    self.mrange=10
    self.xyFile='/home/research/ehowell/SCRATCH/166439/03300_q104_flowtesting/n1_run1/orginal_exb/300000/xy_slice.bin'
    self.homeDir = os.environ['HOME']
    self.parseInput(argv)

  def parseInput(self,argv):
    ''' This function parses the command line input for surfmn2.py'''
    try:
      opts, args = getopt.getopt(argv,"hs:x:y:p:m:")
    except getopt.GetoptError:
      print('syrfmn2.py -s <xy_slice> -x <mx> -y <my> -p <pd> -m <mrange>')
      sys.exit(2)
    for opt,arg in opts:
      if opt == '-h':
        print('syrfmn2.py -s <xy_slice> -x <mx> -y <my> -p <pd> -m <mrange>')
        sys.exit(0)
      if opt == '-s':
        self.xyFile = arg
      if opt =='-x': 
        self.mx = int(arg)
      if opt =='-y': 
        self.my = int(arg)
      if opt =='-p': 
        self.pd = int(arg)
      if opt =='-m': 
        self.mrange = int(arg)
    if self.mx<0:
      sys.exit("specify mx with -mx")
    if self.my<0:
      sys.exit("specify my with -my")
    if self.pd<0:
      sys.exit("specify pd with -pd")

  def calcBmn(self,xyFields):
    self.Z0 = xyFields.Z[0,0] # assumes flux aligned mesh
    self.R0 = xyFields.R[0,0]
    rMinor = np.sqrt((xyFields.R-self.R0)**2+(xyFields.Z-self.Z0)**2)

    #calculate length along poloidal extent of poloidal flux contour
    s=np.zeros_like(xyFields.R)
    ds=np.zeros_like(xyFields.R)
    for i in range(len(s[:,0])):
      for j in range(len(s[0,:])-1):
        ds[i,j]=math.sqrt((xyFields.R[i,j+1]-xyFields.R[i,j])**2+(xyFields.Z[i,j+1]-xyFields.Z[i,j])**2)
        s[i,j+1]=s[i,j]+ds[i,j]

    #calculate equilibrium poloidal field for [r,pol] locations
    self.B0P=np.sqrt(xyFields.B0R**2+xyFields.B0Z**2)
    if xyFields.J0T[0,0]>0: self.B0P=-self.B0P

    dqds=(-xyFields.B0T)/(2*math.pi*self.B0P*xyFields.R**2)
    q1=np.trapz(dqds,s,axis=1)
    jac=q1[:,None]*xyFields.R**3*self.B0P/(-xyFields.B0T) 

    #calculate straight-field line theta (PEST coordinate, Jim's derivation)
    theta_str=np.zeros_like(xyFields.R)
    dtheta_str=np.zeros_like(xyFields.R)
    for i in range(len(theta_str[:,0])):
      for j in range(len(theta_str[0,:])-1):
        theta_str[i,j+1]=theta_str[i,j]+1./(q1[i]+1.0e-11)*(ds[i,j]*(-xyFields.B0T[i,j])/(self.B0P[i,j]*xyFields.R[i,j]**2))
        dtheta_str[i,j]=1./(q1[i]+1.0e-11)*(ds[i,j]*(-xyFields.B0T[i,j])/(self.B0P[i,j]*xyFields.R[i,j]**2))
#      for j in range(len(theta_str[0,:])): #ECH This shift is needed
#        theta_str[i,j]=theta_str[i,j]-theta_str[i,xyFields.pd]


    dFSAAdth=2*math.pi*jac
    FSArea=np.trapz(dFSAAdth,theta_str,axis=1)
    drhodth=2*math.pi*jac*rMinor/(FSArea[:,None]+1.0e-11)
    self.rho1=np.trapz(drhodth,theta_str,axis=1)
# todo
    rholcfs=self.rho1[int(len(self.rho1)*.75)] # this should be ~ (1-mvac)/mx   
    self.rho1=self.rho1/rholcfs

    for i in range(len(q1)):
      mid2=(q1[i]+2.)*(q1[i+1]+2.)
      if mid2<0:
        self.irho12=i
      mid3=(q1[i]+3.)*(q1[i+1]+3.)
      if mid3<0:
        self.irho13=i
      mid4=(q1[i]+4.)*(q1[i+1]+4.)
      if mid4<0:
        self.irho14=i
      mid5=(q1[i]+5.)*(q1[i+1]+5.)
      if mid5<0:
        self.irho15=i
      mid6=(q1[i]+6.)*(q1[i+1]+6.)
      if mid6<0:
        self.irho16=i
      mid10=(q1[i]+10.)*(q1[i+1]+10.)
      if mid10<0:
        irho110=i
        break


    mmax=self.mrange
    mmin=-self.mrange
    m=np.linspace(mmin,mmax,mmax-mmin+1)

# ech not sure about dim
    #bcnm should have dim[nm,nx]
    bcnm=np.zeros([mmax-mmin+1,xyFields.R.shape[0]])
    bsnm=np.zeros([mmax-mmin+1,xyFields.R.shape[0]])


    for k in range(bcnm.shape[0]): #loop over m
      dbcnmdth=2*np.pi/(FSArea[:,None]+1e-11)*jac*(xyFields.BRr*np.cos((mmin+k)*theta_str)-xyFields.BRi*np.sin((mmin+k)*theta_str))
      dbsnmdth=2*np.pi/(FSArea[:,None]+1e-11)*jac*(-xyFields.BRr*np.sin((mmin+k)*theta_str)-xyFields.BRi*np.cos((mmin+k)*theta_str))
      bcnm[k,:]=np.trapz(dbcnmdth,theta_str,axis=1)
      bsnm[k,:]=np.trapz(dbsnmdth,theta_str,axis=1)

    self.bnm=np.sqrt(bcnm**2+bsnm**2)

  def plotBmn(self):

    fig,ax=plt.subplots(figsize=(8,6))
    ax.plot(self.rho1,self.bnm[9,:].real,color='m',label='m=-1',lw=3)
    ax.plot(self.rho1,self.bnm[8,:].real,color='r',label='m=-2',lw=3)
    ax.plot(self.rho1,self.bnm[7,:].real,color='b',label='m=-3',lw=3)
#    ax.plot(self.rho1,self.bnm[6,:].real,color='g',label='m=-4',lw=3)
#    ax.plot(self.rho1,self.bnm[5,:].real,color='y',label='m=-5',lw=3)
#    ax.plot(self.rho1,self.bnm[4,:].real,color='lime',label='m=-6',lw=3)
    ax.plot(self.rho1,self.bnm[11,:].real,color='cyan',label='m=1',lw=3)
    ax.axvline(x=self.rho1[self.irho12],lw=3,ls=(0,(3,2)),c='r',label=r'$q=-2$')
    ax.axvline(x=self.rho1[self.irho13],lw=3,ls=(0,(3,2)),c='b',label=r'$q=-3$')
    ax.axvline(x=self.rho1[self.irho14],lw=3,ls=(0,(3,2)),c='g',label=r'$q=-4$')
    ax.axvline(x=self.rho1[self.irho15],lw=3,ls=(0,(3,2)),c='y',label=r'$q=-5$')
    ax.axvline(x=self.rho1[self.irho16],lw=3,ls=(0,(3,2)),c='lime',label=r'$q=-5$')
    
    ax.legend(loc=1,ncol=2,fontsize=14)
    
    ax.yaxis.major.formatter._useMathText = True
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.locator_params(axis='x',nbins=5)

    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$<r_m>$',fontsize=24)
    ax.set_ylabel(r'$B_m$',fontsize=24)
    plt.show()