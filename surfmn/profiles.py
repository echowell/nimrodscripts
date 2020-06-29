#!/usr/bin/env python3
#
#profiles is a class for calculating 1D profiles
# using the flux surface integration
#
#
import f90nml
from eval_nimrod import *
from field_class import *
from fsa import *
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d,splev,UnivariateSpline
import os
import h5py
import sys
import numpy as np
import pickle

class Profiles:
  def __init__(self,dumpfile,nimrodin):
    self.dumpfile=dumpfile
    self.nimrodin=nimrodin
    self.setprofiles=False
    self.psin=np.empty([1])
    self.psi=np.empty([1])
    self.rhon=np.empty([1])
    self.bigr=np.empty([1])
    self.nd=np.empty([1])
    self.p=np.empty([1])
    self.pe=np.empty([1])
    self.ti=np.empty([1])
    self.te=np.empty([1])
    self.q=np.empty([1])
    self.jpar=np.empty([1])
    self.kpol=np.empty([1])
    self.omegator=np.empty([1])
    self.vprime=np.empty([1])
    self.f=np.empty([1])
    self.invr2=np.empty([1])
  def dump(self,file):
    pickle.dump(self.dumpfile,file)
    pickle.dump(self.nimrodin,file)
    pickle.dump(self.setprofiles,file)
    if self.setprofiles==True:
      pickle.dump(self.psin,file)
      pickle.dump(self.psi,file)
      pickle.dump(self.rhon,file)
      pickle.dump(self.bigr,file)
      pickle.dump(self.nd,file)
      pickle.dump(self.p,file)
      pickle.dump(self.pe,file)
      pickle.dump(self.ti,file)
      pickle.dump(self.te,file)
      pickle.dump(self.q,file)
      pickle.dump(self.jpar,file)
      pickle.dump(self.kpol,file)
      pickle.dump(self.omegator,file)
      pickle.dump(self.vprime,file)
      pickle.dump(self.f,file)
      pickle.dump(self.invr2,file)
  def load(self,file):
    self.dumpfile=pickle.load(file)
    print(self.dumpfile)
    self.nimrodin=pickle.load(file)
    print(self.nimrodin)
    self.setprofiles=pickle.load(file)
    print(self.setprofiles)
    if self.setprofiles==True:
      self.psin=pickle.load(file)
      self.psi=pickle.load(file)
      self.rhon=pickle.load(file)
      self.bigr=pickle.load(file)
      self.nd=pickle.load(file)
      self.p=pickle.load(file)
      self.pe=pickle.load(file)
      self.ti=pickle.load(file)
      self.te=pickle.load(file)
      self.q=pickle.load(file)
      self.jpar=pickle.load(file)
      self.kpol=pickle.load(file)
      self.omegator=pickle.load(file)
      self.vprime=pickle.load(file)
      self.f=pickle.load(file)
      self.invr2=pickle.load(file)
  def fsaint(self,rzc,dy,evalnimrod,isurf):
    '''
    Integrand for fluxsurface integration
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    dy(0)=dl/deta or d eta/dl
    dy(1)=dr/deta or dr/dl
    dy(2)=1/bdgth
    dy(3)=dq
    '''
    n0 = np.array(evalnimrod.eval_field('n',rzc,dmode=0,eq=3))
    v0 = np.array(evalnimrod.eval_field('v', rzc, dmode=0, eq=3)) #eq + n0
    p0 = np.array(evalnimrod.eval_field('p',rzc,dmode=0,eq=3))
    pe0 = np.array(evalnimrod.eval_field('pe',rzc,dmode=0,eq=3))
    ti0 = np.array(evalnimrod.eval_field('ti',rzc,dmode=0,eq=3))
    te0 = np.array(evalnimrod.eval_field('te',rzc,dmode=0,eq=3))
    b0 = np.array(evalnimrod.eval_field('b', rzc, dmode=0, eq=3))
    j0 = np.array(evalnimrod.eval_field('j',rzc,dmode=0,eq=3))

    # set up useful working variables
    modbp = np.sqrt(np.dot(b0[0:2],b0[0:2]))
    modb=np.sqrt(np.dot(b0,b0))
    jpar=np.dot(j0,b0)/modb #

    _small = 1.0e-12
    #update integrands
    dy[4]=n0*dy[2]
    dy[5]=np.dot(v0[0:2],b0[0:2])*dy[2]/max(modbp**2,_small)#Kpol
    dy[6]=v0[2]*dy[2]/rzc[0] #omega_t
    dy[7]=p0*dy[2] #p0
    dy[8]=pe0*dy[2] #pe
    dy[9]=ti0*dy[2] #ti
    dy[10]=te0*dy[2]#te
    dy[11]=jpar*dy[2] #jpar
    dy[12]=rzc[0]*b0[2]*dy[2] #F
    dy[13]=dy[2]/rzc[0]**2 #1/R^2
    #update bigr
    self.bigr[isurf]=max(self.bigr[isurf],rzc[0])
    return dy

  def get_omega_exb(self, n, rhon=None):
    omega_exb=n*(4*np.pi**2*self.q*self.kpol/self.vprime-self.omegator)
    if rhon==None:
      return omega_exb
    return interp1d(self.rhon, omega_exb, kind='cubic')(rhon)

  def get_rho_q(self,q):
    try:
      return interp1d(self.q,self.rhon, kind='cubic',fill_value="extrapolate")(q)
    except:
      print(f"The safety factor {q} is not it the domain")
      raise

  def get_field_rho(self,field,rhon):
    try:
      return interp1d(self.rhon,field, kind='cubic')(rhon)
    except:
      print(f"Problem evaluitng field at rhon={rhon}")
      raise

  def calculate(self,rzo=None,rzx=None,**kwargs):
    mi=kwargs.get("mi",3.3435860e-27)
    qe=kwargs.get("qe",1.609e-19)
    nsurf=kwargs.get("nsurf",150)

    self.setprofiles=True
    self.bigr=np.full([nsurf],-np.inf,dtype=np.float64)

    evalnimrod=EvalNimrod(self.dumpfile,fieldlist='nvptbj')
    dvar, yvar, contours = FSA(evalnimrod, rzo, self.fsaint, 10, nsurf=nsurf, \
                                depvar='eta', dpow=0.5, rzx=rzx,normalize=True)

    iend=-1
    temp = np.where(yvar[0,:]==np.nan)
    print(temp)
    for ix,ivar in enumerate(yvar[0,:]):
      print(ix,ivar)
    print(np.isnan(yvar[:,:iend]).any())
    while np.isnan(yvar[:,iend]).any():
        iend -= 1
    iend += yvar.shape[1]+1
    rhomin=np.min(dvar[1,:iend])
    rhomax=np.max(dvar[1,:iend])
    self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)
    #dvars
    tempbigr=np.copy(self.bigr)
    self.psin=interp1d(dvar[1,:iend], dvar[0,:iend], kind='cubic')(self.rhon)
    self.psi=interp1d(dvar[1,:iend], dvar[2,:iend], kind='cubic')(self.rhon)
    self.bigr=interp1d(dvar[1,:iend], tempbigr[:iend], kind='cubic')(self.rhon)
    self.vprime=np.pi*2*interp1d(dvar[1,:iend], dvar[6,:iend], kind='cubic')(self.rhon)
    self.q=interp1d(dvar[1,:iend], dvar[7,:iend], kind='cubic')(self.rhon)

  # yvars
    self.nd=interp1d(dvar[1,:iend], yvar[0,:iend], kind='cubic')(self.rhon)
    self.kpol=interp1d(dvar[1,:iend], yvar[1,:iend], kind='cubic')(self.rhon)
    self.omegator=interp1d(dvar[1,:iend], yvar[2,:iend], kind='cubic')(self.rhon)
    self.p=interp1d(dvar[1,:iend], yvar[3,:iend], kind='cubic')(self.rhon)
    self.pe=interp1d(dvar[1,:iend], yvar[4,:iend], kind='cubic')(self.rhon)
    self.ti=interp1d(dvar[1,:iend], yvar[5,:iend], kind='cubic')(self.rhon)
    self.te=interp1d(dvar[1,:iend], yvar[6,:iend], kind='cubic')(self.rhon)
    self.jpar=interp1d(dvar[1,:iend], yvar[7,:iend], kind='cubic')(self.rhon)
    self.f=interp1d(dvar[1,:iend], yvar[8,:iend], kind='cubic')(self.rhon)
    self.invr2=interp1d(dvar[1,:iend], yvar[9,:iend], kind='cubic')(self.rhon)
