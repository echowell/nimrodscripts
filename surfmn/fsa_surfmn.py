#!/usr/bin/env python3
#
#profiles is a class for calculating 1D profiles
# using the flux surface integration
#
#
import f90nml
import eval_comp_nimrod as ceval
#from field_class import *
#import comp_fsa as cfsa
import comp_fsa as cfsa
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d,splev,UnivariateSpline
import os
import h5py
import sys
import numpy as np
import pickle

class fsasurfmn:
  def __init__(self,dumpfile,nimrodin):
    #start with dump file and nimrod.in info
    self.dumpfile=dumpfile
    self.nimrodin=nimrodin
    self.time=None
    self.step=None
    # next include info on how fsa's were performed
    self.mmax=None
    self.ifour=[]
    self.nfour=None
    self.setprofiles=False
    #finally end with the profiles
    self.psin=np.empty([1])
    self.psi=np.empty([1])
    self.rhon=np.empty([1])
    self.q=np.empty([1])
    self.qpsi=np.empty([1])
    self.vprime=np.empty([1])
    self.bmn=np.empty([1])
    self.bcmn=np.empty([1])
    self.bsmn=np.empty([1])
  def dump(self,file):
    pickle.dump(self.dumpfile,file)
    pickle.dump(self.nimrodin,file)
    pickle.dump(self.time,file)
    pickle.dump(self.step,file)
    # next include info on how fsa's were performed
    pickle.dump(self.mmax,file)
    pickle.dump(self.ifour,file)
    pickle.dump(self.nfour,file)
    pickle.dump(self.setprofiles,file)
    #finally end with the profiles
    if self.setprofiles==True:
      pickle.dump(self.psin,file)
      pickle.dump(self.psi,file)
      pickle.dump(self.rhon,file)
      pickle.dump(self.q,file)
      pickle.dump(self.vprime,file)
      pickle.dump(self.bmn,file)
      pickle.dump(self.bcmn,file)
      pickle.dump(self.bsmn,file)
  def load(self,file):
    self.dumpfile=pickle.load(file)
    self.nimrodin=pickle.load(file)
    self.time=pickle.load(file)
    self.step=pickle.load(file)
    # next include info on how fsa's were performed
    self.mmax=pickle.load(file)
    self.ifour=pickle.load(file)
    self.nfour=pickle.load(file)
    self.setprofiles=pickle.load(file)
    #finally end with the profiles
    if self.setprofiles==True:
      self.psin=pickle.load(file)
      self.psi=pickle.load(file)
      self.rhon=pickle.load(file)
      self.q=pickle.load(file)
      self.vprime=pickle.load(file)
      self.bmn=pickle.load(file)
      self.bcmn=pickle.load(file)
      self.bsmn=pickle.load(file)

  def get_m_index(self,m):
    '''Return the index for the given m
      Return None if m is out of range'''
    if self.mmax==None:
      write("mmax has not be set in get_m_index")
      raise
    else:
      if m>self.mmax:
        return None
      elif m<-1*self.mmax:
        return None
      else:
        return m+self.mmax
  def dummy_fsa(self,rzc,y,dy,evalnimrod,flag,fargs):
    '''
    Dummy integrand for complex fsa, this is used to get v' and q
    without running a true fsa
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    dy(0)=dl/deta or d eta/dl
    dy(1)=dr/deta or dr/dl
    dy(2)=1/bdgth :v'
    dy(3)=dq: q
    '''
    dy[4]=1.0
    return dy

  def surfmn_int(self,rzc,y,dy,evalnimrod,flag,fargs):
    '''
    Integrand for fluxsurface integration
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    dy(0)=dl/deta or d eta/dl
    dy(1)=dr/deta or dr/dl
    dy(2)=1/bdgth
    dy(3)=dq
    dy(4)=dtheta
    '''

    #self.mmax=fargs.get("mmax")


    b0=np.array(cfsa.get_b0(evalnimrod,rzc,flag))
    b = evalnimrod.eval_field('b', rzc, dmode=0)


    rr =rzc[0]
    q=self.qpsi(fargs.get('psi'))
    jac=rr*q/b0[2]
    dy[4]=dy[2]/jac #dtheta
    for ii, im in enumerate(self.ifour):
      oset = ii * (4*self.mmax+1)
      reb=np.real(b[:,im+1])
      imb=np.imag(b[:,im+1])
      rBePsi=rr*(reb[1]*b0[0]-reb[0]*b0[1])
      iBePsi=rr*(imb[1]*b0[0]-imb[0]*b0[1])
      for im in range(self.mmax):
        nmth=-(self.mmax-im)*y[4] #negative m theta
        pmth=(im+1)*y[4] #positive m theta
        dy[5+im+oset]=(rBePsi*np.cos(nmth)-iBePsi*np.sin(nmth))*dy[2]
        dy[6+self.mmax+im+oset]=(rBePsi*np.cos(pmth)-iBePsi*np.sin(pmth))*dy[2]
        dy[6+2*self.mmax+im+oset]=-(rBePsi*np.sin(nmth)+iBePsi*np.cos(nmth))*dy[2]
        dy[6+3*self.mmax+im+oset]=-(rBePsi*np.sin(pmth)+iBePsi*np.cos(pmth))*dy[2]
      dy[5+self.mmax+oset]=rBePsi*dy[2]
    return dy


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

  def calculate(self,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
    print("in calculate")
    mi=kwargs.get("mi",3.3435860e-27)
    qe=kwargs.get("qe",1.609e-19)
    self.ifour=fargs.get("ifour")
    self.mmax=fargs['mmax']
    self.nfour=len(fargs['ifour'])
    self.setprofiles=True

    #first call to fsa is to calcualte q
    print("before ceval")
    cevalnimrod=ceval.EvalCompNimrod(self.dumpfile,fieldlist='nvptbj')
    print("after ceval")
    dvar, yvar, contours = cfsa.FSA(cevalnimrod, rzo, self.dummy_fsa, 1, \
      nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
      fargs=fargs)


    iend=-1
    while np.isnan(yvar[:,iend]).any():
        iend -= 1
    iend += yvar.shape[1]+1
    #unevaluated interpoate
    self.qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic')

    #second call to fsa is to calcualte b_ms ans psi_mn
    neq=1+self.nfour*(4*self.mmax+1)
    dvar,yvar,contours = cfsa.FSA(cevalnimrod, rzo, self.surfmn_int, neq, \
      nsurf=nsurf,depvar='eta', dpow=0.5,rzx=rzx,flag=eqflag,normalize=False,\
      fargs=fargs)

    iend=-1
    while np.isnan(yvar[:,iend]).any():
        iend -= 1
    iend += yvar.shape[1]+1

    bmn=np.zeros([self.nfour,2*self.mmax+1,iend])
    bcmn=np.zeros([self.nfour,2*self.mmax+1,iend])
    bsmn=np.zeros([self.nfour,2*self.mmax+1,iend])
    for ii in range(self.nfour):
      oset = ii * (4*self.mmax+1)
      bcmn[ii,:,:]= yvar[1+oset:2*self.mmax+2+oset,:iend]*(np.pi*2.0)
      bsmn[ii,0:self.mmax,:]=\
        yvar[2+2*self.mmax+oset:2+3*self.mmax+oset,:iend]*(np.pi*2.0)
      bsmn[ii,self.mmax+1:2*self.mmax+1,:]=\
        yvar[2+3*self.mmax+oset:2+4*self.mmax+oset,:iend]*(np.pi*2.0)
    bmn=np.sqrt(np.square(bcmn)+np.square(bsmn))
    rhomin=np.min(dvar[1,:iend])
    rhomax=np.max(dvar[1,:iend])
    self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)
    #dvars
    self.psin=interp1d(dvar[1,:iend], dvar[0,:iend], kind='cubic')(self.rhon)
    self.psi=interp1d(dvar[1,:iend], dvar[2,:iend], kind='cubic')(self.rhon)
    self.vprime=np.pi*2*interp1d(dvar[1,:iend], dvar[6,:iend], kind='cubic')(self.rhon)
    self.q=interp1d(dvar[1,:iend], dvar[7,:iend], kind='cubic')(self.rhon)

    self.bcmn=interp1d(dvar[1,:iend],bcmn, kind='cubic')(self.rhon)
    self.bsmn=interp1d(dvar[1,:iend],bsmn, kind='cubic')(self.rhon)
    self.bmn =interp1d(dvar[1,:iend],bmn , kind='cubic')(self.rhon)

  def get_resonance(self,field,nn,mm):
      ''' Evaluate the resonant component of a field at the given resonces'''
      if nn<1:
        print("nn must be positive by convention in get_resonance")
        raise ValueError
      ndex=nn-1 #todo check
      mdex=self.get_m_index(mm)
      if ndex==None:
        print(f"{nn} is not a n number in surfmn file")
        raise ValueError
      if mdex==None:
        print(f"{mm} is not an m number in surfmn file")
        raise ValueError
      qres=mm/nn
      #if qres<self.qmin or qres>self.qmax:
      #    print(qres,self.qmin,self.qmax)
      #    print(f"The q value {qres} is not resonant")
      #    raise ValueError
      resfield=interp1d(self.rhon,self.bmn[ndex,mdex,:])
      return resfield(self.get_rho_q(qres))

  def plot(self,pargs={}):
    for im,imode in enumerate(self.ifour):
       self.plot_radial(im,imode,pargs)
       self.plot_surfmn(im,imode,pargs)

  def plot_radial(self,ii,imode,pargs={}):
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111)
    title=f"$\psi$(n={int(imode)}) at {self.time*1000:.3f}ms"
    ylabel=f"$\psi_m$ [mWb]"
    colorlist = list(mcolors.TABLEAU_COLORS)
    xlabel=r'$\rho_N$'
    fontsize=18
    if imode==1:
      mlist=range(-4,1)
    elif imode==2:
      mlist=range(-6,-1)
    else:
      mstart=-2*imode
      mlist=range(mstart,mstart+imode+1)
    if 'mlists' in pargs:
      if ii<len(pargs['mlists'][ii]):
        mlist=pargs['mlists'][ii]
    rhomax=np.max(self.rhon)
    for im,this_m in enumerate(mlist):
      this_i = self.get_m_index(this_m)
      if this_i!= None:
        mlbl = "m = " + str(this_m)
        tc=colorlist[im%len(colorlist)]
        ax.plot(self.rhon/rhomax,self.bmn[ii,this_i,:]*1000, color=tc, label=mlbl)
    try:
      qlist=pargs['qlists'][ii]
    except:
      if imode==1:
        qlist=[-4,-3,-2]
      elif imode==2:
        qlist=[-4,-3,-2.5,-2,-1.5]
      elif imode==3:
        qlist=[-3,-2.33, -2,-1.67,-1.33]
      elif imode==4:
        qlist=[-3,-2,-1.75,-1.5,-1.25]
      elif imode==5:
        qlist=[-3,-2,-1.8,-1.6,-1.4,-1.2]
      else:
        qlist=[-4,-3,-2]

    for iq,qq in enumerate(qlist):
      try:
        irho = self.get_rho_q(qq)
        qlbl = f"q = {qq:.2f}"
        tc=colorlist[iq]
        ax.axvline(irho/rhomax,ls=':',color=tc, label=qlbl)
      except:
        print(f"q={qq:.2f} is not in the domain")
    ax.axhline(0,ls='-',c='k')
    ax.legend(loc=0,frameon=True,fontsize=fontsize)
    plt.title(title,fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.ylabel(ylabel,fontsize=fontsize)
    plt.tight_layout()
    plt.show()

  def plot_surfmn(self,im,imode,surfmn,pargs={}):
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111)
    # Set titles and labels
    title=f"$\psi$(n={int(imode)}) at {self.time*1000:.3f}ms"
    # set contour levels, i could generalize this further if needed
    levels=301
    vmax=np.amax(self.bmn[im,:,:])*1000
    levels=np.linspace(0,vmax,301)
    cbar_ticks=np.linspace(0,vmax,11)
    # Update plot based on keys in kwargs
    xlabel="Poloidal Mode Number m"
    fontsize=18
    # set up mrange()
    qmin=np.amin(self.q)
    qmax=np.amax(self.q)
    mrange=np.linspace(qmin,qmax)
    #create the surfmn plot
    plt.set_cmap('nipy_spectral')
    m=range(-self.mmax,self.mmax+1)
    mv, rv = np.meshgrid(m, self.rhon, sparse=False, indexing='ij')
    conf=plt.contourf(mv,rv,np.clip(self.bmn[im,:,:]*1000,0,None),levels=levels,vmax=vmax)
    plt.plot(imode*mrange,self.get_rho_q(mrange),c='w')
    plt.title(title,fontsize=fontsize)
    plt.ylabel(r'$\rho_N$',fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    cbar=fig.colorbar(conf,ticks=cbar_ticks)
    plt.xlim(-self.mmax,self.mmax)
    plt.show()

  def get_dumptime(self):
    ''' Open the hdf5 dumpfile read the dump time and dumpstep
    '''
    with h5py.File(self.dumpfile, 'r') as h5file:
      try:
        self.time=h5file["dumpTime"].attrs['vsTime']
        self.step=int(h5file["dumpTime"].attrs['vsStep'])
      except:
        print(f"Error reading time or step in {self.dumpfile}")
        raise
