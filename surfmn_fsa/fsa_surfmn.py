#!/usr/bin/env python3
#
# fsa_surfmn is a class for calculating psi_mn and b_mn
# note that the differs from previous versions, psimn was previous called bmn
#
# b_mn is calcualted to match the DIII-D surfmn calculation
# M.J. Schaffer NF 48 (2008), 024004 eqn A.15
#
# the value of b_mn depends on the corrodinate system.
# see eg., J Park Pop 15 (2008)
# the value of psi_mn is more robust
# for any set of Pest-like corrodinate systems psi_mn will be independant
# of the raidal dimension
# furthermore for any flux alligned coordinates (e.g. including Hamada, Boozer)
# the value of psi_mn evalatuted at the rational surface will be the same.
# 
# For pest-corrodinates B_mn is psi_mn / A * factors(2 pi)
# where A is the flux-surface surface area calculated using the Jacobian
# J_surfmn = qR^2/F * R^2 B_p
#
import f90nml
import eval_comp_nimrod as ceval
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
    self.psimn=np.empty([1])
    self.psicmn=np.empty([1])
    self.psismn=np.empty([1])
    self.calc_surfmn = False
    self.s_surfmn = np.empty([1])
    self.phasemn = None
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
      pickle.dump(self.psimn,file)
      pickle.dump(self.psicmn,file)
      pickle.dump(self.psismn,file)
    pickle.dump(self.calc_surfmn,file)
    if self.calc_surfmn == True:
      pickle.dump(self.s_surfmn,file)

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
      self.psimn=pickle.load(file)
      self.psicmn=pickle.load(file)
      self.psismn=pickle.load(file)
    try:
      self.calc_surfmn=pickle.load(file)
    except:
      self.calc_surfmn=False
    if self.calc_surfmn == True:
      self.s_surfmn=pickle.load(file)

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
    dy(4)= dummy
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
    dy(5)=J_smn/J_pest * 1/bdgth
    '''
    #self.mmax=fargs.get("mmax")
    b0=np.array(cfsa.get_b0(evalnimrod,rzc,flag))
    b = evalnimrod.eval_field('b', rzc, dmode=0)
    rr =rzc[0]
    q=self.qpsi(fargs.get('psi'))
    jac=rr*q/b0[2]
    dy[4] = dy[2]/jac #dtheta
    dy[5] = rr * np.sqrt(b0[0]**2 + b0[1]**2) * dy[2] #add surfmn J
    for ii, im in enumerate(self.ifour):
      oset = ii * (4*self.mmax+1)
      ndx = evalnimrod.modelist.index(im)
      reb=np.real(b[:,ndx]) #im+1]) #im+1 is true for nonlinear only
      imb=np.imag(b[:,ndx]) #im+1]) #todo test nonlinear
      rBePsi=rr*(reb[1]*b0[0]-reb[0]*b0[1])
      iBePsi=rr*(imb[1]*b0[0]-imb[0]*b0[1])
      for im in range(self.mmax):
        nmth=-(self.mmax-im)*y[4] #negative m theta
        pmth=(im+1)*y[4] #positive m theta
        dy[6+im+oset]=(rBePsi*np.cos(nmth)-iBePsi*np.sin(nmth))*dy[2]
        dy[7+self.mmax+im+oset]=(rBePsi*np.cos(pmth)-iBePsi*np.sin(pmth))*dy[2]
        dy[7+2*self.mmax+im+oset]=-(rBePsi*np.sin(nmth)+iBePsi*np.cos(nmth))*dy[2]
        dy[7+3*self.mmax+im+oset]=-(rBePsi*np.sin(pmth)+iBePsi*np.cos(pmth))*dy[2]
      dy[6+self.mmax+oset]=rBePsi*dy[2]
    return dy


  def get_rho_q(self,q):
    try:
      return interp1d(self.q, self.rhon,
                      kind='cubic', fill_value="extrapolate")(q)
    except:
      print(f"The safety factor {q} is not it the domain")
      raise

  def get_psi_q(self,q):
    try:
      return interp1d(self.q, self.psi,
                      kind='cubic', fill_value="extrapolate")(q)
    except:
      print(f"The safety factor {q} is not it the domain")
      raise

  def get_field_rho(self,field,rhon):
    try:
      return interp1d(self.rhon, field, kind='cubic')(rhon)
    except:
      print(f"Problem evaluitng field at rhon={rhon}")
      raise

  def compute_phase(self):
    if self.phasemn is None:
      self.phasemn = np.arctan2(self.psismn, self.psicmn)


  def calculate(self,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
    self.ifour=fargs.get("ifour")
    self.mmax=fargs['mmax']
    self.nfour=len(fargs['ifour'])
    self.setprofiles=True

    #first call to fsa is to calcualte q
    cevalnimrod=ceval.EvalCompNimrod(self.dumpfile,fieldlist='nvptbj')
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
    neq=2+self.nfour*(4*self.mmax+1)
    dvar,yvar,contours = cfsa.FSA(cevalnimrod, rzo, self.surfmn_int, neq, \
      nsurf=nsurf,depvar='eta', dpow=0.5,rzx=rzx,flag=eqflag,normalize=False,\
      fargs=fargs)

    iend=-1
    while np.isnan(yvar[:,iend]).any():
        iend -= 1
    iend += yvar.shape[1]+1

    psimn=np.zeros([self.nfour,2*self.mmax+1,iend])
    psicmn=np.zeros([self.nfour,2*self.mmax+1,iend])
    psismn=np.zeros([self.nfour,2*self.mmax+1,iend])
    for ii in range(self.nfour):
      oset = ii * (4*self.mmax+1)
      psicmn[ii,:,:]= yvar[2+oset:2*self.mmax+3+oset,:iend]*(np.pi*2.0)
      psismn[ii,0:self.mmax,:]=\
        yvar[3+2*self.mmax+oset:3+3*self.mmax+oset,:iend]*(np.pi*2.0)
      psismn[ii,self.mmax+1:2*self.mmax+1,:]=\
        yvar[3+3*self.mmax+oset:3+4*self.mmax+oset,:iend]*(np.pi*2.0)
    psimn=np.sqrt(np.square(psicmn)+np.square(psismn))

    s_surfmn = yvar[1,:iend]*np.pi*2
    
    rhomin=np.min(dvar[1,:iend])
    rhomax=np.max(dvar[1,:iend])
    self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)
    #dvars
    self.psin=interp1d(dvar[1,:iend], dvar[0,:iend], kind='cubic')(self.rhon)
    self.psi=interp1d(dvar[1,:iend], dvar[2,:iend], kind='cubic')(self.rhon)
    self.vprime=np.pi*2*interp1d(dvar[1,:iend], dvar[6,:iend], kind='cubic')(self.rhon)
    self.q=interp1d(dvar[1,:iend], dvar[7,:iend], kind='cubic')(self.rhon)

    self.psicmn=interp1d(dvar[1,:iend],psicmn, kind='cubic')(self.rhon)
    self.psismn=interp1d(dvar[1,:iend],psismn, kind='cubic')(self.rhon)
    self.psimn =interp1d(dvar[1,:iend],psimn , kind='cubic')(self.rhon)

    self.s_surfmn = interp1d(dvar[1,:iend],s_surfmn, kind='cubic')(self.rhon)
    self.calc_surfmn = True

  def get_resonance(self,nn,mm,b_flag="psi"):
      ''' Evaluate the resonant component of a b or psi at the given q=m/n'''
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
      if b_flag=="psi":
        resfield=interp1d(self.rhon,self.psimn[ndex,mdex,:])
      elif b_flag=="bmn":
        if self.calc_surfmn==False:
          raise ValueError("Can not calculate b_mn, no S_surfmn")
        else:
          resfield=interp1d(self.rhon,self.psimn[ndex,mdex,:]/self.s_surfmn*2e4)
      else:
        raise ValueError("b_flag not reconized")

      return resfield(self.get_rho_q(qres))

  def get_phase(self,nn,mm,):
      ''' Evaluate the phase of psi at the given resonces'''
      if self.phasemn is None: self.compute_phase()
      if nn<1:
        print("nn must be positive by convention in get_phase")
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
      return interp1d(self.rhon,self.phasemn[ndex,mdex,:])(self.get_rho_q(qres))

  def plot(self,pargs={}):
    if "field" in pargs:
      field=pargs["field"]
    else:
      field='p'
    qmin = np.amin(self.q)
    qmax = np.amax(self.q)
    if (qmin < 0.0 and qmax < 0.0):
      self.signq = -1
    elif (qmin > 0.0 and qmax > 0.0):
      self.signq = +1
    else:
      self.signq = None
    if field=='b':
      if self.calc_surfmn==True:
        for im,imode in enumerate(self.ifour):
          self.plot_radial_bmn(im,imode,pargs)
          self.plot_surfmn_bmn(im,imode,pargs)
      else:
        print("S no caluted")
        raise ValueError
    else:    
      for im,imode in enumerate(self.ifour):
        self.plot_radial(im,imode,pargs)
        self.plot_surfmn(im,imode,pargs)

  def plot_radial(self,ii,imode,pargs={}):

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    if abs(self.time - 0.0005413792981815724) < 1.0e-6:
        title=f"Peak Pulse (0.5ms)\n $\psi$(n={int(imode)})"
    elif abs(self.time - 0.00500030382482843) < 1.0e-6:
        title=f"Early Slow Growth (5ms)\n $\psi$(n={int(imode)})"
    elif abs(self.time - 0.002969027824826107) < 1.0e-6:
      title=f"Earlier in Time (3ms)\n $\psi$(n={int(imode)})"
    elif abs(self.time - 0.008385763824834448) < 1.0e-6:
      title=f"Late Slow Growth (8ms)\n $\psi$(n={int(imode)})"
    else:
        title=f"$\psi$(n={int(imode)}) at {self.time*1000:.3f}ms"
    ylabel=f"$\psi_m$ [mWb]"
    colorlist = list(mcolors.TABLEAU_COLORS)
    xlabel=r'$\rho_N$'
    fontsize=24
    legendfont=20
    if imode==1:
      if self.signq == -1:
        mlist=range(-4,0)#mlist=range(-4,1)
      else:
        mlist=range(0,4)#mlist=range(-4,1)
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
        mlbl = "m = " + str(abs(this_m))
        tc=colorlist[im%len(colorlist)]
        ax.plot(self.rhon/rhomax,self.psimn[ii,this_i,:]*1000, color=tc, label=mlbl)
    try:
      qlist=pargs['qlists'][ii]
    except:
      if imode==1:
        if self.signq == -1:
          qlist=range(-4,-1)
        else:
          qlist=range(1,5)
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
    print(qlist)
    print(pargs)
    for iq,qq in enumerate(qlist):
      try:
        irho = self.get_rho_q(qq)
        qlbl = f"q = {abs(qq):.1f}"
        tc=colorlist[iq]
    #    if (qq in [-3,-2]): #aps/paper
    #        ax.axvline(irho/rhomax,ls=':',color=tc, label=qlbl)
        #if (qq in [-3,-2]):
        ax.axvline(irho/rhomax,ls=':',color=tc, label=qlbl)
      except:
        print(f"q={qq:.2f} is not in the domain")
    ax.axhline(0,ls='-',c='k')
    ax.legend(loc=0,
              frameon=True,
              fontsize=legendfont,
              ncol=1,
              handlelength=1,
              handletextpad=0.4)
    plt.title(title,fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.ylabel(ylabel,fontsize=fontsize)
    plt.tight_layout()
    plt.show()
  
  def plot_radial_bmn(self,ii,imode,pargs={}):
    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    title=f"$B$(n={int(imode)}) at {self.time*1000:.3f}ms"
    title=f"$B$(n={int(imode)})"
    ylabel=f"$B_r$ [G]"
    colorlist = list(mcolors.TABLEAU_COLORS)
    xlabel=r'$\rho_N$'
    fontsize=24
    legendfont=20
    if imode==1:
      if self.signq == -1:
        mlist=range(-4,0)
      else:
        mlist=range(1,5)
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
        mlbl = "m = " + str(abs(this_m))
        tc=colorlist[im%len(colorlist)]
        ax.plot(self.rhon / rhomax,
                2 * self.psimn[ii,this_i,:] / np.abs(self.s_surfmn) * 10000,
                color=tc,
                label=mlbl)
    try:
      qlist=pargs['qlists'][ii]
    except:
      if imode==1:
        if self.signq == -1:
          qlist=range(-4,-1)
        else:
          qlist=range(2,5)
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
    print(qlist)
    print(pargs)
    for iq,qq in enumerate(qlist):
      try:
        irho = self.get_rho_q(qq)
        qlbl = f"q = {abs(qq):.1f}"
        tc=colorlist[iq]
        if qq in [2,-2]:
          tc = colorlist[2]
        if qq in [3,-3]:
          tc = colorlist[3]
        if qq in [4,-4]:
          tc = colorlist[4]
    #    if (qq in [-3,-2]): #aps/paper
    #        ax.axvline(irho/rhomax,ls=':',color=tc, label=qlbl)
        #if (qq in [-3,-2]):
        #ax.axvline(irho/rhomax,ls=':',color=tc, label=qlbl)
        ax.axvline(irho/rhomax,ls=':',color=tc)
      except:
        print(f"q={qq:.2f} is not in the domain")
    ax.axhline(0,ls='-',c='k')
    ax.legend(loc=0,
              frameon=True,
              fontsize=legendfont,
              ncol=1,
              handlelength=1,
              handletextpad=0.4)
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
    vmax=np.amax(self.psimn[im,:,:])*1000
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
    rhomax=np.max(self.rhon)
    mv, rv = np.meshgrid(m, self.rhon/rhomax, sparse=False, indexing='ij')
    conf=plt.contourf(mv,rv,np.clip(self.psimn[im,:,:]*1000,0,None),levels=levels,vmax=vmax)
    plt.plot(imode*mrange,self.get_rho_q(mrange)/rhomax,c='w')
    plt.title(title,fontsize=fontsize)
    plt.ylabel(r'$\rho_N$',fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    cbar=fig.colorbar(conf,ticks=cbar_ticks)
    plt.xlim(-self.mmax,self.mmax)
    plt.show()

  def plot_surfmn_bmn(self,im,imode,surfmn,pargs={}):
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111)
    # Set titles and labels
    title=f"$B_r$(n={int(imode)}) at {self.time*1000:.3f}ms"
    title=f"$B_r$(n={int(imode)}) [G]"
    # set contour levels, i could generalize this further if needed
    levels=301
    vmax=2 * np.amax(self.psimn[im,:,:] / np.abs(self.s_surfmn))*10000
    #vmax=4.2
    #vmax=4.8
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
    psimax=np.max(self.psi)
    mv, rv = np.meshgrid(m, self.psi/psimax, sparse=False, indexing='ij')
    conf=plt.contourf(mv,
                      rv,
                      np.clip(2 * self.psimn[im,:,:] / np.abs(self.s_surfmn) * 10000, 0, None),
                      levels=levels,
                      vmax=vmax)
    plt.plot(imode*mrange,self.get_psi_q(mrange)/psimax,c='w')
    #plt.plot(imode*mrange+1,self.get_psi_q(mrange)/psimax,c='k',ls='--')
    plt.plot(imode*mrange+2,self.get_psi_q(mrange)/psimax,c='r',ls='-.')
    #plt.plot(imode*mrange+3,self.get_psi_q(mrange)/psimax,c='k',ls=':')
    plt.title(title,fontsize=fontsize)
    plt.ylabel(r'$\psi_N$',fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    cbar=fig.colorbar(conf,ticks=cbar_ticks)
    plt.xlim(-self.mmax,self.mmax)
    plt.show()

  def get_dumptime(self):
    ''' Open the hdf5 dumpfile read the dump time and dumpstep '''
    with h5py.File(self.dumpfile, 'r') as h5file:
      try:
        self.time=h5file["dumpTime"].attrs['vsTime']
        self.step=int(h5file["dumpTime"].attrs['vsStep'])
      except:
        print(f"Error reading time or step in {self.dumpfile}")
        raise
