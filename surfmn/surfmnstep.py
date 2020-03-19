#!/usr/bin/env python3
#
#surfmnstep is a container class for storing data useful for analyzing surfmn
#data
#
#
import h5py
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class SurfmnStep:
  def __init__(self,surfmn_file,dumpfile,stepnumber,steptime):
    self.surfmn_file=surfmn_file
    self.dumpfile=dumpfile
    self.step=stepnumber
    self.time=steptime
    self.mmax=None
    self.surfmn_data=False #set to true surfmn file has been read
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
  def read_surfmn(self):
    '''Reads a surfmn h5 file
       Stores profiles in np.arrays
       Calculates psimn from bmn
       sets mmax (no consistancy checks across bmn)
    '''
    self.nlist=[]
    self.ndict={}
    self.bmnlist=[]
    self.psimnlist=[]
    index=0
    with h5py.File(self.surfmn_file,'r') as fc:
      self.surfmn_data=True
      self.rho = np.array(fc['rho'][:])
      profs = fc['prof'][:]
      self.vprime = -1.0*np.array(profs[0])
      self.q = np.array(profs[1])
      self.mr=np.array(fc['surfmnGrid'][:])
      self.psi_q = interp.interp1d(self.q,self.rho)
      self.qmin=np.amin(self.q)
      self.qmax=np.amax(self.q)
      for key in fc.keys():
        if key.startswith("Bmn"):
          self.nlist.append(int(key[3:]))
          self.ndict[int(key[3:])]=index
          index+=1
          thisbmn=np.array(fc[key][:])
          if not(self.mmax):
            self.mmax=int((thisbmn.shape[1]-1)/2)
          self.bmnlist.append(thisbmn)
          thispsi=thisbmn
          for ix,iv in enumerate(self.vprime):
            thispsi[ix,:]=iv*thispsi[ix,:]
          self.psimnlist.append(thispsi)
    if not(self.surfmn_data):
      print("Can not read surfmn file in read_surfmn")
      raise
  def plot_surfmn(self,field,nn,**kwargs): #(*args,**kwargs)
    ''' To do, I want to clean up the api'''
# make sure surfmn has been read and preform a consistancy check
# on the input
    if nn<1:
      print("nn must be positive by convention in plot_surfmn")
      raise ValueError
    if not(self.surfmn_data):
      self.read_surfmn()
    ndex=self.ndict[nn]
    showplot=False
# Create a new figure if an axis object is not included in kwargs
    if 'axis' in kwargs:
      ax = kwargs['axis']
    else:
      showplot=True
      if 'figsize' in kwargs:
        fig = plt.figure(figsize=kwargs['figsize'])
      else:
        fig = plt.figure(figsize=(6,5))
      ax=fig.add_subplot(111)
# Determine which field to plot and set titles and labels
    if (field=='b'):
      fmn=np.copy(self.bmnlist[ndex])
      title=f"$b$ for n={int(nn)} at {self.time*1000:.3f}ms"
    elif (field=='psi'):
      fmn=np.copy(self.psimnlist[ndex])
      title=f"$\psi$ for n={int(nn)} at {self.time*1000:.3f}ms"
    if "scale" in kwargs:
      fmn*=kwargs["scale"]
# set contour levels, i could generalize this further if needed
    levels=kwargs.get("levels",301)
    vmax=kwargs.get("vmax",np.amax(fmn))
    levels=np.linspace(0,vmax,301)
    cbar_ticks=np.linspace(0,vmax,11)
# Update plot based on keys in kwargs
    title=kwargs.get("title",title)
    ylabel=kwargs.get("ylabel",r'<r> [m]')#todo this should be rho
    xlabel=kwargs.get("xlabel","Poloidal Mode number m")
    fontsize=kwargs.get("fontsize",16)
    showplot=kwargs.get("plot",showplot) #orverides showplot logic
# set up mrange()
    mrange=np.linspace(self.qmin,self.qmax)
#create the surfmn plot
    plt.set_cmap('nipy_spectral')
    conf=plt.contourf(self.mr[0,:,:],self.mr[1,:,:],fmn,levels=levels,vmax=vmax)
    plt.plot(nn*mrange,self.psi_q(mrange),c='w')
    plt.title(title,fontsize=fontsize)
    plt.ylabel(r'<r>',fontsize=fontsize)
    plt.xlabel('m',fontsize=fontsize)
    cbar=fig.colorbar(conf,ticks=cbar_ticks)
    if showplot:
      plt.show()

    # Do awesome stuff
  def plot_radial(self,field,nn,mlist,**kwargs):
    ''' To do, I want to clean up the api'''
    if not(self.surfmn_data):
      self.read_surfmn()
    if nn<1:
      print("nn must be positive by convention in plot_radial")
      raise ValueError
    ndex=self.ndict[nn]
    showplot=False
# Create a new figure if an axis object is not included in kwargs
    if 'axis' in kwargs:
      ax = kwargs['axis']
    else:
      showplot=True
      if 'figsize' in kwargs:
        fig = plt.figure(figsize=kwargs['figsize'])
      else:
        fig = plt.figure(figsize=(6,5))
      ax=fig.add_subplot(111)
    if (field=='b'):
      fmn=np.copy(self.bmnlist[ndex])
      title=f"$b_m$ for n={int(nn)} at {self.time*1000:.3f}ms"
      ylable="b"
    elif (field=='psi'):
      fmn=np.copy(self.psimnlist[ndex])
      title=f"$\psi_m$ for n={int(nn)} at {self.time*1000:.3f}ms"
      ylabel=f"$\psi$"
    if "scale" in kwargs:
      fmn*=kwargs["scale"]
    colorlist = list(mcolors.TABLEAU_COLORS)
# Update plot based on keys in kwargs
    title=kwargs.get("title",title)
    ylabel=kwargs.get("ylabel",ylabel)#todo this should be rho
    xlabel=kwargs.get("xlabel",r'<r> [m]')
    fontsize=kwargs.get("fontsize",16)
    showplot=kwargs.get("plot",showplot) #orverides showplot logic
    qlist=kwargs.get("qlist",[])
    for im,this_m in enumerate(mlist):
      this_i = self.get_m_index(this_m)
      if this_i!= None:
        mlbl = "m = " + str(this_m)
        tc=colorlist[im]
        ax.plot(self.mr[1,:,1],fmn[:,this_i], color=tc, label=mlbl)
    for iq,qq in enumerate(qlist):
      try:
        irho = self.psi_q(qq)
        qlbl = f"q = {qq:.2f}"
        tc=colorlist[iq]
        ax.axvline(irho,ls=':',color=tc, label=qlbl)
      except:
        print(f"q={qq:.2f} is not in the domain")
    ax.axhline(0,ls='-',c='k')
    ax.legend(loc=0)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    if showplot:
      plt.show()

  def get_resonance(self,field,nn,mm):
    ''' Evaluate the resonant component of a field at the given resonces'''
    if nn<1:
      print("nn must be positive by convention in get_resonance")
      raise ValueError
    ndex=self.ndict[nn]
    mdex=self.get_m_index(mm)
    if ndex==None:
      print(f"{nn} is not a n number in surfmn file")
      raise ValueError
    if mdex==None:
      print(f"{mm} is not an m number in surfmn file")
      raise ValueError
    qres=mm/nn
    if qres<self.qmin or qres>self.qmax:
      print(qres,self.qmin,self.qmax)
      print(f"The q value {qres} is not resonant")
      raise ValueError
    if (field=='b'):
      resfield=interp.interp1d(self.rho,self.bmnlist[ndex][:,mdex])
    elif (field=='psi'):
      resfield=interp.interp1d(self.rho,self.psimnlist[ndex][:,mdex])
    else:
      print(f"Field {field} is not reconized")
      raise
    return resfield(self.psi_q(qres))
