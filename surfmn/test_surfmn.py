#!/usr/bin/env python3

import fsa_surfmn as surfmn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

case =2
if case==0:
  dumpfile = "dumpgll.10000.h5"
  nimrodin = "nimrod.in"
  time=0.541
elif case ==1:
  dumpfile = "dumpgll.28000.h5"
  nimrodin = "nimrod.in"
  time=1.784
elif case ==2:
  dumpfile = "dumpgll.66000.h5"
  nimrodin = "nimrod.in"
  time=5.000
else:
  print(f"Case {case} not recognized")
  raise

surf=surfmn.fsasurfmn(dumpfile,nimrodin)
rzo=np.array([1.768,-0.018831,0.0])
fargs={}
eqflag=1
mmax=15
n_max=5
fargs['ifour']=list(range(1,n_max+1))
fargs['mmax']=mmax
surf.calculate(rzo=rzo,nsurf=150,eqflag=eqflag,fargs=fargs)


q1list=[-4,-3,-2]
m1list=list(range(-4,0))
q2list=[-4,-3,-2.5,-2,-1.5]
m2list=list(range(-8,0))
q3list=[-3,-2.33, -2,-1.67,-1.33]
m3list=list(range(-12,-3))
q4list=[-3,-2,-1.75,-1.5,-1.25]
m4list=list(range(-15,-6))
q5list=[-3,-2,-1.8,-1.6,-1.4,-1.2]
m5list=list(range(-15,-6))

qlists=[q1list,q2list,q3list,q4list,q5list]
mlists=[m1list,m2list,m3list,m4list,m5list]

#fig=plt.figure(figsize=(8,8))
#ax=fig.add_subplot(111)
#plt.plot(surf.rhon,-1.0*surf.q)
#plt.title(r"Safety factor",fontsize=16)
#plt.ylabel(r'|q|',fontsize=16)
#plt.xlabel(r'$\rho_N$',fontsize=16)
#plt.tight_layout()
#plt.show()

for ii,imode in enumerate(fargs['ifour']):
  fig = plt.figure(figsize=(10,8))
  ax=fig.add_subplot(111)
  title=f"$\psi$(n={int(imode)}) at {time:.3f}ms"
  ylabel=f"$\psi_m$ [mWb]"
  colorlist = list(mcolors.TABLEAU_COLORS)
  xlabel=r'$\rho_N$'
  fontsize=18
  for im,this_m in enumerate(mlists[ii]):
    this_i = surf.get_m_index(this_m)
    if this_i!= None:
      mlbl = "m = " + str(this_m)
      tc=colorlist[im%len(colorlist)]
      ax.plot(surf.rhon,surf.bmn[ii,this_i,:]*1000, color=tc, label=mlbl)
  for iq,qq in enumerate(qlists[ii]):
    try:
      irho = surf.get_rho_q(qq)
      qlbl = f"q = {qq:.2f}"
      tc=colorlist[iq]
      ax.axvline(irho,ls=':',color=tc, label=qlbl)
    except:
      print(f"q={qq:.2f} is not in the domain")
  ax.axhline(0,ls='-',c='k')
  ax.legend(loc=0,frameon=True,fontsize=fontsize)
  plt.title(title,fontsize=fontsize)
  plt.xlabel(xlabel,fontsize=fontsize)
  plt.ylabel(ylabel,fontsize=fontsize)
  plt.tight_layout()
  plt.show()

  if True:
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111)
    # Determine which field to plot and set titles and labels
    #fmn=np.copy(surf.psimnlist[ndex])
    title=f"$\psi$(n={int(imode)}) at {time:.3f}ms"
    # set contour levels, i could generalize this further if needed
    levels=301
    vmax=np.amax(surf.bmn[ii,:,:])*1000
    levels=np.linspace(0,vmax,301)
    cbar_ticks=np.linspace(0,vmax,11)
    # Update plot based on keys in kwargs
    xlabel="Poloidal Mode Number m"
    fontsize=18
    # set up mrange()
    qmin=np.amin(surf.q)
    qmax=np.amax(surf.q)
    mrange=np.linspace(qmin,qmax)
    #create the surfmn plot
    plt.set_cmap('nipy_spectral')
    m=range(-mmax,mmax+1)
    mv, rv = np.meshgrid(m, surf.rhon, sparse=False, indexing='ij')
    conf=plt.contourf(mv,rv,np.clip(surf.bmn[ii,:,:]*1000,0,None),levels=levels,vmax=vmax)
    plt.plot(imode*mrange,surf.get_rho_q(mrange),c='w')
    plt.title(title,fontsize=fontsize)
    plt.ylabel(r'$\rho_N$',fontsize=fontsize)
    plt.xlabel(xlabel,fontsize=fontsize)
    cbar=fig.colorbar(conf,ticks=cbar_ticks)
    plt.xlim(-mmax,mmax)
    plt.show()
