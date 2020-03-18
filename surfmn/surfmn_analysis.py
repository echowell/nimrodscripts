#!/usr/bin/env python3
#
#
# Do neat stuff

import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import scipy.interpolate as interp
import surfmnstep 
import matplotlib.colors as mcolors

def get_m_index(m, m_max):
  return m+m_max

def surfmn_plot(mr,fmn,ffac,levels,vmax,m_range,q,title,cbar_ticks,reset=False):
  if reset:
    levels=301
    vmax=np.amax(fmn)*ffac
    levels=np.linspace(0,vmax,301)
    cbar_ticks = np.linspace(0,vmax,11)

  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.set_cmap('nipy_spectral')
  conf=plt.contourf(mr[0,:,:],mr[1,:,:],fmn*ffac,levels=levels,vmax=vmax)
  plt.plot(m_range,q,c='w')
  plt.title(title,fontsize=16)
  plt.ylabel(r'<r>',fontsize=16)
  plt.xlabel('m',fontsize=16)
  cbar=fig.colorbar(conf,ticks=cbar_ticks)
  plt.show()

def radial_plot(mr,fmn,ffac,mlist,m_max,qlist,qlabel,title,ylabel):
  colorlist = list(mcolors.TABLEAU_COLORS)
  fig = plt.figure(figsize=(6,5))
  ax= fig.add_subplot(111)
  for im,this_m in enumerate(mlist):
    this_i = get_m_index(this_m,m_max)
    plt_lbl = "m = " + str(this_m)
    tc=colorlist[im]
    ax.plot(mr[1,:,1],fmn[:,this_i]*ffac, color=tc, label=plt_lbl)
  for iq,this_q in enumerate(qlist):
    this_rho = psi_of_q(this_q) 
    this_lbl = qlabel[iq]
    tc=colorlist[iq]
    ax.axvline(this_rho,ls=':',color=tc, label=this_lbl)

  ax.axhline(0,ls='-',c='k')
  ax.legend(loc=0)
  plt.title(title)
  plt.xlabel(r'<r>')
  plt.ylabel(ylabel)
  plt.tight_layout()
  plt.show()
