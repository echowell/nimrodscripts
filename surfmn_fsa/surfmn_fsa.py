#!/usr/bin/env python3
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
from shutil import copy2

import fsa_surfmn as surfmn
import matplotlib.colors as mcolors


''' This is a generic script for analysing a multiple dumpfile/timestep
    using the surfmn fsa scripts '''
def pickle_sort(file):
  return int(file.split('.')[1])

def dump_sort(file):
  return int(file.split('.')[1])

def time_hist(steplist,pargs={}):
    if "field" in pargs:
      f=pargs["field"]
    else:
      f='p'
    if f in ['b']:
      field_flag = 'bmn'
      title=r"$B_r$"
      ytitle=r"$B_r [G]$"
      scale=1
    else:
      field_flag = 'psi'
      title=r"Evolution of $\psi$"
      ytitle=r"$\psi$ [mWb]"
      scale=1e3


    time=np.zeros(len(steplist))
    psi21=np.zeros(len(steplist))
    psi42=np.zeros(len(steplist))
    psi52=np.zeros(len(steplist))
    psi62=np.zeros(len(steplist))
    psi31=np.zeros(len(steplist))
    psi41=np.zeros(len(steplist))
    psi32=np.zeros(len(steplist))
    psi43=np.zeros(len(steplist))
    psi54=np.zeros(len(steplist))
    psi65=np.zeros(len(steplist))
    phase21=np.zeros(len(steplist))
    phase42=np.zeros(len(steplist))
    phase52=np.zeros(len(steplist))
    phase62=np.zeros(len(steplist))
    phase31=np.zeros(len(steplist))
    phase41=np.zeros(len(steplist))
    phase32=np.zeros(len(steplist))
    phase43=np.zeros(len(steplist))
    phase54=np.zeros(len(steplist))
    phase65=np.zeros(len(steplist))
    times=[]
    mlist=[-1,-2,-3,-4]
    qlist=[-2,-3,-4]

    for istep,step in enumerate(steplist):
        print(istep,step.step, step.time)
        time[istep]=step.time
        psi21[istep]=step.get_resonance(1,-2,b_flag=field_flag)*scale
        psi31[istep]=step.get_resonance(1,-3,b_flag=field_flag)*scale
        psi41[istep]=step.get_resonance(1,-4,b_flag=field_flag)*scale
        psi32[istep]=step.get_resonance(2,-3,b_flag=field_flag)*scale
        psi42[istep]=step.get_resonance(2,-4,b_flag=field_flag)*scale
        psi52[istep]=step.get_resonance(2,-5,b_flag=field_flag)*scale
        psi43[istep]=step.get_resonance(3,-4,b_flag=field_flag)*scale
        psi54[istep]=step.get_resonance(4,-5,b_flag=field_flag)*scale
        psi65[istep]=step.get_resonance(5,-6,b_flag=field_flag)*scale
        phase21[istep]=step.get_phase(1,-2)
        phase31[istep]=step.get_phase(1,-3)
        phase41[istep]=step.get_phase(1,-4)

        if step.step==100000:
            eq_q2 = step.get_rho_q(q=-2)
            eq_q3 = step.get_rho_q(q=-3)

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,np.abs(psi21),label="2/1")
    plt.plot(time*1000,np.abs(psi31),label="3/1")
    plt.plot(time*1000,np.abs(psi41),label="4/1")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=20)
    #plt.title(r"$B_r$",fontsize=16)
    #plt.title(r"Vacuum $B_r$",fontsize=16)
    #ax.fill([0,1.0,1.0,0],[0,0,100,100],'gray',alpha=0.2)
    ax.fill([0,1.0,1.0,0],[0,0,40,40],'gray',alpha=0.2)
    plt.title(title,fontsize=24)
    plt.ylabel(ytitle,fontsize=24)
    plt.xlabel(r't [ms]',fontsize=24)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,np.abs(psi21),label="2/1")
    plt.plot(time*1000,np.abs(psi31),label="3/1")
    plt.plot(time*1000,np.abs(psi41),label="4/1")
    plt.plot(time*1000,np.abs(psi32),label="3/2")
    plt.plot(time*1000,np.abs(psi42),label="4/2")
    #plt.plot(time*1000,np.abs(psi52),label="5/2")
    plt.plot(time*1000,np.abs(psi43),label="4/3")
    #plt.plot(time*1000,np.abs(psi54),label="5/4")
    #plt.plot(time*1000,np.abs(psi65),label="6/5")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=20)
    #plt.title(r"$B_r$",fontsize=16)
    #plt.title(r"Vacuum $B_r$",fontsize=16)
    plt.title(title,fontsize=24)
    plt.ylabel(ytitle,fontsize=24)
    plt.xlabel(r't [ms]',fontsize=24)
    plt.tight_layout()
    #ax.fill([0,1.0,1.0,0],[0,0,100,100],'gray',alpha=0.2)
    ax.fill([0,1.0,1.0,0],[0,0,40,40],'gray',alpha=0.2)
    plt.show()

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,phase21,ls=':',label="2/1",alpha=0.4)
    plt.plot(time*1000,phase31,label="3/1")
    plt.plot(time*1000,phase41,label="4/1")
    ax.legend(loc='upper left',frameon=True,ncol=3,fontsize=20)
    plt.title(r"Evolution of $\psi$ Phase",fontsize=24)
    plt.ylabel(r'$Phase$',fontsize=24)
    plt.xlabel(r't [ms]',fontsize=24)
    ax.fill([0,1.0,1.0,0],[-np.pi,-np.pi,np.pi,np.pi],'gray',alpha=0.2)
    plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],
               [r"$-\pi$",r"$-\pi/2$",0,r"$\pi/2$",r"$\pi$"])
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,phase21,label="2/1")
    plt.plot(time*1000,phase31,label="3/1")
#    plt.plot(time*1000,phase41,label="4/1")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=14)
    plt.title(r"Phase",fontsize=16)
    plt.ylabel(r'$Phase$',fontsize=16)
    plt.xlabel(r't [ms]',fontsize=16)
    plt.tight_layout()
    plt.show()


    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,np.abs(psi21),label="2/1")
    plt.plot(time*1000,np.abs(psi31),label="3/1")
    plt.plot(time*1000,np.abs(psi41),label="4/1")
    plt.plot(time*1000,np.abs(psi42),label="4/2")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=14)
    maxt=[2.18, 2.60, 2.94, 3.30, 3.66, 4.01, 4.38, 4.72, 5.08, 5.41]
    for t in maxt:
      plt.axvline(x=t,color='k',ls=':')
    #plt.title(r"$B_r$",fontsize=16)
    #plt.title(r"Vacuum $B_r$",fontsize=16)
    plt.title(title,fontsize=16)
    plt.ylabel(ytitle,fontsize=16)
    plt.xlabel(r't [ms]',fontsize=16)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(8,6))
    rel31=phase31 - phase21
    rel31 = np.where(rel31 < np.pi, rel31, rel31 - np.pi)
    rel31 = np.where(rel31 > -np.pi, rel31, rel31 + np.pi)
    rel41=phase41 - phase21
    rel41 = np.where(rel41 < np.pi, rel41, rel41 - np.pi)
    rel41 = np.where(rel41 > -np.pi, rel41, rel41 + np.pi)

    ax=fig.add_subplot(111)
    #plt.plot(time*1000,phase21,label="2/1")
    plt.plot(time*1000,rel31,label="3/1")
    plt.plot(time*1000,rel41,label="4/1")
    for t in maxt:
      plt.axvline(x=t,color='k',ls=':')
    plt.axhline(y=0,color='k')
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=14)
    plt.title(r"Phase relative to 2/1",fontsize=16)
    plt.ylabel(r'$Phase$',fontsize=16)
    plt.xlabel(r't [ms]',fontsize=16)
    plt.tight_layout()
    plt.show()

    fig = plt.figure(figsize=(8,6))
    rel34= phase31 - phase41
    rel34 = np.where(rel34 < np.pi, rel34, rel34 - np.pi)
    rel34 = np.where(rel34 > -np.pi, rel34, rel34 + np.pi)

    ax=fig.add_subplot(111)
    #plt.plot(time*1000,phase21,label="2/1")
    plt.plot(time*1000,rel34,label="3/1")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=14)
    plt.title(r"Phase relative to 4/1",fontsize=16)
    plt.ylabel(r'$Phase$',fontsize=16)
    plt.xlabel(r't [ms]',fontsize=16)
    plt.tight_layout()
    plt.show()

def surfmn_runner(show_plot=True,pickle_data=False,\
                  read_pickle=False,args=[]):
  ''' Perform surfmn analysie based on the the options
      If read pickle then it will look for pickle files else it will
      look for dump files
  '''
  dump_pre=["dumpgll","dump"]
  dump_suf=["h5"]
  pickle_pre=["surfmn"]
  pickle_suf=["pickle"]
  nimrodin="nimrod.in"
  steplist=[]
  read_new = True
  if read_pickle:
    pickle_list=glob.glob("surfmn*.pickle")
    pickle_list.sort(key=pickle_sort)
    print(pickle_list)
    if len(pickle_list)<=0:
      print("No pickle files found")
      raise IOError
    for iobj in pickle_list:
      with open(iobj,'rb') as file:
        surf=surfmn.fsasurfmn(None,None)
        surf.load(file)
        steplist.append(surf)
  #        steplist.append(pickle.load(open(iobj, "rb" )))
  else: #read from dump files and calculate fsa
    dumplist=[]
    if args['folder']:
      workdir=os.getcwd()
      listobjs = os.listdir(workdir)
      listobjs.sort()
      for iobj in listobjs:
        if os.path.isdir(iobj):
          thisdir=workdir+'/'+iobj+'/dumpgll*.h5'
          dumpfiles=glob.glob(thisdir)
          for file in dumpfiles:
            dumplist.append(file)
    else:
      dumplist=glob.glob("dumpg*.h5")
    dumplist.sort(key=dump_sort)
    if not os.path.isfile(nimrodin):
      print(f"nimrod.in not found")
      raise IOError

    rzo=np.array(args.get("rzo",[1.768,-0.018831,0.0]))
    nsurf=args.get("nsurf",150)
    fargs={}
    eqflag=args.get("eqflag",1)
    mmax=args.get("mmax",10)
    nmax=args.get("nmax",5)
    nmodes=args.get('nmodes',-1)
    if nmodes <1 or nmodes>nmax :
      fargs['ifour']=list(range(1,nmax+1))
    else:
      start=nmax-nmodes+1
      fargs['ifour']=list(range(start,nmax+1))
    fargs['mmax']=mmax

    for file_name in dumplist:
      print(file_name)
      print(nimrodin)
      surf=surfmn.fsasurfmn(file_name,nimrodin)
      surf.get_dumptime()
      surf.calculate(rzo=rzo,nsurf=nsurf,eqflag=eqflag,fargs=fargs)
      steplist.append(surf)
      if pickle_data:
        pfile=pickle_pre[0]+'.'+str(surf.step).zfill(5)+'.'+pickle_suf[0]
        print(f"writing file {pfile}")
        with open(pfile,'wb') as file:
          surf.dump(file)
  if show_plot:
    time_hist(steplist,pargs=args)
    plotlist=[100000,120000, 200000,300000,400000,steplist[-1]]
    args = {}
    #args['mlists'] = [[0,1,2,3,4,5,6]]
    #args['qlists'] = [[2,3,4]]
    for step in steplist:
        if step.step in plotlist:
            step.plot(pargs=args)
    step=steplist[-1]
    step.plot(pargs=args)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Surfmn runner.')
  parser.add_argument('--plot', action='store_true',help='shows plots')
  parser.add_argument('--folder', action='store_true',help='indicates the the dump files are stored in folders')
  parser.add_argument('--pickle', action='store_true',help='pickle data')
  parser.add_argument('--read', '-r', action='store_true',help='read pickled data')
  parser.add_argument('--mmax',type=int,default=15,help="max poloidal mode number")
  parser.add_argument('--nmax',type=int,default=5, help="max toroidal mode number")
  parser.add_argument('--nmodes', type=int, default=-1, \
    help="number of toroidal modes, defualt goes from 1 to nmax")
  parser.add_argument('--rzo',type=float, nargs=3, default=[1.768,-0.018831,0.0], help="intial guess for o-point")
  parser.add_argument('--nsurf', type=int, default=150, help="number of surfaces")
  parser.add_argument('--eqflag', type=int, default=1, help="flag to add n=0 perturbation to eq")
  parser.add_argument('--field', default='p', help="plot psi[p] or br[b]")
  args = vars(parser.parse_args())
  print(args)
  surfmn_runner(show_plot=args['plot'],\
                pickle_data=args['pickle'],read_pickle=args['read'],args=args)
