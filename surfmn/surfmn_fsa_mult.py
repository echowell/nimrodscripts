#!/usr/bin/env python3
import os
import h5py
import surfmnstep
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

def time_hist(steplist):
    print(len(steplist))
    time=np.zeros(len(steplist))
    psi21=np.zeros(len(steplist))
    psi31=np.zeros(len(steplist))
    psi41=np.zeros(len(steplist))
    psi32=np.zeros(len(steplist))
    psi43=np.zeros(len(steplist))
    psi54=np.zeros(len(steplist))
    psi65=np.zeros(len(steplist))
    times=[]
    mlist=[-1,-2,-3,-4]
    qlist=[-1,-2,-3,-4]

    for istep,step in enumerate(steplist):
        print(istep,step.step, step.time)
        time[istep]=step.time
        #if step.surfmn_data==False:
        #    step.read_surfmn()
        psi21[istep]=step.get_resonance("psi",1,-2)
        psi31[istep]=step.get_resonance("psi",1,-3)
        psi41[istep]=step.get_resonance("psi",1,-4)
        psi43[istep]=step.get_resonance("psi",3,-4)
        psi32[istep]=step.get_resonance("psi",2,-3)
        psi54[istep]=step.get_resonance("psi",4,-5)
        psi65[istep]=step.get_resonance("psi",5,-6)
        if step.step==00000:
            #print(step.mr.shape,step.q.shape)
            eq_q2 = step.get_rho_q(q=-2)
            eq_q3 = step.get_rho_q(q=-3)
            eq_q65 = step.get_rho_q(q=-1.2)
            eq_q54 = step.get_rho_q(q=-1.25)
            eq_q43 = step.get_rho_q(q=-4./3.)
            eq_q32 = step.get_rho_q(q=-3./2.)

    fig = plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    plt.plot(time*1000,psi21*1000,label="2/1")
    plt.plot(time*1000,psi31*1000,label="3/1")
    plt.plot(time*1000,psi32*1000,label="3/2")
    plt.plot(time*1000,psi43*1000,label="4/3")
    plt.plot(time*1000,psi54*1000,label="5/4")
    plt.plot(time*1000,psi65*1000,label="6/5")
    plt.plot(time*1000,psi41*1000,label="4/1")
    ax.legend(loc='best',frameon=True,ncol=2,fontsize=14)
    plt.title(r"$\psi$",fontsize=16)
    plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
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
      print("After init")
      surf.get_dumptime()
      print("after dumptime")
      surf.calculate(rzo=rzo,nsurf=nsurf,eqflag=eqflag,fargs=fargs)
      print("after calculate")
      steplist.append(surf)
      if pickle_data:
        pfile=pickle_pre[0]+'.'+str(surf.step).zfill(5)+'.'+pickle_suf[0]
        print(f"writing file {pfile}")
        with open(pfile,'wb') as file:
          surf.dump(file)
  if show_plot:
    time_hist(steplist)
    plotlist=[54000]
    for step in steplist:
        if step.step in plotlist:
            step.plot()
    pass




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
  args = vars(parser.parse_args())
  print(args)
  surfmn_runner(show_plot=args['plot'],\
                pickle_data=args['pickle'],read_pickle=args['read'],args=args)
