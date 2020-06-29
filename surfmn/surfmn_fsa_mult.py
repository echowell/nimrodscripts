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
#    time_hist(steplist)
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
