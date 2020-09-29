#!/usr/bin/env python3
import os
import h5py
#import surfmnstep
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
from shutil import copy2
import ntm_step as step

#import fsa_surfmn as surfmn
import matplotlib.colors as mcolors


''' This is a generic script for analysing a single dumpfile/timestep
    using the surfmn fsa scripts '''

def ntm_runner(file_name=None,show_plot=True,pickle_data=False,\
                  read_pickle=False,args=[]):
  ''' Perform surfmn analysie based on the the options
      If file is a dump file, then it will perform the surfmn fsa
      If it is a surfmn.pickle file then it will read the fsa info from the pickle
      file.
  '''

  #check to see if file exists
  if not os.path.isfile(file_name):
    print(f"File {file_name} not found")
    raise IOError

  dump_pre=["dumpgll","dump"]
  dump_suf=["h5"]
  pickle_pre=["ntm"]
  pickle_suf=["pickle"]
  nimrodin="nimrod.in"
  pre=file_name.split('.')[0]
  if pre in dump_pre:
    print(f"Performing ntm analysis from dump file")
    # check for nimrod.in and hdf5 format
    if not os.path.isfile(nimrodin):
      print(f"nimrod.in not found")
      raise IOError
    if not file_name.split('.')[-1] in dump_suf:
      print(f"dump file is not hdf5 format")
      raise IOError
    ntm=step.ntmstep(file_name,nimrodin)
    ntm.get_domain()
    # initalize surfmn object
#    surf=surfmn.fsasurfmn(file_name,nimrodin)
#    surf.get_dumptime()
#    rzo=np.array(args.get("rzo",[1.768,-0.018831,0.0]))
#    nsurf=args.get("nsurf",150)
#    fargs={}
#    eqflag=args.get("eqflag",1)
#    mmax=args.get("mmax",10)
#    nmax=args.get("nmax",5)
#    nmodes=args.get('nmodes',-1)
#    if nmodes <1 or nmodes>nmax :
#      fargs['ifour']=list(range(1,nmax+1))
#    else:
#      start=nmax-nmodes+1
#      fargs['ifour']=list(range(start,nmax+1))
#    fargs['mmax']=mmax
#    surf.calculate(rzo=rzo,nsurf=nsurf,eqflag=eqflag,fargs=fargs)

    #pickle data here
    if args['pickle']:
      pfile=pickle_pre[0]+'.'+str(surf.step).zfill(5)+'.'+pickle_suf[0]
      print(f"writing file {pfile}")
      with open(pfile,'wb') as file:
        surf.dump(file)
  elif pre in pickle_pre:
    print("pickle_pre")
#    with open(file_name,'rb') as file:
#      surf=surfmn.fsasurfmn(None,None)
#      surf.load(file)
  else:
    print(f"File {file_name} is not a recognized file type")
    raise IOError

  #plot data here
  if args['plot']:
    pass
#    surf.plot()

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='NTM runner.')
  parser.add_argument('file',help='file name')
  parser.add_argument('--plot', action='store_true',help='shows plots')
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
  ntm_runner(file_name=args['file'],show_plot=args['plot'],\
                pickle_data=args['pickle'],read_pickle=args['read'],args=args)
