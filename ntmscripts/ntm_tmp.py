#!/usr/bin/env python3
#
#profiles is a class for calculating 1D profiles
# using the flux surface integration
#
#
import os
import h5py
#import surfmnstep
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
from shutil import copy2
import ntm_step_real as step
import nim_timer as timer

#import fsa_surfmn as surfmn
import matplotlib.colors as mcolors




def ntm_tmp(file_name=None,show_plot=True,pickle_data=False,\
                  read_pickle=False,args=[]):
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
      ### todo
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
  #    timer=nim_timer.nimTimer()
      ntm=step.ntmstep(file_name,nimrodin)
      print("After ntmstep")
      ntm.get_dumptime()
      print("After dumptime")
      ntm.fields.set_method("induction")
      method='rad'
      if method=='rad':
          ntm.analyze_radial()
      else:
          ntm.analyze()


  #    ntm.calculate_psi(rzo=np.array(args['rzo']),nsurf=nsurf,fargs=fargs)
  #    raise
 #

#      ntm.calculate_induction(rzo=np.array(args['rzo']),nsurf=nsurf,fargs=fargs)
#      timer.timer.print_times()
#      ntm.plot_fsa_phase(key=None)
      #raise
      #ntm.analyze()

      #raise ValueError
      #ntm.get_domain()
      #ntm.eval_plot()
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
        pfile=pickle_pre[0]+'.'+str(ntm.step).zfill(5)+'.'+pickle_suf[0]
        print(f"writing file {pfile}")
        with open(pfile,'wb') as file:
          ntm.dump(file)
    elif pre in pickle_pre:
      print("pickle_pre")
      ntm=step.ntmstep(None,None)
      ntm.load(file_name)
      print(f"Time: {ntm.time}" )
      #ntm.plot_fsa_phase(key=None)
  #    with open(file_name,'rb') as file:
  #      surf=surfmn.fsasurfmn(None,None)
  #      surf.load(file)
    else:
      print(f"File {file_name} is not a recognized file type")
      raise IOError

    #plot data here
    if args['plot']:
      ntm.plot_fsa_phase(key=None,time=ntm.time)
  #    surf.plot()
    pass

















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
  ntm_tmp(file_name=args['file'],show_plot=args['plot'],\
                pickle_data=args['pickle'],read_pickle=args['read'],args=args)
