#!/usr/bin/env python3
#
#
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
from shutil import copy2
import hcStepFsa as step
import nim_timer as nimtime
import matplotlib.colors as mcolors


def hcrunner(file_name=None,
             pickle_data=False,
             read_pickle=False,
             args={}):

    if not os.path.isfile(file_name):
        print(f"File {file_name} not found")
        raise IOError

    dump_pre=["dumpgll","dump"]
    dump_suf=["h5"]
    pickle_suf=["pickle"]
    pickle_pre=["fsapower"]
    nimrodin="nimrod.in"
    pre=file_name.split('.')[0]
    if pre in dump_pre:
      print(f"Performing hc analysis from dump file")
      # check for nimrod.in and hdf5 format
      if not os.path.isfile(nimrodin):
        print(f"nimrod.in not found")
        raise IOError
      if not file_name.split('.')[-1] in dump_suf:
        print(f"dump file is not hdf5 format")
        raise IOError

      hc=step.hcstepfsa(file_name,nimrodin,args['lphi'])
      hc.get_dumptime()
      hc.calculate_power_fsa(nsurf=args['npts'],dpow=args['dpow'])
#      hc.analyze_power(npts=args['npts'],plot=True)
#      hc.analyze_power_adv(npts=args['npts'],plot=True)
#      hc.print_integrals()
      nimtime.timer.print_times()
      print(hc.step, hc.time)

      #pickle data here
      if args['pickle']:
        pfile=pickle_pre[0]+'.'+str(hc.step).zfill(5)+'.'+pickle_suf[0]
        print(f"writing file {pfile}")
        with open(pfile,'wb') as file:
          hc.dump(file)
    elif pre in pickle_pre:
      print("pickle_pre")
      hc=step.hcstepfsa(None,None,None)
      hc.load(file_name)
      print(f"Time: {hc.time}" )
     # hc.print_integrals()
    else:
      print(f"File {file_name} is not a recognized file type")
      raise IOError

    #plot data here
    if args['plot']:
        hc.interpolate_fsa(radial='rhon',npts=200)
        hc.default_plot()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Ho-Cradick FSA runner.')
    parser.add_argument('file',help='file name')
    parser.add_argument('--plot', action='store_true',help='shows plots')
    parser.add_argument('--pickle', action='store_true',help='pickle data')
    parser.add_argument('--npts', '-n', type=int, default=100,help='number of surfaces')
    parser.add_argument('--lphi', '-l', type=int, default=5, help='lphi')
    parser.add_argument('--dpow', '-d', type=float, default=0.5, help='dpow')
    parser.add_argument('--read', '-r', action='store_true',help='read pickled data')
    args = vars(parser.parse_args())
    print(args)
    hcrunner(file_name=args['file'],args=args)#\
                #pickle_data=args['pickle'],read_pickle=args['read'],args=args)
