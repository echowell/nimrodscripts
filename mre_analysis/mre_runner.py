#!/usr/bin/env python3
"""Modified Rutherford Equation Runner

This script allows the user to perform a Modified Rutherford Equation
type analysis of the supplied NIMROD dump file or saved pickle file.

The script handles the IO file operations for a single time step. The
actual analysis is performed in the module mre_step.py.

The input file is either a nimrod h5 dumpfile or are mre.pickle file.
If a dumpfile is supplied, then the correspond "nimrod.in" also needs
to be present.

If the pickle command line argument is specified, then the calculated
data will be saved in a pickle for further reuse.
"""


import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
from shutil import copy2
import mre_step
import nim_timer as nimtime
import matplotlib.colors as mcolors


def mre_main(file_name=None, args={}):
    """ Gets the file and calls the MRE analysis if not pickled

        Parameters
        ----------
        file_name : str
            The file to be read. It can either be a nimrod dumpfile or 
            a pickle file
        
        Returns
        ----------
          None
    """

    if not os.path.isfile(file_name):
        print(f"File {file_name} not found")
        raise IOError
    dump_pre = ["dumpgll", "dump"]
    dump_suf = ["h5"]
    pickle_suf = ["pickle"]
    pickle_pre = ["mre"]
    nimrodin = "nimrod.in"
    pre = file_name.split('.')[0]
    if pre in dump_pre:
        print(f"Performing hc analysis from dump file")
        # check for nimrod.in and hdf5 format
        if not os.path.isfile(nimrodin):
            print(f"nimrod.in not found")
            raise IOError
        if not file_name.split('.')[-1] in dump_suf:
            print(f"dump file is not hdf5 format")
            raise IOError
        mre = mre_step.MreStep(file_name, nimrodin)        
        mre.read_dumptime()
        mre.mre_analysis()
        #mre.calculate_power_fsa(nsurf=args['npts'],**args)
        nimtime.timer.print_times()
        time, step = mre.get_time()
        print(time, step)
        #pickle data here
        if args['pickle']:
            pfile = pickle_pre[0] + '.' + str(step).zfill(5) \
                    + '.' + pickle_suf[0]
            print(f"Writing file {pfile}")
            with open(pfile,'wb') as file:
                mre.dump(file)
    elif pre in pickle_pre:
      print("pickle_pre")
      mre = mre_step.MreStep(None, None)
      mre.load(file_name)
      print(f"Time: {mre.time}" )
    else:
      print(f"File {file_name} is not a recognized file type")
      raise IOError
    #plot data here
    if args['plot']:
        mre.interpolate_fsa(radial='rhon', npts=200, fsa=False)
        mre.default_plot()
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'MRE analysis')
    parser.add_argument('file', help = 'file name')
    parser.add_argument('--plot', action = 'store_true',
                        help = 'shows plots')
    parser.add_argument('--pickle', action = 'store_true',
                        help = 'pickle data')
    parser.add_argument('--npts', '-n', type = int, 
                        default = 100, help = 'number of surfaces')
    parser.add_argument('--dpow', '-d', type = float, default = 0.5,
                        help = 'Controls radial distribtuion of surfaces')
    args = vars(parser.parse_args())
    print(args)
    mre_main(file_name = args['file'], args = args)#\
                #pickle_data=args['pickle'],read_pickle=args['read'],args=args)
