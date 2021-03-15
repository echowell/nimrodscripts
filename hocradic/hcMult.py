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
import hcStep as step
import nim_timer as nimtime
import matplotlib.colors as mcolors


def pickle_sort(file):
  print(file[5:])
  return int(file.split[0][5:])

def find_files(thisdir):
  ''' This fuction finds the hc, dump files in the directory
      input: thisdir
      output: hc filename, dumpfilename, stepnumber, step time
      the outputs return None, if a file does not exists
  '''
  dumpfile=None
  nimrodin=None
  listobjs = os.listdir(thisdir)
  for iobj in listobjs:
    wordlist = iobj.split('.')
    if (wordlist[0].lower()=='dumpgll' and wordlist[-1]=='h5'):
      if (dumpfile==None):
        dumpfile=iobj
      else:
        print(f"Multiple dumpfiles in directory {thisdir}")
        raise
    elif (iobj=='nimrod.in'):
      nimrodin=iobj

  return dumpfile, nimrodin

def hcmult(args):
    dump_pre=["dumpgll","dump"]
    dump_suf=["h5"]
    pickle_suf=["pickle"]
    pickle_pre=["power"]
    steplist = []
    read_new = True
    if args['read']:
        pickle_list=glob.glob("power*")
        pickle_list.sort(key=pickle_sort)
        if len(pickle_list)>0:
            read_new=False
            for iobj in pickle_list:
                with open(iobj,'rb') as file:
                    this=step.hcstep(None,None)
                    this.load(file)
                    steplist.append(this)
    if read_new==True:
        workdir=os.getcwd()
        listobjs = os.listdir(workdir)
        listobjs.sort()
        for iobj in listobjs:
            if os.path.isdir(iobj):
                thisdir=workdir+'/'+iobj
                dump, nimrodin = find_files(thisdir)
                if dump == None:
                    continue
                try:
                    os.mkdir('tempdir')
                except:
                    print("tempdir exists")
                copy2(thisdir + '/' + dump,'./tempdir')
                copy2('nimrod.in','./tempdir')
                os.chdir('tempdir')
                this=step.hcstep(dump,nimrodin)
                this.get_dumptime()
                this.analyze_power(npts=args['npts'])
                for iobj in os.listdir('.'):
                    os.remove(iobj)
                os.chdir('../')
                os.rmdir('tempdir')
                this.print_integrals()
                steplist.append(this)
                if args['pickle']:
                    pfile=pickle_pre[0]+'.'+str(this.step).zfill(5)+ \
                        '.'+pickle_suf[0]
                    print(f"writing file {pfile}")
                    with open(pfile,'wb') as file:
                        this.dump(file)

    for this in steplist:
        this.print_integrals()
        nimtime.timer.print_times()
        print(hc.step, hc.time)
    #plot data here
    if args['plot']:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Ho-Cradick runner.')
    parser.add_argument('--plot', action='store_true',help='shows plots')
    parser.add_argument('--pickle', action='store_true',help='pickle data')
    parser.add_argument('--read', '-r', action='store_true',help='read pickled data')
    parser.add_argument('--npts', '-n', type=int, default=512,help='number of points in 1D')
    args = vars(parser.parse_args())
    print(args)
    hcmult(args=args)
