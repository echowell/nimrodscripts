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
import sys

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
    pickle_pre=["power","poweradv"]
    steplist = []
    read_new = True
    if args['merge']:
        pickle_list=glob.glob("power.*")
        if len(pickle_list)>0:
            for iobj in pickle_list:
                this=step.hcstep(None,None)
                this.load(iobj)
                advpickle=pickle_pre[1]+'.'+str(this.step).zfill(5)+ \
                    '.'+pickle_suf[0]
                print(this.step,advpickle)
                try:
                    this2=step.hcstep(None,None)
                    this2.load(advpickle)
                    for key, field in this2.powerDict.items():
                        print(key)
                        if key not in this.powerDict:
                           print(f"adding key {key} to power dict")
                           this.powerDict[key]=field
                    print("Found file")
                    if args['pickle']:
                        print(f"writing file {iobj}")
                        with open(iobj,'wb') as file:
                            this.dump(file)
                        print('pickle')
                except:
                    print("File not found")
        sys.exit(0)
    if args['read']:
        pickle_list=glob.glob("power.*")
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
                if args['mode']==0:
                    this.analyze_power(npts=args['npts'])
                    this.analyze_power_adv(npts=args['npts'])
                elif args['mode']==1:
                    this.analyze_power_adv(npts=args['npts'])
                else:
                    print(f"mode {args['mode']} is not valid")
                    raise ValueError
                for iobj in os.listdir('.'):
                    os.remove(iobj)
                os.chdir('../')
                os.rmdir('tempdir')
                this.print_integrals()
                steplist.append(this)
                if args['pickle']:
                    pfile=pickle_pre[args['mode']]+'.'+str(this.step).zfill(5)+ \
                        '.'+pickle_suf[0]
                    print(f"writing file {pfile}")
                    with open(pfile,'wb') as file:
                        this.dump(file)
                this.clean_up()
    for this in steplist:
        this.print_integrals()
        nimtime.timer.print_times()
        print(this.step, this.time)
    #plot data here
    if args['plot']:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Ho-Cradick runner.')
    parser.add_argument('--plot', action='store_true',help='shows plots')
    parser.add_argument('--pickle', action='store_true',help='pickle data')
    parser.add_argument('--read', '-r', action='store_true',help='read pickled data')
    parser.add_argument('--npts', '-n', type=int, default=512,help='number of points in 1D')
    parser.add_argument('--mode', '-m', type=int, default=0,help='0 standard, 1 advect ')
    parser.add_argument('--merge', action='store_true', help='copy advection to power dict ')
    args = vars(parser.parse_args())
    print(args)
    hcmult(args=args)
