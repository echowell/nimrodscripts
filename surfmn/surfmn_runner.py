#!/usr/bin/env python3
import os
import h5py
import surfmnstep
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob

''' This is a generic runner for surfmn. It loops over a bunch of directorys
in the current directory, searchs for a surfmn file, and a dump file. If the
dump file exists, then record the time. If not only record step number. It
then calls surfmn routines to plot the data and record the reconnected flux'''

def pickle_sort(file):
  print(file[6:])
  return int(file[6:])

def get_dumptime(thisfile):
  ''' Open an hdf5 file and read the dump time
      if I can't open file return None
  '''
  time=None
  with h5py.File(thisfile, 'r') as h5file:
    try:
      time=h5file["dumpTime"].attrs['vsTime']
    except:
      print(f"No dumpTime in dumpfile {thisfile}")
      raise
  return time

def find_files(thisdir):
  ''' This fuction finds the surfmn, dump files in the directory
      input: thisdir
      output: surfmn filename, dumpfilename, stepnumber, step time
      the outputs return None, if a file does not exists
  '''
  surfmn_file=None
  dumpfile=None
  stepnumber=None
  steptime=None
  listobjs = os.listdir(thisdir)
  for iobj in listobjs:
    wordlist = iobj.split('.')
    if (wordlist[0].lower()=='dumpgll'):
      if (dumpfile==None):
        dumpfile=thisdir+'/'+iobj
        thisstep=int(wordlist[1])
        if (stepnumber==None):
          stepnumber=thisstep
        elif (stepnumber!=thisstep):
          print(f"Dump step does not match surfmn step")
          raise
        steptime=get_dumptime(dumpfile)
      else:
        print(f"Multiple dumpfiles in directory {thisdir}")
        raise
    elif (wordlist[0].lower()=='surfmn'):
      if (surfmn_file==None):
        surfmn_file=thisdir+'/'+iobj
        thisstep=int(wordlist[1])
        if (stepnumber==None):
          stepnumber=thisstep
        elif (stepnumber!=thisstep):
          print(f"Surfmn step does not match dump step")
          raise
      else:
        print(f"Multiple surfmn files in directory {thisdir}")
        raise

  return surfmn_file, dumpfile, stepnumber, steptime

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
  mlist=[-1,-2,-3,-4]
  qlist=[-1,-2,-3,-4]
  for istep,step in enumerate(steplist):
    print(istep,step.step, step.time)
    time[istep]=step.time
    if step.surfmn_data==False:
      step.read_surfmn()

    psi21[istep]=step.get_resonance("psi",1,-2)
    psi31[istep]=step.get_resonance("psi",1,-3)
    psi41[istep]=step.get_resonance("psi",1,-4)
    psi43[istep]=step.get_resonance("psi",3,-4)
    psi32[istep]=step.get_resonance("psi",2,-3)
    psi54[istep]=step.get_resonance("psi",4,-5)
    psi65[istep]=step.get_resonance("psi",5,-6)
    if step.step in []:#[10000,18680,28000,66000,96000]:
#    if step.time>0: #avoids issues if no pertubation
      step.plot_surfmn("psi",1,**{"scale":1000})
      step.plot_surfmn("psi",2,**{"scale":1000})
      step.plot_surfmn("psi",3,**{"scale":1000})
      step.plot_surfmn("psi",4,**{"scale":1000})
      step.plot_surfmn("psi",5,**{"scale":1000})
      step.plot_radial("psi",1,mlist,**{"scale":1,"qlist":qlist})
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi21*1000)
  plt.title(r"$\psi$ 2/1",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi31*1000)
  plt.title(r"$\psi$ 3/1",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi41*1000)
  plt.title(r"$\psi$ 4/1",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi43*1000)
  plt.title(r"$\psi$ 4/3",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi32*1000)
  plt.title(r"$\psi$ 3/2",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi54*1000)
  plt.title(r"$\psi$ 5/4",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi65*1000)
  plt.title(r"$\psi$ 6/5",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,psi21*1000,label="2/1")
  plt.plot(time*1000,psi31*1000,label="3/1")
  plt.plot(time*1000,psi32*1000,label="3/2")
  plt.plot(time*1000,psi43*1000,label="4/3")
  plt.plot(time*1000,psi54*1000,label="5/4")
  plt.plot(time*1000,psi65*1000,label="6/5")
  ax.legend(loc=0)
  plt.title(r"$\psi$",fontsize=16)
  plt.ylabel(r'$\psi$ [mWb] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.show()

#  for istep in steplist:

def surfmn_runner(show_plot=True,pickle_data=False,read_pickle=False):
  ''' main runner for surfmn
      loops over all objects in a directory
      checks to see if the objects are a directory
      if so, then searches that directoy for a dump file and surfmn file'''
  steplist=[]
  read_new = True
  if read_pickle:
    pickle_list=glob.glob("pickle*")
    pickle_list.sort(key=pickle_sort)
    if len(pickle_list)>0:
      read_new=False
      for iobj in pickle_list:
        with open(iobj,'rb') as file:
          step=surfmnstep.SurfmnStep(None, None, None, None)
          step.load(file)
          steplist.append(step)
#        steplist.append(pickle.load(open(iobj, "rb" )))
  if read_new==True:
    workdir=os.getcwd()
    listobjs = os.listdir(workdir)
    listobjs.sort()
    for iobj in listobjs:
      if os.path.isdir(iobj):
        thisdir=workdir+'/'+iobj
        surfmn_file, dump, step, time = find_files(thisdir)
        steplist.append(surfmnstep.SurfmnStep(surfmn_file, dump, step, time))
  if show_plot:
    time_hist(steplist)
  if pickle_data:
    for step in steplist:
      if step==None:
        continue
      if step.surfmn_data==False:
        step.read_surfmn()
      filename="pickle"+str(step.step).zfill(5)
      with open(filename,'wb') as file:
        step.dump(file)
#      pickle.dump(step,open(filename,'wb'))

#  steplist[1].read_surfmn()
#  print(steplist[1].get_resonance("psi",1,-2))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Surfmn runner.')
  parser.add_argument('--plot', action='store_true',help='shows plots')
  parser.add_argument('--pickle', action='store_true',help='pickle data')
  parser.add_argument('--read', '-r', action='store_true',help='read pickled data')
  args = vars(parser.parse_args())
  surfmn_runner(show_plot=args['plot'],pickle_data=args['pickle'],read_pickle=args['read'])
