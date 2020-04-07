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
  nimrodin=None
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
    elif (iobj=='nimrod.in'):
      nimrodin=thisdir+'/'+iobj


  return surfmn_file, dumpfile, stepnumber, steptime, nimrodin

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
  exb21=np.zeros(len(steplist))
  exb31=np.zeros(len(steplist))
  exb32=np.zeros(len(steplist))
  exb43=np.zeros(len(steplist))
  mlist=[-1,-2,-3,-4]
  qlist=[-1,-2,-3,-4]
  for istep,step in enumerate(steplist):
    print(istep,step.step, step.time)
    time[istep]=step.time
    if step.surfmn_data==False:
      step.read_surfmn()
    if step.profdata==False:
      try:
        os.mkdir('tempprofile')
      except:
        print("tempprofile directoy exists")
      copy2(step.dumpfile,'./tempprofile')
      copy2(step.nimrodin,'./tempprofile')
      os.chdir('tempprofile')
      step.get_profiles()
      for iobj in os.listdir('.'):
        os.remove(iobj)
      os.chdir('../')
      os.rmdir('tempprofile')
#    fig = plt.figure(figsize=(6,5))
#    ax=fig.add_subplot(111)
#    plt.plot(step.profs.rhon,step.profs.omegator)
#    plt.title(r"fsa",fontsize=16)
#    plt.ylabel(r'nd ',fontsize=16)
#    plt.xlabel(r'rho',fontsize=16)
#    plt.tight_layout()
#    plt.show()
    psi21[istep]=step.get_resonance("psi",1,-2)
    psi31[istep]=step.get_resonance("psi",1,-3)
    psi41[istep]=step.get_resonance("psi",1,-4)
    psi43[istep]=step.get_resonance("psi",3,-4)
    psi32[istep]=step.get_resonance("psi",2,-3)
    psi54[istep]=step.get_resonance("psi",4,-5)
    psi65[istep]=step.get_resonance("psi",5,-6)
    this_q=step.profs.get_rho_q(q=-2)
    exb21[istep]=step.profs.get_omega_exb(n=1,rhon=this_q)/(2*np.pi)
    this_q=step.profs.get_rho_q(q=-3)
    exb31[istep]=step.profs.get_omega_exb(n=1,rhon=this_q)/(2*np.pi)
    this_q=step.profs.get_rho_q(q=-1.5)
    exb32[istep]=step.profs.get_omega_exb(n=2,rhon=this_q)/(2*np.pi)
    this_q=step.profs.get_rho_q(q=-4/3)
    exb43[istep]=step.profs.get_omega_exb(n=3,rhon=this_q)/(2*np.pi)
    if step.step in [2000000]:#[10000,18680,28000,66000,96000]:
#    if step.time>0: #avoids issues if no pertubation
      this_exb=step.profs.get_omega_exb(n=1)
      this_rho_21=step.profs.get_rho_q(q=-2)
      fig = plt.figure(figsize=(6,5))

      fig = plt.figure(figsize=(6,5))
      ax=fig.add_subplot(111)
      plt.plot(step.profs.rhon,this_exb/(2*np.pi))
      plt.title(r"fsa",fontsize=16)
      plt.ylabel(r'f_exb ',fontsize=16)
      plt.xlabel(r'rho',fontsize=16)
      ax.axvline(this_rho_21,ls=':')
      plt.tight_layout()
      plt.show()
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
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,np.abs(exb21)/1000)
  plt.title(r"$f_{eb}$ 2/1",fontsize=16)
  plt.ylabel(r'f [kHz] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.tight_layout()
  plt.show()
  fig = plt.figure(figsize=(6,5))
  ax=fig.add_subplot(111)
  plt.plot(time*1000,np.abs(exb43)/1000)
  plt.plot(time*1000,np.abs(exb32)/1000)
  plt.plot(time*1000,np.abs(exb31)/1000)
  plt.plot(time*1000,np.abs(exb21)/1000)
  plt.title(r"$f_{eb}$ ",fontsize=16)
  plt.ylabel(r'f [kHz] ',fontsize=16)
  plt.xlabel(r't [ms]',fontsize=16)
  plt.tight_layout()
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
          step=surfmnstep.SurfmnStep(None, None, None, None, None)
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
        surfmn_file, dump, step, time, nimrodin = find_files(thisdir)
        steplist.append(surfmnstep.SurfmnStep(surfmn_file, dump, step, time,nimrodin))
  if show_plot:
    time_hist(steplist)
  if pickle_data:
    for step in steplist:
      if step==None:
        continue
      if step.surfmn_data==False:
        step.read_surfmn()
      if step.profdata==False:
        try:
          os.mkdir('tempprofile')
        except:
          print("tempprofile directoy exists")
        copy2(step.dumpfile,'./tempprofile')
        copy2(step.nimrodin,'./tempprofile')
        os.chdir('tempprofile')
        step.get_profiles()
        for iobj in os.listdir('.'):
          os.remove(iobj)
        os.chdir('../')
        os.rmdir('tempprofile')
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
