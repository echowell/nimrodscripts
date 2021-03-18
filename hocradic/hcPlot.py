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
  #print(file[5:],file.split('.'))
  return int(file.split('.')[1])

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

def hcplot(args):
    pickle_suf=["pickle"]
    pickle_pre=["power"]
    steplist = []
    nmode_list = [1,2,3,4,5]
    times = None
    if args['read']:
        pickle_list=glob.glob("power.*")
        pickle_list.sort(key=pickle_sort)
        time_steps = len(pickle_list)
        times = np.zeros(time_steps)
        nmodes = len(nmode_list)
        mag_energy = np.zeros([time_steps,nmodes])
        kin_energy = np.zeros([time_steps,nmodes])
        lin_power = np.zeros([time_steps,nmodes])
        qua_power = np.zeros([time_steps,nmodes])
        non_power = np.zeros([time_steps,nmodes])
        ohm_power = np.zeros([time_steps,nmodes])
        visc_power = np.zeros([time_steps,nmodes])
        neoi_power = np.zeros([time_steps,nmodes])
        neoe_power = np.zeros([time_steps,nmodes])
        poyn_power = np.zeros([time_steps,nmodes])
        press_power = np.zeros([time_steps,nmodes])
        total_power = np.zeros([time_steps,nmodes])
#        adv_power = np.zeros([time_steps,nmodes])
        if len(pickle_list)>0:
            for ii, file in enumerate(pickle_list):
                print(file)
                this=step.hcstep(None,None)
                this.load(file)
                steplist.append(this)
                times[ii] = this.time
                for jj, nn in enumerate(nmode_list):
                    mag_energy[ii,jj]=this.energyDict['bpert'][nn]
                    kin_energy[ii,jj]=this.energyDict['vpert'][nn]
                    lin_power[ii,jj]=this.powerDict['vxbeq'][nn] + \
                                     this.powerDict['jxbeq'][nn]
                    qua_power[ii,jj]=this.powerDict['vxbn0'][nn] + \
                                     this.powerDict['jxbn0'][nn]
                    non_power[ii,jj]=this.powerDict['vxbp'][nn] + \
                                     this.powerDict['jxbp'][nn]
                    ohm_power[ii,jj]=this.powerDict['etajp'][nn]
                    visc_power[ii,jj]=this.powerDict['divpip'][nn]
                    neoi_power[ii,jj]=this.powerDict['divPii'][nn]
                    neoe_power[ii,jj]=this.powerDict['divPie'][nn]
                    poyn_power[ii,jj]=this.powerDict['poynting'][nn]
                    press_power[ii,jj]=this.powerDict['ngpp'][nn]
                    #adv_power[ii,jj]=this.powerDict['ngpp'][nn]
                    for key in this.powerDict:
                        total_power[ii,jj]+=this.powerDict[key][nn]

    else:
        print("Hc plot can not read dumpfiles")
        raise KeyError
    figsize=[8,6]
    #plot magnetic energy
    fig, ax = plt.subplots(figsize=figsize)
    for jj, nn in enumerate(nmode_list):
        label = f"n = {nn}"
        plt.plot(times*1000,mag_energy[:,jj],label=label)
    plt.legend()
    plt.title("Magnetic Energy")
    plt.ylabel('Energy [J]')
    plt.xlabel('Time [ms]')
    plt.tight_layout()
    plt.show()

    #plot kinetic energy
    fig, ax = plt.subplots(figsize=figsize)
    for jj, nn in enumerate(nmode_list):
        label = f"n = {nn}"
        plt.plot(times*1000,kin_energy[:,jj],label=label)
    plt.legend()
    plt.title("Kinetic Energy")
    plt.ylabel('Energy [J]')
    plt.xlabel('Time [ms]')
    plt.tight_layout()
    plt.show()

    #plot power flux
    for jj, nn in enumerate(nmode_list):
        fig, ax = plt.subplots(figsize=figsize)
        Title = f"n = {nn} Power"
        plt.plot(times*1000,lin_power[:,jj]/1e3,label="Linear")
        plt.plot(times*1000,qua_power[:,jj]/1e3,label="Quasi-linear")
        plt.plot(times*1000,non_power[:,jj]/1e3,label="Nonlinear")
        plt.plot(times*1000,press_power[:,jj]/1e3,label="Grad P")
        plt.legend()
        plt.axhline(0,color='k')
        plt.title(Title)
        plt.ylabel('Power [kW]')
        plt.xlabel('Time [ms]')
        plt.tight_layout()
        plt.show()

    #n=1 power flux
    jj=0
    nn=1
    fig, ax = plt.subplots(figsize=figsize)
    Title = f"n = {nn} Power"
    plt.plot(times*1000,lin_power[:,jj]/1e3,label="L")
    plt.plot(times*1000,qua_power[:,jj]/1e3,label="QL")
    plt.plot(times*1000,non_power[:,jj]/1e3,label="NL")
    plt.plot(times*1000,press_power[:,jj]/1e3,label="GP")
    plt.plot(times*1000,ohm_power[:,jj]/1e3,label="etaJ")
#    plt.plot(times*1000,visc_power[:,jj]/1e3,label="Viscous")
#    plt.plot(times*1000,neoi_power[:,jj]/1e3,label="Neoclassical ion")
#    plt.plot(times*1000,neoe_power[:,jj]/1e3,label="Neoclassical electron")
#    plt.plot(times*1000,poyn_power[:,jj]/1e3,label="Poynting")
    plt.legend()
    plt.title(Title)
    plt.axhline(0,color='k')
    plt.ylim([-60,20])
    plt.ylabel('Power [kW]')
    plt.xlabel('Time [ms]')
    plt.tight_layout()
    plt.show()

    jj=1
    nn=2
    fig, ax = plt.subplots(figsize=figsize)
    Title = f"n = {nn} Power"
    plt.plot(times*1000,lin_power[:,jj]/1e3,label="L")
    plt.plot(times*1000,qua_power[:,jj]/1e3,label="QL")
    plt.plot(times*1000,non_power[:,jj]/1e3,label="NL")
    plt.plot(times*1000,press_power[:,jj]/1e3,label="GP")
    plt.plot(times*1000,ohm_power[:,jj]/1e3,label="etaJ")
    #    plt.plot(times*1000,visc_power[:,jj]/1e3,label="Viscous")
    #    plt.plot(times*1000,neoi_power[:,jj]/1e3,label="Neoclassical ion")
    #    plt.plot(times*1000,neoe_power[:,jj]/1e3,label="Neoclassical electron")
    #    plt.plot(times*1000,poyn_power[:,jj]/1e3,label="Poynting")
    plt.legend(loc=2)
    plt.title(Title)
    plt.axhline(0,color='k')
    plt.ylim([-200,250])
    plt.ylabel('Power [kW]')
    plt.xlabel('Time [ms]')
    plt.tight_layout()
    plt.show()

    jj=2
    nn=3
    fig, ax = plt.subplots(figsize=figsize)
    Title = f"n = {nn} Power"
    plt.plot(times*1000,lin_power[:,jj]/1e3,label="L")
    plt.plot(times*1000,qua_power[:,jj]/1e3,label="QL")
    plt.plot(times*1000,non_power[:,jj]/1e3,label="NL")
    plt.plot(times*1000,press_power[:,jj]/1e3,label="GP")
    plt.plot(times*1000,ohm_power[:,jj]/1e3,label="etaJ")
    #    plt.plot(times*1000,visc_power[:,jj]/1e3,label="Viscous")
    #    plt.plot(times*1000,neoi_power[:,jj]/1e3,label="Neoclassical ion")
    #    plt.plot(times*1000,neoe_power[:,jj]/1e3,label="Neoclassical electron")
    #    plt.plot(times*1000,poyn_power[:,jj]/1e3,label="Poynting")
    plt.legend(loc=2)
    plt.title(Title)
    plt.axhline(0,color='k')
    #plt.ylim([-200,250])
    plt.ylabel('Power [kW]')
    plt.xlabel('Time [ms]')
    plt.tight_layout()
    plt.show()

    for jj, nn in enumerate(nmode_list):
        fig, ax = plt.subplots(figsize=figsize)
        Title = f"n = {nn} Power"
#        plt.plot(times*1000,lin_power[:,jj]/1e3,label="Linear")
#        plt.plot(times*1000,qua_power[:,jj]/1e3,label="Quasi-linear")
#        plt.plot(times*1000,non_power[:,jj]/1e3,label="Nonlinear")
#        plt.plot(times*1000,ohm_power[:,jj]/1e3,label="Ohmic")
        plt.plot(times*1000,visc_power[:,jj]/1e3,label="Visc")
        plt.plot(times*1000,neoi_power[:,jj]/1e3,label="Neo ion")
        plt.plot(times*1000,neoe_power[:,jj]/1e3,label="Neo e-")
        plt.plot(times*1000,poyn_power[:,jj]/1e3,label="Poyn")
#        plt.plot(times*1000,press_power[:,jj]/1e3,label="Grad P")
        plt.legend()
        plt.title(Title)
        plt.axhline(0,color='k')
        plt.ylabel('Power [kW]')
        plt.xlabel('Time [ms]')
        plt.tight_layout()
        plt.show()

    for jj, nn in enumerate(nmode_list):
        fig, ax = plt.subplots(figsize=figsize)
        Title = f"n = {nn} Power"

#        plt.plot(times*1000,lin_power[:,jj]/1e3,label="Linear")
#        plt.plot(times*1000,qua_power[:,jj]/1e3,label="Quasi-linear")
#        plt.plot(times*1000,non_power[:,jj]/1e3,label="Nonlinear")
#        plt.plot(times*1000,ohm_power[:,jj]/1e3,label="Ohmic")
        plt.plot(times*1000,total_power[:,jj]/1e3,label="Total Power")
        plt.legend()
        plt.title(Title)
        plt.axhline(0,color='k')
        plt.ylabel('Power [kW]')
        plt.xlabel('Time [ms]')
        plt.tight_layout()
        plt.show()

    for this in steplist:
        pass
#        this.print_integrals()
#        nimtime.timer.print_times()
#        print(this.step, this.time)
    #plot data here
    if args['plot']:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Ho-Cradick runner.')
    parser.add_argument('--read', action='store_true',help='read pickled data')
    args = vars(parser.parse_args())
    print(args)
    hcplot(args=args)
