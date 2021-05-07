#!/usr/bin/env python3
import os
import h5py
#import surfmnstep
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import glob
import fsa
from shutil import copy2
#import nim_timer as timer
import f90nml
import eval_nimrod as eval
from scipy.interpolate import interp1d,splev,UnivariateSpline,griddata
import labellines as ll
import matplotlib.ticker as ticker


import matplotlib.colors as mcolors

#uses labelines package: https://github.com/cphyc/matplotlib-label-lines/
#pip3 install matplotlib-label-lines

def dummy_fsa(rzc,y,dy,evalnimrod,fargs):
  '''
  Dummy integrand for fsa, this is used to get and q
  without running a true fsa
  Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
  dy(0)=dl/deta or d eta/dl
  dy(1)=dr/deta or dr/dl
  dy(2)=1/bdgth :v'
  dy(3)=dq: q
  '''
  dy[4]=1.0
  return dy

#@timer.timer_func
def get_qsurfaces(dumpfile,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
    qlist=[-1.25, -1.50, -2.0, -3.0,-4.0]
    qlist2=[-1.125,-1.75,-2.5,-3.5]
    evalNimrod=eval.EvalNimrod(dumpfile,fieldlist='b')
#
    dvar, yvar, contours = fsa.FSA(evalNimrod, rzo, dummy_fsa, 1, \
        nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
        fargs=fargs)

    iend=-1
    while np.isnan(yvar[:,iend]).any():
        iend -= 1
    iend += yvar.shape[1]+1
    #unevaluated interpoate
    qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic',fill_value="extrapolate")
    Rq=interp1d(dvar[7,:iend], dvar[4,:iend], kind='cubic',fill_value="extrapolate")
    Zq=interp1d(dvar[7,:iend], dvar[5,:iend], kind='cubic',fill_value="extrapolate")

    contourdict={}
    for iq in qlist:
        iyvar, icontour = fsa.FSA(evalNimrod,rzo, dummy_fsa, 1, \
            nsurf=1,depvar='eta',dpow=0.5,rzp=[Rq(iq),Zq(iq)],flag=eqflag,normalize=True, \
            fargs=fargs)
        contourdict[iq]=icontour

    contourdict2={}
    for iq in qlist2:
        iyvar, icontour = fsa.FSA(evalNimrod,rzo, dummy_fsa, 1, \
            nsurf=1,depvar='eta',dpow=0.5,rzp=[Rq(iq),Zq(iq)],flag=eqflag,normalize=True, \
            fargs=fargs)
        contourdict2[iq]=icontour

    return contourdict, contourdict2

def q_runner(fileName,plot,pickle,args):
    dump_pre=["dumpgll","dump"]
    dump_suf=["h5"]
    pickle_pre=["ntm"]
    pickle_suf=["pickle"]
    nimrodin="nimrod.in"
    pre=fileName.split('.')[0]
    step=fileName.split('.')[1]
    if pre in dump_pre:
        print(f"Performing ntm analysis from dump file")
        # check for nimrod.in and hdf5 format
        if not os.path.isfile(nimrodin):
            print(f"nimrod.in not found")
            raise IOError
        if not fileName.split('.')[-1] in dump_suf:
            print(f"dump file is not hdf5 format")
            raise IOError
        nsurf=args.get("nsurf",150)
        fargs={}
        eqflag=args.get("eqflag",1)
        print(fileName)
        qsurfaces, qsurfaces2=get_qsurfaces(fileName,rzo=np.array(args['rzo']),nsurf=nsurf,fargs=fargs)
    #  timer.timer.print_times()

    #pickle data here
        if args['pickle']:
            pfile=pickle_pre[0]+'.'+str(step).zfill(5)+'.'+pickle_suf[0]
            print(f"writing file {pfile}")
            with open(pfile,'wb') as file:
                pass
    #            ntm.dump(file)
    elif pre in pickle_pre:
        pass
      #print("pickle_pre")
      #ntm=step.ntmstep(None,None)
      #ntm.load(file_name)
      #print(f"Time: {ntm.time}" )
      #ntm.plot_fsa_phase(key=None)
#    with open(file_name,'rb') as file:
#      surf=surfmn.fsasurfmn(None,None)
#      surf.load(file)
    else:
        #print(f"File {file_name} is not a recognized file type")
        raise IOError

#plot data here
    if plot:
        xvals=[]
        llTune=8
        fig = plt.figure(figsize=(6,5))
        ax=fig.add_subplot(111)
        plt.title("Equilibrium q-surfaces",size=16)
        for ii, [iq, isurface] in enumerate(qsurfaces.items()):
            plt.plot(isurface[0,:],isurface[1,:],label=f"{iq}",ls='-')
            xvals.append(isurface[0,int((3*ii)%llTune*isurface.shape[1]/llTune)])
        for ii, [iq, isurface] in enumerate(qsurfaces2.items()):
            plt.plot(isurface[0,:],isurface[1,:],ls='-',color='k', linewidth=1)
        plt.scatter( 1.76893653, -0.01890057, marker='.',color='k',s=4)
        ll.labelLines(plt.gca().get_lines(),xvals=xvals,fontsize=10)
        ax.set_xlabel(r"R [m]",size=16)
        ax.set_ylabel(r"Z [m]",size=16,rotation=90)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        #ax.set_xlim([1,2.5])
        plt.axis('equal')
        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
      parser = argparse.ArgumentParser(description='plots eq q surfaces runner.')
      parser.add_argument('file',help='file name')
      parser.add_argument('--plot', action='store_true',help='shows plots')
      parser.add_argument('--pickle', action='store_true',help='pickle data')
      parser.add_argument('--nsurf', type=int, default=150, help="number of surfaces")
      parser.add_argument('--eqflag', type=int, default=1, help="flag to add n=0 perturbation to eq")
      parser.add_argument('--rzo',type=float, nargs=3, default=[1.768,-0.018831,0.0], help="intial guess for o-point")
      args = vars(parser.parse_args())
      print(args)
      q_runner(fileName=args['file'],plot=args['plot'],\
                    pickle=args['pickle'],args=args)
