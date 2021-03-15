#!/usr/bin/env python3
import eval_nimrod as eval
import plot_nimrod as pn
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import h5py
import numpy as np
import pickle
import hcFields as hc
import nim_timer as timer

class hcstep:
    def __init__(self,dumpfile,nimrodin):
        self.dumpfile=dumpfile
        self.nimrodin=nimrodin
        self.fields=hc.hcfields(dumpfile,nimrodin)
        self.time=None
        self.step=None
        self.eval=None
        self.energyDict={}
        self.powerDict={}
        self.volume=None
        self.grid=None
        self.intSet=False
        self.intWeight=np.empty([1])

    def dump(self,file):
        pickle.dump(self.dumpfile,file)
        pickle.dump(self.nimrodin,file)
        pickle.dump(self.time,file)
        pickle.dump(self.step,file)
        pickle.dump(self.energyDict,file)
        pickle.dump(self.powerDict,file)
        pickle.dump(self.volume,file)

    def load(self,file):
        with open(file,'rb') as pickle_file:
            self.dumpfile=pickle.load(pickle_file)
            self.nimrodin=pickle.load(pickle_file)
            self.time=pickle.load(pickle_file)
            self.step=pickle.load(pickle_file)
            self.energyDict=pickle.load(pickle_file)
            self.powerDict=pickle.load(pickle_file)
            self.volume=pickle.load(pickle_file)

    def set_3dgrid(self,rmin,rmax,zmin,zmax,nr,nz,lphi,nonlin_order=2,debug=0):
        '''sets up a 3d grid, using non to determine the number of phi planes
           based on
        '''
        self.nmodes=self.calc_nmodes(lphi)
        self.nmax=self.nmodes-1
        nphi=self.calc_nplanes(lphi)
        phimax=np.pi*2*(nphi-1)/nphi
        p1 = np.array([rmin, zmin, 0.0])
        p2 = np.array([rmax, zmin, 0.0])
        p3 = np.array([rmin, zmax, 0.0])
        p4 = np.array([rmin, zmin, phimax])
        rzp3d = pn.PlotNimrod.grid_3d_gen(p1, p2, p3, p4, nr, nz,nphi)
        self.fields.grid=eval.EvalGrid(rzp3d)
        self.fields.grid.set_debug(debug)
        self.fields.grid.set_3d_symm()

    @staticmethod
    def calc_nmodes(lphi):
        nmodes=int(2**lphi/3)+1
        return nmodes

    @staticmethod
    def calc_nplanes(lphi):
        nplanes = 2**lphi
        return nplanes

    def get_dumptime(self):
        ''' Open the hdf5 dumpfile read the dump time and dumpstep
        '''
        with h5py.File(self.dumpfile, 'r') as h5file:
            try:
                self.time=h5file["dumpTime"].attrs['vsTime']
                self.step=int(h5file["dumpTime"].attrs['vsStep'])
            except:
              print(f"Error reading time or step in {self.dumpfile}")
              raise

    def setUpIntegrand(self):
        ''' This function calcualtes the weights for integration
            using normal quadrature'''
        ndim=self.fields.grid.rzp.ndim
        if ndim ==3:
            rz=self.fields.grid.rzp
            xy=self.fields.grid.xy
        elif ndim==4:
            if not self.fields.grid.symm_3d:
                print("integration assuming symmetric grid")
                raise ValueError
            rz=self.fields.grid.rzp[:,:,:,0]
            xy=self.fields.grid.xy[:,:,:,0]
        else:
            #first dim is RZ component
            print("integration require ndim=3,4")
            raise ValueError
        self.intWeight=np.zeros_like(rz[0])
        for jj in range(np.shape(self.intWeight)[1]):
            for ii in range(np.shape(self.intWeight)[0]):
                if not np.isnan(xy[:,ii,jj]).any():
                    horz='m'
                    vert='m'
                    if ii==0: #logic ignore degenerate cases /\ and \/
                        horz='l'
                    elif ii== np.shape(self.intWeight)[0]-1:
                        horz='r'
                    elif np.isnan(xy[:,ii-1,jj]).any():
                        horz='l'
                    elif np.isnan(xy[:,ii+1,jj]).any():
                        horz='r'
                    if jj==0:
                        vert='b'
                    elif jj== np.shape(self.intWeight)[1]-1:
                        vert='t'
                    elif np.isnan(xy[:,ii,jj-1]).any():
                        vert='b'
                    elif np.isnan(xy[:,ii,jj+1]).any():
                        vert='t'
                    if horz=='l':
                        if vert=='b':
                            #bottom left
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj])
                        elif vert=='t':
                            #top left
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii,jj])* \
                                (rz[1,ii,jj]-rz[1,ii,jj-1])
                        else:
                            #left middle
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj-1])
                    elif horz=='r':
                        if vert=='b':
                            #bottom right
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj])
                        elif vert=='t':
                            #top right
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj]-rz[1,ii,jj-1])
                        else:
                            #middle right
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj-1])
                    else:
                        if vert=='b':
                            #bottom middle
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj])
                        elif vert=='t':
                            #top middle
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj]-rz[1,ii,jj-1])
                        else:
                            #middle middle
                            self.intWeight[ii,jj]=0.25*rz[0,ii,jj]* \
                                (rz[0,ii+1,jj]-rz[0,ii-1,jj])* \
                                (rz[1,ii,jj+1]-rz[1,ii,jj-1])
        return

    def integrateSurface(self,integrand,twopi=True):
        #first check to see if integrand is set up.
        if not self.intSet:
            self.setUpIntegrand()
        if twopi:
            phifac=np.pi*2.0
        else:
            phifac=1.0
        if type(integrand) == np.ndarray:
            integral=phifac*np.tensordot(integrand,self.intWeight,axes=([0,1],[0,1]))
        else:
            print(f"Integrand of type {type(integrand)} is not supported")
            raise TypeError
        return integral

    def integrateEnergy(self,plot_integrand=False):
        self.fields.energyDensity()
        for key, integrand in self.fields.energyDict.items():
            integral=self.integrateSurface(integrand)
            print(key,integral)
            self.energyDict[key]=integral
            if plot_integrand:
                #todo This needs work
                self.plotIntegrand(integrand,imode=1)

    def integratePowerFlux(self,plot_integrand=False):
        self.fields.powerFlux()
        for key, integrand in self.fields.powerFluxDict.items():
            integral=self.integrateSurface(integrand)
            print(key,integral)
            if plot_integrand:
                #todo This needs work
                self.plotIntegrand(integrand,imode=1)
            self.powerDict[key]=integral


    def plotIntegrand(self,integrand,imode=None,title="Power Transfer"):
        ndim=self.fields.grid.rzp.ndim
        if ndim ==3:
            rz=self.fields.grid.rzp
            xy=self.fields.grid.xy
        elif ndim==4:
            if not self.fields.grid.symm_3d:
                print("integration assumes symmetric grid")
                raise ValueError
            rz=self.fields.grid.rzp[:,:,:,0]
            xy=self.fields.grid.xy[:,:,:,0]
        else:
            #first dim is RZ component
            print("integration require ndim=3,4")
            raise ValueError
        if imode==None:
            pass
        else:
            figsize=[6,6]
            fig, ax = plt.subplots(figsize=figsize)
            ax.set_aspect('equal')
            ax.set(title=title )#r"$v\cdot \nabla p$")
            plt.contourf(rz[0],rz[1],integrand[...,imode])
            plt.show()

    def analyze_power(self,grid='close',npts=512,lphi=5,nonlin_order=2):
        if grid == 'close':
            rmin=1.15
            rmax=2.3
            zmin=-1.25
            zmax=1.0
        elif grid == 'd3d':
            rmin=0.8
            rmax=2.5
            zmin=-1.5
            zmax=1.5
        else:
            print(f"Grid = {grid} is not recognized")
            raise ValueError
        #npts=128 #debug

        self.fields.set_method("powerflux")
        self.set_3dgrid(rmin,rmax,zmin,zmax,npts,npts,lphi,nonlin_order)
        self.integrateEnergy()
        self.integratePowerFlux()
        volumeInt=np.ones_like(self.fields.grid.rzp[0,:,:,0])
        self.volume=self.integrateSurface(volumeInt)

    def print_integrals(self):
        for key, integral in self.energyDict.items():
            print(key, integral)
        for key, integral in self.powerDict.items():
            print(key, integral)
        print('volume ', self.volume)

    def clean_up(self):
        ''' Clean up fields to reduce memory foot print'''
        self.fields.clean_up()
        self.grid=None
        self.intSet=False
        self.intWeight=np.empty([1])


###############################################################################
# TESTING NOTES
#############################
# Test 1, integrate from r=1.5 to 2, z=-1,1
#         everypoint in domain, and volume is 1.75x2xpi
#         v=10.99557428756
#         test gets right answer for 30 pts to high degree

# Test two, integrate volume from r=0.5 to 2.5, z=-1.5 to 1
#           some points are outside domain
#           pts            volume
#           30              35.1007188376975
#           100             37.215318055511986
#           200             37.68293725920001
#           300             37.826641583095295
#           500             37.96294543920083
#           1000            38.05936996...
#           extrapolation yeilds a volume of 38.14
#           300  pts error 0.82%
#           500  pts error 0.46%
#           1000 pts error 0.21% but slow
#
#           based on these results I trust volume integral
#           300x300 seems a good compromise between speed and accuracy
#           base on volume data
#           should verify with magnetic energy (test 3)

#Test two, integrate magentic energy from r=0.5 to 2.5, z=-1.5 to 1
#           some points are outside domain, compare n=0 (excludes eq)
#           n=1 and n=5 energies, Comments these results of n=1 and
#           n=5 divide the energy by 2. But for n=/0 this should not
#           be done.
#
#           pts     E_0             E_1             E_5
#           30      4.77153322e+02  7.88561043e+01  2.11200084e-01
#           50      5.45079600e+02  8.14135131e+01  2.26558515e-01
#           100     5.33019371e+02  8.25029029e+01  2.22909222e-01
#           300     5.36869913e+02  8.26732084e+01  2.25595720e-01
#           500     5.37585636e+02  8.27392466e+01  2.25714039e-01
#
#           Note extrapolation for n=0,5 isn't great (R^2~.5,.6)
#           Errors of (n=0,1,5)
#           pts 30(13%,5%,7%)
#           pts 50(.6%,2%,.5%)
#           pts 100(3%,1%,2%)
#           pts 300(2%,.8%,1%
#           pts 500(2%,.7%,.9%)
#
#           Conclusion, 100-300 points is good enought for 1% error
