#!/usr/bin/env python3
#
#profiles is a class for calculating 1D profiles
# using the flux surface integration
#
#
import f90nml
import eval_nimrod as eval
#import field_class as fc
import plot_nimrod as pn
import fsa
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d,splev,UnivariateSpline,griddata
from scipy.fft import fft, ifft
import os
import h5py
import sys
import numpy as np
import pickle
import ntm_fields
import time
import nim_timer as timer

class ntmstep:
    def __init__(self,dumpfile,nimrodin):
        #start with dump file and nimrod.in info
        self.dumpfile=dumpfile
        self.nimrodin=nimrodin
        self.fields=ntm_fields.ntmfields(dumpfile,nimrodin)
        self.time=None
        self.step=None
        self.eval=None
        self.grid=None
        self.mode=None #grid, fsa
        #field data
        self.fielddict={}
        # next include info on how fsa's were performed
        self.mmax=None
        self.ifour=[]
        self.nfour=None

        self.fsa=False
        self.dvar_dict={}
        self.fsa_dict={}
        self.interp_dict={}
        self.setprofiles=False
        self.raw_bcmn=np.empty([1])
        self.raw_bsmn=np.empty([1])
        self.raw_bmn_amp=np.empty([1])
        self.raw_bmn_phase=np.empty([1])

        #finally end with the profiles
        self.psin=np.empty([1])
        self.psi=np.empty([1])
        self.rhon=np.empty([1])
        self.q=np.empty([1])
        self.qpsi=np.empty([1])
        self.vprime=np.empty([1])
        self.bmn=np.empty([1])
        self.bcmn=np.empty([1])
        self.bsmn=np.empty([1])


    def dump(self,file):
        pickle.dump(self.dumpfile,file)
        pickle.dump(self.nimrodin,file)
        pickle.dump(self.time,file)
        pickle.dump(self.step,file)
        # next include info on how fsa's were performed
        pickle.dump(self.mmax,file)
        pickle.dump(self.ifour,file)
        pickle.dump(self.nfour,file)
        pickle.dump(self.fsa,file)
        pickle.dump(self.dvar_dict,file)
        pickle.dump(self.fsa_dict,file)
        pickle.dump(self.raw_bcmn,file)
        pickle.dump(self.raw_bsmn,file)
        pickle.dump(self.raw_bmn_amp,file)
        pickle.dump(self.raw_bmn_phase,file)


    def load(self,file):
        print(file)
        with open(file,'rb') as pickle_file:
            self.dumpfile=pickle.load(pickle_file)
            self.nimrodin=pickle.load(pickle_file)
            self.time=pickle.load(pickle_file)
            self.step=pickle.load(pickle_file)
            # next include info on how fsa's were performed
            self.mmax=pickle.load(pickle_file)
            self.ifour=pickle.load(pickle_file)
            self.nfour=pickle.load(pickle_file)
            self.fsa=pickle.load(pickle_file)
            self.dvar_dict=pickle.load(pickle_file)
            self.fsa_dict=pickle.load(pickle_file)
            self.raw_bcmn=pickle.load(pickle_file)
            self.raw_bsmn=pickle.load(pickle_file)
            self.raw_bmn_amp=pickle.load(pickle_file)
            self.raw_bmn_phase=pickle.load(pickle_file)

        if self.fsa == True:
            self.interpolate_fsa()

    def interpolate_fsa(self):

        rhomin=np.min(self.dvar_dict['rhon'])
        rhomax=np.max(self.dvar_dict['rhon'])
        print(self.dvar_dict['rhon'].shape,self.raw_bcmn.shape)
        self.rhon = np.linspace(rhomin,rhomax,300,endpoint=True)
        self.psin=interp1d(self.dvar_dict['rhon'],self.dvar_dict['psin'],kind='cubic')(self.rhon)
        self.psi =interp1d(self.dvar_dict['rhon'],self.dvar_dict['psi'] ,kind='cubic')(self.rhon)
        self.vprime=np.pi*2*interp1d(self.dvar_dict['rhon'],self.dvar_dict['vprime'],kind='cubic')(self.rhon)
        self.q=interp1d(self.dvar_dict['rhon'],self.dvar_dict['q'], kind='cubic')(self.rhon)
        self.bcmn=interp1d(self.dvar_dict['rhon'],self.raw_bcmn, kind='cubic')(self.rhon)
        self.bsmn=interp1d(self.dvar_dict['rhon'],self.raw_bsmn, kind='cubic')(self.rhon)
        self.bmn_amp =interp1d(self.dvar_dict['rhon'],self.raw_bmn_amp , kind='cubic')(self.rhon)
        self.bmn_phase =interp1d(self.dvar_dict['rhon'],self.raw_bmn_phase , kind='linear')(self.rhon)

        for ikey in self.fsa_dict:
            this_cmn,this_smn = self.fsa_dict[ikey]
            thismn_amp=np.sqrt(np.square(this_cmn)+np.square(this_smn))
            thismn_phase=np.arctan2(this_smn,this_cmn)

            cosmn=interp1d(self.dvar_dict['rhon'],this_cmn , kind='cubic')(self.rhon)
            sinmn=interp1d(self.dvar_dict['rhon'],this_smn , kind='cubic')(self.rhon)
            amp =interp1d(self.dvar_dict['rhon'],thismn_amp , kind='cubic')(self.rhon)
            phase =interp1d(self.dvar_dict['rhon'],thismn_phase , kind='linear')(self.rhon)
            print(ikey)
            self.interp_dict[ikey]=(cosmn,sinmn,amp,phase)


    def interpolate_psi(self):
        rhomin=np.min(self.dvar_dict['rhon'])
        rhomax=np.max(self.dvar_dict['rhon'])
        print(self.dvar_dict['rhon'].shape,self.raw_bcmn.shape)
        self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)
        self.psin=interp1d(self.dvar_dict['rhon'],self.dvar_dict['psin'],kind='cubic')(self.rhon)
        self.psi =interp1d(self.dvar_dict['rhon'],self.dvar_dict['psi'] ,kind='cubic')(self.rhon)
        self.vprime=np.pi*2*interp1d(self.dvar_dict['rhon'],self.dvar_dict['vprime'],kind='cubic')(self.rhon)
        self.q=interp1d(self.dvar_dict['rhon'],self.dvar_dict['q'], kind='cubic')(self.rhon)
        self.bcmn=interp1d(self.dvar_dict['rhon'],self.raw_bcmn, kind='cubic')(self.rhon)
        self.bsmn=interp1d(self.dvar_dict['rhon'],self.raw_bsmn, kind='cubic')(self.rhon)
        self.bmn_amp =interp1d(self.dvar_dict['rhon'],self.raw_bmn_amp , kind='cubic')(self.rhon)
        self.bmn_phase =interp1d(self.dvar_dict['rhon'],self.raw_bmn_phase , kind='linear')(self.rhon)


    def set_3dgrid(self,rmin,rmax,zmin,zmax,nr,nz,lphi,nonlin_order=2,debug=0):
        '''sets up a 3d grid, using non to determine the number of phi planes
           based on
        '''
        self.nmodes=self.calc_nmodes(lphi)
        self.nmax=self.nmodes-1
        nphi=nonlin_order*self.nmodes
        phimax=np.pi*2*(nphi-1)/nphi
        p1 = np.array([rmin, zmin, 0.0])
        p2 = np.array([rmax, zmin, 0.0])
        p3 = np.array([rmin, zmax, 0.0])
        p4 = np.array([rmin, zmin, phimax])
        rzp3d = pn.PlotNimrod.grid_3d_gen(p1, p2, p3, p4, nr, nz,nphi)
        self.fields.grid=eval.EvalGrid(rzp3d)
        self.fields.grid.set_debug(debug)
        self.fields.grid.set_3d_symm()

    def set_2dgrid(self,start,stop,npts,lphi,nonlin_order=2,debug=0):
        '''sets up a 3d grid, using non to determine the number of phi planes
           based on
        '''
        self.nmodes=self.calc_nmodes(lphi)
        self.nmax=self.nmodes-1
        nphi=nonlin_order*self.nmodes
        phimax=np.pi*2*(nphi-1)/nphi
        p1 = np.array([start[0], start[1], 0.0])
        p2 = np.array([stop[0], stop[1], 0.0])
        p3 = np.array([start[0], start[1], phimax])
        rzp2d = pn.PlotNimrod.grid_2d_gen(p1, p2, p3, npts,nphi)
        self.fields.grid=eval.EvalGrid(rzp2d)
        self.fields.grid.set_debug(debug)

    def set_fsagrid(self,r,z,lphi,nonlin_order):
        self.nmodes=self.calc_nmodes(lphi)
        self.nmax=self.nmodes-1
        nphi=nonlin_order*self.nmodes
        phimax=np.pi*2*(nphi-1)/nphi
        rzp=np.zeros([3,nphi])
        rzp[0,:]=r
        rzp[1,:]=z
        rzp[2,:]=np.linspace(0,phimax,nphi)
        return rzp

    def analyze(self):
        rmin=0.8
        rmax=2.5
        zmin=-1.5
        zmax=1.5
        npts=300
        lphi=5
        nonlin_order=2

        self.set_3dgrid(rmin,rmax,zmin,zmax,npts,npts,lphi,nonlin_order)
        #self.set_fsagrid(1,1,lphi,nonlin_order)
        self.fields.induction()
        divpie_fft=self.fields.fft(self.fields.edict['divpie'])
        figsize=[6,6]
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_aspect('equal')
        ax.set(title=r"$|E_\pi|_R$")
        pn.PlotNimrod.plot_scalar_plane(self.fields.grid.rzp[:,:,:,0], np.abs(divpie_fft[0,:,:,1]),ax=ax)


        fig, ax = plt.subplots(figsize=figsize)
        ax.set_aspect('equal')
        ax.set(title=r"$|E_\pi|_Z$")
        pn.PlotNimrod.plot_scalar_plane(self.fields.grid.rzp[:,:,:,0], np.abs(divpie_fft[1,:,:,1]),ax=ax)
        #  self.continuity()
        #  self.momentum()
        #  self.temperature()
        #  self.induction()

    def find_comp_boundary(self,inpoint,outpoint,tol=1e-8):
        index=0
        maxit=100
        ntest=self.fields.eval.eval_field('n',outpoint,dmode=0,eq=1)
        if ntest ==ntest:
            '''return if outpoint is in domain'''
            return outpoint
        fst=np.copy(inpoint)
        lst=np.copy(outpoint)
        print(fst,inpoint)
        print(lst,outpoint)
        dist=np.linalg.norm(lst-fst,2)
        if dist<tol:
            return fst
        while True:
            tst=(fst+lst)/2.0
            print(tst)
            ntest=self.fields.eval.eval_field('n',tst,dmode=0,eq=1)
            if ntest ==ntest:
                fst=tst
            else:
                lst=tst
            dist=np.linalg.norm(lst-fst,2)
            print(dist)
            if dist<tol:
                break
            index+=1
            if index>maxit:
                print("no convergence")
                break
        print(fst)
        print(lst)
        return fst

    def plot_line(self,line,field,flabel=None,xvar=None,title='Title',
                  xlabel='x', ylabel='y',fig_size = [12,6.75],qlist=None):
        print(np.shape(line))
        if xvar==None:
            xvar = np.sqrt((line[0, :] - line[0, 0]) ** 2 +
                          (line[1, :] - line[1, 0]) ** 2 +
                          (line[2, :] - line[2, 0]) ** 2)

        fig, ax = plt.subplots(figsize=fig_size)

        ax.plot(xvar, field, alpha=0.7, label=flabel)
        if qlist is not None:
            for ii in qlist:
                ax.axvline(ii[1],label=ii[0],ls=':')
            ax.legend(loc=1)

        ax.set(ylabel=ylabel, xlabel=xlabel,title=title)


        ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
                    useOffset=None, useLocale=None, useMathText=True)

        plt.tight_layout()
        plt.show()
        return

    def analyze_radial(self,dir='mid'):
        npts=500
        lphi=5
        nonlin_order=2
        rzo=np.array([1.76821,-0.0188439,0.0])
        oPoint=self.find_pf_null(self.fields.eval, rzo, flag=0)
        rzx=np.array([1.27,-1.14,0.0])
        xPoint=self.find_pf_null(self.fields.eval, rzx, flag=0)
        #find rmax
        rzmax=np.copy(oPoint)
        rzmax[0]=3.0
        rzmax=self.find_comp_boundary(oPoint,rzmax)
        self.set_2dgrid(oPoint,rzmax,npts,lphi,nonlin_order)
        self.fields.induction()
        divpie_fft=self.fields.fft(self.fields.edict['divpie'])
        divpie0_fft=self.fields.fft(self.fields.edict['divpie0'])
        divpie1_fft=self.fields.fft(self.fields.edict['divpie1'])
        divpieb_fft=self.fields.fft(self.fields.edict['divpieb'])



        pn.PlotNimrod.plot_scalar_line(self.fields.grid.rzp[:,:,0], self.fields.fielddict['fsa_beq2'].data[:,0])


        q=np.empty(npts)
        q[:]=np.NaN
        length=np.sqrt(
        (self.fields.grid.rzp[0,:,0] - self.fields.grid.rzp[0,0,0]) ** 2 +
        (self.fields.grid.rzp[1,:,0] - self.fields.grid.rzp[1,0,0]) ** 2 +
        (self.fields.grid.rzp[2,:,0] - self.fields.grid.rzp[2,0,0]) ** 2)
        for ii in range(npts):
            yvar, contours = fsa.FSA(self.fields.eval, self.fields.grid.rzp[:,ii,0], self.q_fsa, neq=1, \
            nsurf=1,flag=1,normalize=False,rzp=self.fields.grid.rzp[:,ii,0])
            if yvar!=yvar:
                break
            q[ii]=yvar/(np.pi*2.0)

        iend=npts-1
        while np.isnan(q[iend]):
            iend -= 1
        lenofq=interp1d(q[:iend], length[:iend], kind='cubic',fill_value="extrapolate")
        qlist=[]
        qlist.append(('q=-2',lenofq(-2)))
        qlist.append(('q=-3',lenofq(-3)))
        qlist.append(('q=-4',lenofq(-4)))

        self.plot_line(self.fields.grid.rzp[:,:,0], np.abs(divpie_fft[0,:,1]),
        title="Radial n=1 neoclassical E",xlabel="Distance from o-point", ylabel=r"$|E_\Pi|_R$",
        qlist=qlist)

        self.plot_line(self.fields.grid.rzp[:,:,0], np.abs(divpie_fft[1,:,1]),
        title="Radial n=1 neoclassical E",xlabel="Distance from o-point", ylabel=r"$|E_\Pi|_Z$",
        qlist=qlist)

        self.plot_line(self.fields.grid.rzp[:,:,0], np.abs(divpie0_fft[1,:,0]),
        title="Vertical n=1 E factor",xlabel="Distance from o-point", ylabel=r"$|F_E|_Z$",
        qlist=qlist)

        self.plot_line(self.fields.grid.rzp[:,:,0], np.abs(divpie1_fft[1,:,0]),
        title="Vertical n=1 neoclassical stress",xlabel="Distance from o-point [m]", ylabel=r"$ |\Pi_e|_Z$",
        qlist=qlist)

        self.plot_line(self.fields.grid.rzp[:,:,0], np.abs(divpieb_fft[1,:,0]),
        title="Vertical n=1 neoclassical factor",xlabel="Distance from o-point", ylabel=r"$|Fl_\Pi|_Z$",
        qlist=qlist)

    @staticmethod
    def calc_nmodes(lphi):
        nmodes=int(2**lphi/3)+1
        return nmodes

    def get_m_index(self,m):
        '''Return the index for the given m
          Return None if m is out of range'''
        if self.mmax==None:
            write("mmax has not be set in get_m_index")
            raise
        else:
            if m>self.mmax:
                return None
            elif m<-1*self.mmax:
                return None
            else:
                return m+self.mmax

    def dummy_fsa(self,rzc,y,dy,evalnimrod,fargs):
        '''
        Dummy integrand for complex fsa, this is used to get v' and q
        without running a true fsa
        Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
        dy(0)=dl/deta or d eta/dl
        dy(1)=dr/deta or dr/dl
        dy(2)=1/bdgth :v'
        dy(3)=dq: q
        '''
        dy[4]=1.0
        return dy

    def q_fsa(self,rzc,y,dy,evalnimrod,fargs):
        '''
        Simple integrand for this is used to get q when running with 1 rzp
        Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
        dy(0)=dl/deta or d eta/dl
        dy(1)=dr/deta or dr/dl
        dy(2)=1/bdgth :v'
        dy(3)=dq: q
        '''
        dy[4]=dy[3]
        return dy

    def surfmn_int(self,rzc,y,dy,evalnimrod,fargs):
        '''
        Integrand for fluxsurface integration
        Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
        dy(0)=dl/deta or d eta/dl
        dy(1)=dr/deta or dr/dl
        dy(2)=1/bdgth
        dy(3)=dq
        dy(4)=dtheta
        '''
        addpert=fdict.get("addpert",True)
        #self.mmax=fargs.get("mmax")
        b0=np.array(fsa.get_b0(evalnimrod,rzc,addpert))
        b = evalnimrod.eval_field('b', rzc, dmode=0)
        rr =rzc[0]
        q=self.qpsi(fargs.get('psi'))
        jac=rr*q/b0[2]
        dy[4]=dy[2]/jac #dtheta
        for ii, im in enumerate(self.ifour):
            oset = ii * (4*self.mmax+1)
            reb=np.real(b[:,im+1])
            imb=np.imag(b[:,im+1])
            rBePsi=rr*(reb[1]*b0[0]-reb[0]*b0[1])
            iBePsi=rr*(imb[1]*b0[0]-imb[0]*b0[1])
            for im in range(self.mmax):
                nmth=-(self.mmax-im)*y[4] #negative m theta
                pmth=(im+1)*y[4] #positive m theta
                dy[5+im+oset]=(rBePsi*np.cos(nmth)-iBePsi*np.sin(nmth))*dy[2]
                dy[6+self.mmax+im+oset]=(rBePsi*np.cos(pmth)-iBePsi*np.sin(pmth))*dy[2]
                dy[6+2*self.mmax+im+oset]=-(rBePsi*np.sin(nmth)+iBePsi*np.cos(nmth))*dy[2]
                dy[6+3*self.mmax+im+oset]=-(rBePsi*np.sin(pmth)+iBePsi*np.cos(pmth))*dy[2]
            dy[5+self.mmax+oset]=rBePsi*dy[2]
        return dy

    @timer.timer_func
    def induction_int(self,rzc,y,dy,evalnimrod,fargs):
        ifour=fargs.get("ifour")
        mmax=fargs['mmax']
        nfour=len(fargs['ifour'])
        addpert=fargs.get("addpert",True)
        lphi=5
        nonlin_order=2
        grid=self.set_fsagrid(rzc[0],rzc[1],lphi,nonlin_order)
        self.fields.fielddict={}
        self.fields.ndict={}
        self.fields.edict={}
        self.fields.eval_b(grid=grid,fft=True)
        self.fields.induction(grid=grid)
        if addpert:
            b0=self.fields.fielddict['b0'].data[:,0]
        else:
            b0=self.fields.fielddict['beq'].data[:,0]
        rr =rzc[0]
        #print(fargs.get('psi'))
        q=self.qpsi(fargs.get('psi'))
        jac=rr*q/b0[2]
        dy[4]=dy[2]/jac
        nfour=ifour[0]
        #get br first
        ii=0
        #    for ii, im in enumerate(ifour):
        oset = ii * (4*mmax+1)
        #todo test
#        reb=np.real(self.fields.fielddict['bpert'].data[:,nfour])
#        imb=np.imag(self.fields.fielddict['bpert'].data[:,nfour])
        reb=np.real(self.fields.fielddict['bfour'][:,nfour])
        imb=np.imag(self.fields.fielddict['bfour'][:,nfour])

        rBePsi=rr*(reb[1]*b0[0]-reb[0]*b0[1])
        iBePsi=rr*(imb[1]*b0[0]-imb[0]*b0[1])
        for im in range(mmax):
            nmth=-(mmax-im)*y[4] #negative m theta
            pmth=(im+1)*y[4] #positive m theta
            dy[5+im+oset]=(rBePsi*np.cos(nmth)-iBePsi*np.sin(nmth))*dy[2]
            dy[6+mmax+im+oset]=(rBePsi*np.cos(pmth)-iBePsi*np.sin(pmth))*dy[2]
            dy[6+2*mmax+im+oset]=-(rBePsi*np.sin(nmth)+iBePsi*np.cos(nmth))*dy[2]
            dy[6+3*mmax+im+oset]=-(rBePsi*np.sin(pmth)+iBePsi*np.cos(pmth))*dy[2]
        dy[5+mmax+oset]=rBePsi*dy[2]
        #print(self.fields.dbdtdict.keys())
        for ii, ikey in enumerate(self.fields.dbdtdict):
            oset = (ii+1) * (4*mmax+1)
            thisFi=self.fields.fft(self.fields.dbdtdict[ikey])[:,nfour]
            reFi=np.real(thisFi)
            imFi=np.imag(thisFi)
            #        print(ikey,reFi,imFi)
            rFiPsi=rr*(reFi[1]*b0[0]-reFi[0]*b0[1])
            iFiPsi=rr*(imFi[1]*b0[0]-imFi[0]*b0[1])
            #        print(rFiPsi,iFiPsi)
            for im in range(mmax):
                nmth=-(mmax-im)*y[4] #negative m theta
                pmth=(im+1)*y[4] #positive m theta
                dy[5+im+oset]=(rFiPsi*np.cos(nmth)-iFiPsi*np.sin(nmth))*dy[2]
                dy[6+mmax+im+oset]=(rFiPsi*np.cos(pmth)-iFiPsi*np.sin(pmth))*dy[2]
                dy[6+2*mmax+im+oset]=-(rFiPsi*np.sin(nmth)+iFiPsi*np.cos(nmth))*dy[2]
                dy[6+3*mmax+im+oset]=-(rFiPsi*np.sin(pmth)+iFiPsi*np.cos(pmth))*dy[2]
            dy[5+mmax+oset]=rFiPsi*dy[2]
        return dy

    @timer.timer_func
    def psi_int(self,rzc,y,dy,evalnimrod,fargs):
        ifour=fargs.get("ifour")
        mmax=fargs['mmax']
        nfour=len(fargs['ifour'])
        addpert=fargs.get("addpert",True)
        lphi=5
        nonlin_order=2
        grid=self.set_fsagrid(rzc[0],rzc[1],lphi,nonlin_order)
        self.fields.fielddict={}
        self.fields.ndict={}
        self.fields.edict={}
        self.fields.eval_b(grid=grid,fft=True)
        if addpert:
            b0=self.fields.fielddict['b0'].data[:,0]
        else:
            b0=self.fields.fielddict['beq'].data[:,0]
        rr =rzc[0]
        #print(fargs.get('psi'))
        q=self.qpsi(fargs.get('psi'))
        jac=rr*q/b0[2]
        dy[4]=dy[2]/jac
        nfour=ifour[0]

        ii=0
        #    for ii, im in enumerate(ifour):
        oset = ii * (4*mmax+1)
        #print(self.fields.fielddict['bfour'])
        #print(self.fields.fielddict['bfour'].shape)
        reb=np.real(self.fields.fielddict['bfour'][:,nfour])
        #print(self.fields.fielddict['bfour'].data[:,:])
        imb=np.imag(self.fields.fielddict['bfour'][:,nfour])
        rBePsi=rr*(reb[1]*b0[0]-reb[0]*b0[1])
        iBePsi=rr*(imb[1]*b0[0]-imb[0]*b0[1])
        for im in range(mmax):
            nmth=-(mmax-im)*y[4] #negative m theta
            pmth=(im+1)*y[4] #positive m theta
            dy[5+im+oset]=(rBePsi*np.cos(nmth)-iBePsi*np.sin(nmth))*dy[2]
            dy[6+mmax+im+oset]=(rBePsi*np.cos(pmth)-iBePsi*np.sin(pmth))*dy[2]
            dy[6+2*mmax+im+oset]=-(rBePsi*np.sin(nmth)+iBePsi*np.cos(nmth))*dy[2]
            dy[6+3*mmax+im+oset]=-(rBePsi*np.sin(pmth)+iBePsi*np.cos(pmth))*dy[2]
        dy[5+mmax+oset]=rBePsi*dy[2]
        return dy

    def get_rho_q(self,q):
        try:
            return interp1d(self.q,self.rhon, kind='cubic',fill_value="extrapolate")(q)
        except:
            print(f"The safety factor {q} is not it the domain")
            raise

    def get_field_rho(self,field,rhon):
        try:
            return interp1d(self.rhon,field, kind='cubic')(rhon)
        except:
            print(f"Problem evaluitng field at rhon={rhon}")
            raise

    @timer.timer_func
    def calculate_induction(self,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
        self.fields.set_method("induction")
        ifour=fargs.get("ifour")
        mmax=fargs['mmax']
        nfour=len(fargs['ifour'])
        self.mmax=mmax
        self.ifour=ifour
        self.nfour=nfour
        self.fsa=True
        fargs['rtol']=1.e-8
      #
        dvar, yvar, contours = fsa.FSA(self.fields.eval, rzo, self.dummy_fsa, 1, \
          nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
          fargs=fargs)

        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1
        #unevaluated interpoate
        self.qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic',fill_value="extrapolate")
        #call induction at opoint to get number of fields
        self.fields.induction(grid=rzo)
        nterms=len(self.fields.dbdtdict)+1
        #neq 1 for dtheta
        #4*self.mmax+1 for each field per toroidal mode
        #to keep things simple let's assume 1 toroidal mode
        #include 1 term for B itself
        neq =1+nterms*1*(4*mmax+1)
        dvar,yvar,contours = fsa.FSA(self.fields.eval, rzo, self.induction_int, neq, \
          nsurf=nsurf,depvar='eta', dpow=0.5,rzx=rzx,flag=eqflag,normalize=False,\
          **fargs)
        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1

        rhomin=np.min(dvar[1,:iend])
        rhomax=np.max(dvar[1,:iend])
        self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)

        bmn_amp=np.zeros([1,2*mmax+1,iend])
        bmn_phase=np.zeros([1,2*mmax+1,iend])
        bcmn=np.zeros([1,2*mmax+1,iend])
        bsmn=np.zeros([1,2*mmax+1,iend])

        ii=0
        oset = ii * (4*mmax+1)
        bcmn[ii,:,:]= yvar[1+oset:2*mmax+2+oset,:iend]*(np.pi*2.0)
        bsmn[ii,0:mmax,:]=\
            yvar[2+2*mmax+oset:2+3*mmax+oset,:iend]*(np.pi*2.0)
        bsmn[ii,mmax+1:2*mmax+1,:]=\
            yvar[2+3*mmax+oset:2+4*mmax+oset,:iend]*(np.pi*2.0)
        bmn_amp=np.sqrt(np.square(bcmn)+np.square(bsmn))
        bmn_phase=np.arctan2(bsmn,bcmn)
        this_dict={}
        self.fsa_dict={}
        for ii, ikey in enumerate(self.fields.dbdtdict):
            oset = (ii+1) * (4*mmax+1)
            #oset = (ii+1)
            this_amp=np.zeros([1,2*mmax+1,iend])
            this_phase=np.zeros([1,2*mmax+1,iend])
            this_cmn=np.zeros([1,2*mmax+1,iend])
            this_smn=np.zeros([1,2*mmax+1,iend])
            this_cmn[0,:,:]= yvar[1+oset:2*mmax+2+oset,:iend]*(np.pi*2.0)
            this_smn[0,0:mmax,:]=\
              yvar[2+2*mmax+oset:2+3*mmax+oset,:iend]*(np.pi*2.0)
            this_smn[0,mmax+1:2*mmax+1,:]=\
              yvar[2+3*mmax+oset:2+4*mmax+oset,:iend]*(np.pi*2.0)
            print(ikey)
            self.fsa_dict[ikey]=(this_cmn,this_smn)

        #dvars
        self.dvar_dict={}
        self.dvar_dict['psin']=dvar[0,:iend]
        self.dvar_dict['rhon']=dvar[1,:iend]
        self.dvar_dict['psi']=dvar[2,:iend]
        self.dvar_dict['vprime']=dvar[6,:iend]
        self.dvar_dict['q']=dvar[7,:iend]

        self.raw_bcmn=bcmn
        self.raw_bsmn=bsmn
        self.raw_bmn_amp=bmn_amp
        self.raw_bmn_phase =bmn_phase

        self.interpolate_fsa()

        # neq=1+nterm*self.nfour*(4*self.mmax+1)
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(self.rhon,self.q)
        plt.show()

        fig =plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(self.rhon,self.bmn_amp[0,2,:])
        plt.show()

        fig =plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(self.rhon,self.bmn_phase[0,2,:])
        plt.show()

        return None


    @timer.timer_func
    def calculate_psi(self,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
        self.fields.set_method("induction")

        ifour=fargs.get("ifour")
        mmax=fargs['mmax']
        nfour=len(fargs['ifour'])
        self.mmax=mmax
        self.ifour=ifour
        self.nfour=nfour

        self.fsa=True
        fargs['rtol']=1.e-8
    #
        dvar, yvar, contours = fsa.FSA(self.fields.eval, rzo, self.dummy_fsa, 1, \
            nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
            fargs=fargs)

        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1
        #unevaluated interpoate
        self.qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic',fill_value="extrapolate")
        #call induction at opoint to get number of fields
        self.fields.induction(grid=rzo)
        nterms=len(self.fields.dbdtdict)+1
        #neq 1 for dtheta
        #4*self.mmax+1 for each field per toroidal mode
        #to keep things simple let's assume 1 toroidal mode
        #include 1 term for B itself
        neq =1+(4*mmax+1)
        #  neq=1+(4*mmax+1)
        dvar,yvar,contours = fsa.FSA(self.fields.eval, rzo, self.psi_int, neq, \
            nsurf=nsurf,depvar='eta', dpow=0.5,rzx=rzx,flag=eqflag,normalize=False,\
            **fargs)

        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1

        rhomin=np.min(dvar[1,:iend])
        rhomax=np.max(dvar[1,:iend])
        self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)

        bmn_amp=np.zeros([1,2*mmax+1,iend])
        bmn_phase=np.zeros([1,2*mmax+1,iend])
        bcmn=np.zeros([1,2*mmax+1,iend])
        bsmn=np.zeros([1,2*mmax+1,iend])

        ii=0
        oset = ii * (4*mmax+1)
        bcmn[ii,:,:]= yvar[1+oset:2*mmax+2+oset,:iend]*(np.pi*2.0)
        bsmn[ii,0:mmax,:]=\
            yvar[2+2*mmax+oset:2+3*mmax+oset,:iend]*(np.pi*2.0)
        bsmn[ii,mmax+1:2*mmax+1,:]=\
            yvar[2+3*mmax+oset:2+4*mmax+oset,:iend]*(np.pi*2.0)
        bmn_amp=np.sqrt(np.power(bcmn,2)+np.power(bsmn,2))
        bmn_phase=np.arctan2(bsmn,bcmn)

        #dvars
        self.dvar_dict={}
        self.dvar_dict['psin']=dvar[0,:iend]
        self.dvar_dict['rhon']=dvar[1,:iend]
        self.dvar_dict['psi']=dvar[2,:iend]
        self.dvar_dict['vprime']=dvar[6,:iend]
        self.dvar_dict['q']=dvar[7,:iend]

        self.raw_bcmn=bcmn
        self.raw_bsmn=bsmn
        self.raw_bmn_amp=bmn_amp
        self.raw_bmn_phase =bmn_phase

        self.interpolate_psi()

        # neq=1+nterm*self.nfour*(4*self.mmax+1)
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(self.rhon,self.q)
        plt.show()

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(dvar[1,:iend],bmn_amp[0,3,:])
        plt.show()

        fig =plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        for im in range(self.bmn_amp.shape[1]):
            conf=plt.plot(self.rhon,self.bmn_amp[0,im,:],label=f"index={im}")
        plt.legend()
        plt.show()

        fig =plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        for im in range(self.bmn_amp.shape[1]):
            conf=plt.plot(self.rhon,self.bcmn[0,im,:],label=f"cos index={im}")
        plt.legend()
        plt.show()

        fig =plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        conf=plt.plot(self.rhon,self.bmn_amp[0,3,:],label=f"index={3}")
        conf=plt.plot(self.rhon,self.bmn_amp[0,7,:],label=f"index={7}")
        plt.legend()
        plt.show()


        return None

    def calculate(self,rzo=None,rzx=None,nsurf=150,eqflag=0,fargs={},**kwargs):
        mi=kwargs.get("mi",3.3435860e-27)
        qe=kwargs.get("qe",1.609e-19)
        self.ifour=fargs.get("ifour")
        self.mmax=fargs['mmax']
        self.nfour=len(fargs['ifour'])
        self.setprofiles=True

        #first call to fsa is to calcualte q
        evalnimrod=eval.EvalNimrod(self.dumpfile,fieldlist='nvptbj')
        dvar, yvar, contours = fsa.FSA(evalnimrod, rzo, self.dummy_fsa, 1, \
          nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
          fargs=fargs)

        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1
        #unevaluated interpoate
        self.qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic')

        #second call to fsa is to calcualte b_ms ans psi_mn
        neq=1+self.nfour*(4*self.mmax+1)
        dvar,yvar,contours = fsa.FSA(evalnimrod, rzo, self.surfmn_int, neq, \
            nsurf=nsurf,depvar='eta', dpow=0.5,rzx=rzx,flag=eqflag,normalize=False,\
            fargs=fargs)

        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1

        bmn=np.zeros([self.nfour,2*self.mmax+1,iend])
        bcmn=np.zeros([self.nfour,2*self.mmax+1,iend])
        bsmn=np.zeros([self.nfour,2*self.mmax+1,iend])
        for ii in range(self.nfour):
            oset = ii * (4*self.mmax+1)
            bcmn[ii,:,:]= yvar[1+oset:2*self.mmax+2+oset,:iend]*(np.pi*2.0)
            bsmn[ii,0:self.mmax,:]=\
                yvar[2+2*self.mmax+oset:2+3*self.mmax+oset,:iend]*(np.pi*2.0)
            bsmn[ii,self.mmax+1:2*self.mmax+1,:]=\
                yvar[2+3*self.mmax+oset:2+4*self.mmax+oset,:iend]*(np.pi*2.0)
            bmn=np.sqrt(np.square(bcmn)+np.square(bsmn))
        rhomin=np.min(dvar[1,:iend])
        rhomax=np.max(dvar[1,:iend])
        self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)
        #dvars
        self.psin=interp1d(dvar[1,:iend], dvar[0,:iend], kind='cubic')(self.rhon)
        self.psi=interp1d(dvar[1,:iend], dvar[2,:iend], kind='cubic')(self.rhon)
        self.vprime=np.pi*2*interp1d(dvar[1,:iend], dvar[6,:iend], kind='cubic')(self.rhon)
        self.q=interp1d(dvar[1,:iend], dvar[7,:iend], kind='cubic')(self.rhon)

        self.bcmn=interp1d(dvar[1,:iend],bcmn, kind='cubic')(self.rhon)
        self.bsmn=interp1d(dvar[1,:iend],bsmn, kind='cubic')(self.rhon)
        self.bmn =interp1d(dvar[1,:iend],bmn , kind='cubic')(self.rhon)

    def get_b0(self,eval_nimrod,rzn,flag,abort=False):
        """
        Find b at a given point
        :param eval_nimrod: eval_nimrod class instance
        :param rzn: initial guess for poloidal field null
        :param flag: if 0 only use eq, if 1 add n=0 to eq
        :param abort: raise an exception if true and can't find b
        """
        b = eval_nimrod.eval_field('b', rzn, dmode=0)
        b0=np.real(b[:])
        if (abort and np.isnan(b0).any()):
            print(b)
            raise Exception('FSA_find_pf_null: Hit wall')
        return b0

    def find_pf_null(self,eval_nimrod, rzn, flag=0):
        """
        Find a poloidal field null
        :param eval_nimrod: eval_nimrod class instance
        :param rzn: initial guess for poloidal field null
        :param flag: if 0 only use eq, if 1 add n=0 to eq
        """
        rzn = np.array(rzn)
        maxsteps=1000
        it=0
        rtol=1.e-8
        drz0=0.125*rzn[0]
        while True:
            b = self.get_b0(eval_nimrod,rzn,flag,abort=False)
            norm0=np.sqrt(b[0]**2+b[1]**2)
            rvn=-rzn[0]*b[1]/norm0
            zvn= rzn[0]*b[0]/norm0
            drz=drz0*(1.0-float(it)/maxsteps)+rtol*rzn[0]
            while True:
                rr=rzn[0]+rvn*drz
                zz=rzn[1]+zvn*drz
                rzng=np.array([rr, zz, 0.0])
                b = self.get_b0(eval_nimrod,rzng,flag,abort=False)
                if not np.isnan(b).any():
                    norm=np.sqrt(b[0]**2+b[1]**2)
                    if (norm < norm0):
                        rzn[:]=rzng[:]
                        break
                rr=rzn[0]-rvn*drz
                zz=rzn[1]-zvn*drz
                rzng=np.array([rr, zz, 0.0])
                b = self.get_b0(eval_nimrod,rzng,flag,abort=False)
                if not np.isnan(b).any():
                    norm=np.sqrt(b[0]**2+b[1]**2)
                    if (norm < norm0):
                        rzn[:]=rzng[:]
                        break
                drz=drz/2.0
                if (drz/rzn[0] < rtol):
                    return rzn # done
            it=it+1
            if it>=maxsteps:
                raise Exception('FSA find_pf_null: No convergence')
                return None

    def plot(self,pargs={}):
        for im,imode in enumerate(self.ifour):
            self.plot_radial(im,imode,pargs)
            self.plot_surfmn(im,imode,pargs)

    def plot_radial(self,ii,imode,pargs={}):
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        title=f"$\psi$(n={int(imode)}) at {self.time*1000:.3f}ms"
        ylabel=f"$\psi_m$ [mWb]"
        colorlist = list(mcolors.TABLEAU_COLORS)
        xlabel=r'$\rho_N$'
        fontsize=18
        if imode==1:
            mlist=range(-4,1)
        elif imode==2:
            mlist=range(-6,-1)
        else:
            mstart=-2*imode
            mlist=range(mstart,mstart+imode+1)
        if 'mlists' in pargs:
            if ii<len(pargs['mlists'][ii]):
                mlist=pargs['mlists'][ii]

        for im,this_m in enumerate(mlist):
            this_i = self.get_m_index(this_m)
            if this_i!= None:
                mlbl = "m = " + str(this_m)
                tc=colorlist[im%len(colorlist)]
                ax.plot(self.rhon,self.bmn[ii,this_i,:]*1000, color=tc, label=mlbl)
        try:
            qlist=pargs['qlists'][ii]
        except:
          if imode==1:
            qlist=[-4,-3,-2]
          elif imode==2:
            qlist=[-4,-3,-2.5,-2,-1.5]
          elif imode==3:
            qlist=[-3,-2.33, -2,-1.67,-1.33]
          elif imode==4:
            qlist=[-3,-2,-1.75,-1.5,-1.25]
          elif imode==5:
            qlist=[-3,-2,-1.8,-1.6,-1.4,-1.2]
          else:
            qlist=[-4,-3,-2]

        for iq,qq in enumerate(qlist):
          try:
            irho = self.get_rho_q(qq)
            qlbl = f"q = {qq:.2f}"
            tc=colorlist[iq]
            ax.axvline(irho,ls=':',color=tc, label=qlbl)
          except:
            print(f"q={qq:.2f} is not in the domain")
        ax.axhline(0,ls='-',c='k')
        ax.legend(loc=0,frameon=True,fontsize=fontsize)
        plt.title(title,fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        plt.ylabel(ylabel,fontsize=fontsize)
        plt.tight_layout()
        plt.show()

    def plot_surfmn(self,im,imode,surfmn,pargs={}):
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        # Set titles and labels
        title=f"$\psi$(n={int(imode)}) at {self.time*1000:.3f}ms"
        # set contour levels, i could generalize this further if needed
        levels=301
        vmax=np.amax(self.bmn[im,:,:])*1000
        levels=np.linspace(0,vmax,301)
        cbar_ticks=np.linspace(0,vmax,11)
        # Update plot based on keys in kwargs
        xlabel="Poloidal Mode Number m"
        fontsize=18
        # set up mrange()
        qmin=np.amin(self.q)
        qmax=np.amax(self.q)
        mrange=np.linspace(qmin,qmax)
        #create the surfmn plot
        plt.set_cmap('nipy_spectral')
        m=range(-self.mmax,self.mmax+1)
        mv, rv = np.meshgrid(m, self.rhon, sparse=False, indexing='ij')
        conf=plt.contourf(mv,rv,np.clip(self.bmn[im,:,:]*1000,0,None),levels=levels,vmax=vmax)
        plt.plot(imode*mrange,self.get_rho_q(mrange),c='w')
        plt.title(title,fontsize=fontsize)
        plt.ylabel(r'$\rho_N$',fontsize=fontsize)
        plt.xlabel(xlabel,fontsize=fontsize)
        cbar=fig.colorbar(conf,ticks=cbar_ticks)
        plt.xlim(-self.mmax,self.mmax)
        plt.show()

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

    def plot_scalar(self,rr,zz,field):
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.set_cmap('nipy_spectral')
        conf=plt.contourf(rr,zz,field,corner_mask=False)
        plt.show()

    @staticmethod
    def phase_shift(array,shift):
        array=array-shift
        array=np.where(array<-np.pi,array+np.pi,array)
        array=np.where(array>np.pi,array-np.pi,array)
        return array

    @staticmethod
    def phase_shift2(fcos,fsin,phase):
        ncosphase=np.cos(-phase)
        nsinphase=np.sin(-phase)
        inphase=fcos*ncosphase-fsin*nsinphase
        outphase=fcos*nsinphase+fsin*ncosphase
        return inphase,outphase

    def plot_fsa(self, key=None, **kwargs):
        qlist=[-2,-3,-4]
        fontsize=18
        colorlist = list(mcolors.TABLEAU_COLORS)
        titledict={}
        titledict['veqbeq']=r"$\nabla \times \left(V_0 \times B_0\right)_{2/1}$"
        titledict['vpbeq']=r"$\nabla \times \left(\tilde V \times B_0\right)_{2/1}$"
        titledict['veqbp']=r"$\nabla \times \left(V_0 \times \tilde B\right)_{2/1}$"
        titledict['vpbp']=r"$\nabla \times \left(\tilde V \times \tilde B\right)_{2/1}$"
        titledict['etajpert']=r"$\nabla \times \left(\eta \tilde J\right)_{2/1}$"
        titledict['divpie']=r"$\nabla \times \left(\frac{1}{ne}\nabla \cdot \pi_e\right)_{2/1}$"

        phase_21=interp1d(self.rhon,self.bmn_phase[0,3,:], kind='linear')(self.get_rho_q(-2))
        print(phase_21)
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        tc=colorlist[0%len(colorlist)]
        plt.title(r"$\psi_{2/1}$ amplitude")
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|\psi_{2/1}|$ [mWb/s]',fontsize=fontsize)
        ax.plot(self.rhon,self.bmn_amp[0,3,:]*1000, color=tc)
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.2f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.show()

        if key is None:
            for ikey in self.interp_dict:
                this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
        #'amp'
                fig = plt.figure(figsize=(10,8))
                ax=fig.add_subplot(111)
                plt.title(titledict.get(ikey,ikey))
                plt.plot(self.rhon,self.this_amp[0,3,:])
                plt.xlabel(r'$\rho_N$',fontsize=fontsize)
                plt.ylabel(r'$|F_{m/n}|$ [Wb/s]',fontsize=fontsize)
                for iq,qq in enumerate(qlist):
                    try:
                        irho = self.get_rho_q(qq)
                        qlbl = f"q = {qq:.2f}"
                        tc=colorlist[iq]
                        ax.axvline(irho,ls=':',color=tc, label=qlbl)
                    except:
                        print("q not found")
                plt.tight_layout()
                plt.show()
        else:
            pass
        return None


    def plot_fsa_phase(self, key=None, **kwargs):
        qlist=[-1.5,-2,-3,-4]
        fontsize=18
        colorlist = list(mcolors.TABLEAU_COLORS)
        titledict={}
        titledict={}
        titledict['veqbeq']=r"$\nabla \times \left(V_0 \times B_0\right)_{2/1}$"
        titledict['vpbeq']=r"$\nabla \times \left(\tilde V \times B_0\right)_{2/1}$"
        titledict['veqbp']=r"$\nabla \times \left(V_0 \times \tilde B\right)_{2/1}$"
        titledict['vpbp']=r"$\nabla \times \left(\tilde V \times \tilde B\right)_{2/1}$"
        titledict['etajpert']=r"$\nabla \times \left(\eta \tilde J\right)_{2/1}$"
        titledict['divpie']=r"$\nabla \times \left(\frac{1}{ne}\nabla \cdot \pi_e\right)_{2/1}$"
        sumlist=['etajpert','divpie']
        phase_21=interp1d(self.rhon,self.bmn_phase[0,3,:], kind='linear')(self.get_rho_q(-2))
        binphase,boutphase=self.phase_shift2(self.bcmn[0,3,:],self.bsmn[0,3,:],phase_21)
        print(phase_21)

        if "time" in kwargs:
            tmp=kwargs["time"]
            time_str = f" at {tmp*1000:.1f} ms"
        else:
            time_str=""
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        tc=colorlist[0%len(colorlist)]
        plt.title(r"$\psi_{m/n}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|\psi_{m/1}|$ [mWb]',fontsize=fontsize)
        thisamp=np.sqrt(np.square(self.bcmn)+np.square(self.bsmn))
        ax.plot(self.rhon,thisamp[0,4,:]*1000,color=colorlist[4],label=r"$m=-1$")
        ax.plot(self.rhon,thisamp[0,3,:]*1000,color=colorlist[1],label=r"$m=-2$")
        ax.plot(self.rhon,thisamp[0,2,:]*1000,color=colorlist[2],label=r"$m=-3$")
        ax.plot(self.rhon,thisamp[0,1,:]*1000,color=colorlist[3],label=r"$m=-4$")
        ax.axhline(0,ls='-',color='k')
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.1f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.legend(ncol=2,fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.show()

        vxb_inphase = np.zeros(self.bmn_amp.shape)
        vxb_outphase = np.zeros(self.bmn_amp.shape)
        for ikey in ['vpbeq','veqbp','vpbp']: #veqbeq is not calculated
            this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
            this_in,this_out=self.phase_shift2(this_cosmn,this_sinmn,phase_21)
            vxb_inphase+=this_in
            vxb_outphase+=this_out

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"$\nabla \times \left(\vec v \times \vec b\right)_{m/n}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'In phase component [Wb/s]',fontsize=fontsize)
        ax.plot(self.rhon,vxb_inphase[0,4,:],color=colorlist[4],label=r"$m=-1$")
        ax.plot(self.rhon,vxb_inphase[0,3,:],color=colorlist[1],label=r"$m=-2$")
        ax.plot(self.rhon,vxb_inphase[0,2,:],color=colorlist[2],label=r"$m=-3$")
        ax.plot(self.rhon,vxb_inphase[0,1,:],color=colorlist[3],label=r"$m=-4$")
        ax.axhline(0,ls='-',color='k')
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.1f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.legend(ncol=2,fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.show()

        reconnect_inphase = np.zeros(self.bmn_amp.shape)
        reconnect_outphase = np.zeros(self.bmn_amp.shape)
        for ikey in ['etajpert','divpie']:
            this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
            this_in,this_out=self.phase_shift2(this_cosmn,this_sinmn,phase_21)
            reconnect_inphase+=this_in
            reconnect_outphase+=this_out

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"-$\nabla \times \left(\eta J + \frac{1}{ne}\nabla \cdot \pi_e\right)_{m/n}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'In phase component [Wb/s]',fontsize=fontsize)
        ax.plot(self.rhon,reconnect_inphase[0,4,:],color=colorlist[4],label=r"$m=-1$")
        ax.plot(self.rhon,reconnect_inphase[0,3,:],color=colorlist[1],label=r"$m=-2$")
        ax.plot(self.rhon,reconnect_inphase[0,2,:],color=colorlist[2],label=r"$m=-3$")
        ax.plot(self.rhon,reconnect_inphase[0,1,:],color=colorlist[3],label=r"$m=-4$")
        ax.axhline(0,ls='-',color='k')
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.1f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.legend(ncol=2,fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.show()

        rhomin=0.5
        rhomax=0.8
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"-$\nabla \times \left(\eta J + \frac{1}{ne}\nabla \cdot \pi_e\right)_{m/n}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'In phase component [Wb/s]',fontsize=fontsize)
        range=np.where((self.rhon>rhomin) & (self.rhon<rhomax))[0]
        ax.plot(self.rhon[range],reconnect_inphase[0,4,range],color=colorlist[4],label=r"$m=-1$")
        ax.plot(self.rhon[range],reconnect_inphase[0,3,range],color=colorlist[1],label=r"$m=-2$")
        ax.plot(self.rhon[range],reconnect_inphase[0,2,range],color=colorlist[2],label=r"$m=-3$")
        ax.plot(self.rhon[range],reconnect_inphase[0,1,range],color=colorlist[3],label=r"$m=-4$")
        ax.axhline(0,ls='-',color='k')

        irho = self.get_rho_q(-2)
        qlbl = f"q = {-2:.1f}"
        tc=colorlist[1]
        ax.axvline(irho,ls=':',color=tc, label=qlbl)

        plt.legend(ncol=2,fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlim(rhomin,rhomax)
        plt.tight_layout()
        plt.show()



        eta_cosmn,eta_sinmn,eta_amp,eta_phase=self.interp_dict['etajpert']
        eta_in,eta_out=self.phase_shift2(eta_cosmn,eta_sinmn,phase_21)
        divpie_cosmn,divpie_sinmn,divpie_amp,divpie_phase=self.interp_dict['divpie']
        divpie_in,divpie_out=self.phase_shift2(divpie_cosmn,divpie_sinmn,phase_21)


        fig= plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"$\delta \psi_{-2/1}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'Re $\delta \psi$ [Wb/s]',fontsize=fontsize)
        irho = self.get_rho_q(-2)
        qlbl = f"q = {-2:.1f}"
        tc=colorlist[1]
        ax.axvline(irho,ls=':',color=tc,label="q=-2.0")
#        plt.text(irho*.985,-4,'q=-2',rotation=90,color=tc,fontsize=fontsize)
        range=np.where((self.rhon>rhomin) & (self.rhon<rhomax))[0]
        ax.plot(self.rhon[range],reconnect_inphase[0,3,range],label=r"$\delta \psi_{rec}$",color=colorlist[0])
        ax.plot(self.rhon[range],eta_in[0,3,range],label=r"$\delta \psi_{\eta} $",color=colorlist[2])
        ax.plot(self.rhon[range],divpie_in[0,3,range],label=r"$\delta \psi_{\pi_e}$",color=colorlist[3])
        ax.axhline(0,ls='-',color='k')



        plt.legend(fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlim(rhomin,rhomax)
        plt.tight_layout()
        plt.show()

        #rhomin=0.60
        #rhomax=0.85
        fig= plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"$\delta \psi_{-2/1}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'Re $\delta \psi$ [Wb/s]',fontsize=fontsize)
        irho = self.get_rho_q(-2)
        qlbl = f"q = {-2:.1f}"
        tc=colorlist[1]
        ax.axvline(irho,ls=':',color=tc,label="q=-2.0")
        range=np.where((self.rhon>rhomin) & (self.rhon<rhomax))[0]
        ax.plot(self.rhon[range],reconnect_inphase[0,3,range],label=r"$\delta \psi_{rec}$",color=colorlist[0])
        ax.plot(self.rhon[range],eta_in[0,3,range],label=r"$\delta \psi_{\eta} $",color=colorlist[2])
        ax.plot(self.rhon[range],divpie_in[0,3,range],label=r"$\delta \psi_{\pi_e}$",color=colorlist[3])
        ax.plot(self.rhon[range],vxb_inphase[0,3,range],label=r"$\delta \psi_{VB}$",color=colorlist[4])
        ax.axhline(0,ls='-',color='k')



        plt.legend(fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.tight_layout()
        plt.show()

        #rhomin=0.60
        #rhomax=0.8
        fig= plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"$\delta \psi_{-2/1}$"+time_str,fontsize=fontsize)
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'Re $\delta \psi$ [Wb/s]',fontsize=fontsize)
        irho = self.get_rho_q(-2)
        qlbl = f"q = {-2:.1f}"
        tc=colorlist[1]
        ax.axvline(irho,ls=':',color=tc,label="q=-2.0")
        range=np.where((self.rhon>rhomin) & (self.rhon<rhomax))[0]
        ax.plot(self.rhon[range],reconnect_inphase[0,3,range],label=r"$\delta \psi_{rec}$",color=colorlist[0])
        ax.plot(self.rhon[range],eta_in[0,3,range],label=r"$\delta \psi_{\eta} $",color=colorlist[2])
        ax.plot(self.rhon[range],divpie_in[0,3,range],label=r"$\delta \psi_{\pi_e}$",color=colorlist[3])
        ax.plot(self.rhon[range],vxb_inphase[0,3,range],label=r"$\delta \psi_{VB}$",color=colorlist[4])
        ax.axhline(0,ls='-',color='k')
        plt.legend(fontsize=fontsize,loc=2,frameon=True,handlelength=.75,handletextpad=0.3)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.ylim(-3,3)
        plt.tight_layout()
        plt.show()



        return None

    def plot_fsa_phase_old(self, key=None, **kwargs):
        qlist=[-1.5, -2,-3,-4]
        fontsize=18
        colorlist = list(mcolors.TABLEAU_COLORS)
        titledict={}
        titledict={}
        titledict['veqbeq']=r"$\nabla \times \left(V_0 \times B_0\right)_{2/1}$"
        titledict['vpbeq']=r"$\nabla \times \left(\tilde V \times B_0\right)_{2/1}$"
        titledict['veqbp']=r"$\nabla \times \left(V_0 \times \tilde B\right)_{2/1}$"
        titledict['vpbp']=r"$\nabla \times \left(\tilde V \times \tilde B\right)_{2/1}$"
        titledict['etajpert']=r"$\nabla \times \left(\eta \tilde J\right)_{2/1}$"
        titledict['divpie']=r"$\nabla \times \left(\frac{1}{ne}\nabla \cdot \pi_e\right)_{2/1}$"
        sumlist=['etajpert','divpie']
        phase_21=interp1d(self.rhon,self.bmn_phase[0,3,:], kind='linear')(self.get_rho_q(-2))
        binphase,boutphase=self.phase_shift2(self.bcmn[0,3,:],self.bsmn[0,3,:],phase_21)
        print(phase_21)

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        tc=colorlist[0%len(colorlist)]
        plt.title(r"$\psi_{2/1}$")
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|\psi_{2/1}|$ [mWb]',fontsize=fontsize)
        #ax.plot(self.rhon,self.bmn_amp[0,0,:]*1000,label='amp 0')
        #ax.plot(self.rhon,self.bmn_amp[0,1,:]*1000,label='amp 1')
        #ax.plot(self.rhon,self.bmn_amp[0,2,:]*1000,label='amp 2')
        ax.scatter(self.dvar_dict['rhon'],self.raw_bmn_amp[0,3,:]*1000,label='amp 3')
        plt.show()
        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        tc=colorlist[0%len(colorlist)]
        plt.title(r"$\psi_{2/1}$")
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|\psi_{2/1}|$ [mWb]',fontsize=fontsize)
        #ax.plot(self.rhon,self.bmn_amp[0,0,:]*1000,label='amp 0')
        #ax.plot(self.rhon,self.bmn_amp[0,1,:]*1000,label='amp 1')
        #ax.plot(self.rhon,self.bmn_amp[0,2,:]*1000,label='amp 2')
        ax.plot(self.rhon,self.bmn_amp[0,3,:]*1000,label='amp 3')
        #ax.plot(self.rhon,self.bmn_amp[0,4,:]*1000,label='amp 4')
        #ax.plot(self.rhon,self.bmn_amp[0,5,:]*1000,label='amp 5')
        #ax.plot(self.rhon,self.bmn_amp[0,6,:]*1000,label='amp 6')
        #ax.plot(self.rhon,self.bmn_amp[0,7,:]*1000,label='amp 7')
        #ax.plot(self.rhon,self.bmn_amp[0,8,:]*1000,label='amp 8')
        #ax.plot(self.rhon,self.bmn_amp[0,9,:]*1000,label='amp 9')

        ax.plot(self.rhon,binphase*1000,label='in phase')
        ax.plot(self.rhon,boutphase*1000,label='90 out of phase')
        plt.legend()
        ax.axhline(0,ls='-',color='k')
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.2f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.show()

        in_sum = np.zeros(self.bmn_amp.shape)
        out_sum = np.zeros(self.bmn_amp.shape)
        if key is None:
            for ikey in self.interp_dict:
                print(ikey)
                this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
                this_in,this_out=self.phase_shift2(this_cosmn,this_sinmn,phase_21)
                if ikey in sumlist:
                    in_sum+=this_in
                    out_sum+=this_out
            #amp
                fig = plt.figure(figsize=(10,8))
                ax=fig.add_subplot(111)
                plt.title(titledict.get(ikey,ikey))
                plt.plot(self.rhon,this_amp[0,0,:],label='amp')
                plt.plot(self.rhon,this_in[0,0,:],label='in phase')
                plt.plot(self.rhon,this_out[0,0,:], label='90 out of phase')
                plt.plot(self.rhon,this_amp[0,1,:],label='amp')
                plt.plot(self.rhon,this_in[0,1,:],label='in phase')
                plt.plot(self.rhon,this_out[0,1,:], label='90 out of phase')
                plt.plot(self.rhon,this_amp[0,2,:],label='amp')
                plt.plot(self.rhon,this_in[0,2,:],label='in phase')
                plt.plot(self.rhon,this_out[0,2,:], label='90 out of phase')
                plt.xlabel(r'$\rho_N$',fontsize=fontsize)
                plt.ylabel(r'$|F_{m/n}|$ [Wb/s]',fontsize=fontsize)
                ax.axhline(0,ls='-',color='k')
                plt.legend()
                for iq,qq in enumerate(qlist):
                    try:
                        irho = self.get_rho_q(qq)
                        qlbl = f"q = {qq:.2f}"
                        tc=colorlist[iq]
                        ax.axvline(irho,ls=':',color=tc, label=qlbl)
                    except:
                        print("q not found")
                plt.tight_layout()
                plt.show()

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        ikey='divpie'
        this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
        this_in,this_out=self.phase_shift2(this_cosmn,this_sinmn,phase_21)
        plt.title(titledict.get(ikey,ikey))
        plt.plot(self.rhon,this_amp[0,3,:],label='amp')
        plt.plot(self.rhon,this_in[0,3,:],label='in phase')
        plt.plot(self.rhon,this_out[0,3,:], label='90 out of phase')
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|F_{m/n}|$ [Wb/s]',fontsize=fontsize)
        ax.axhline(0,ls='-',color='k')
        plt.ylim(-50,50)
        plt.legend()
        for iq,qq in enumerate(qlist):
            try:
                irho = self.get_rho_q(qq)
                qlbl = f"q = {qq:.2f}"
                tc=colorlist[iq]
                ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.tight_layout()
        plt.show()

        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        ikey='etajpert'
        this_cosmn,this_sinmn,this_amp,this_phase=self.interp_dict[ikey]
        this_in,this_out=self.phase_shift2(this_cosmn,this_sinmn,phase_21)
        plt.title(titledict.get(ikey,ikey))
        plt.plot(self.rhon,this_amp[0,3,:],label='amp')
        plt.plot(self.rhon,this_in[0,3,:],label='in phase')
        plt.plot(self.rhon,this_out[0,3,:], label='90 out of phase')
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|F_{m/n}|$ [Wb/s]',fontsize=fontsize)
        ax.axhline(0,ls='-',color='k')
        plt.ylim(-50,50)
        plt.legend()
        for iq,qq in enumerate(qlist):
            try:
              irho = self.get_rho_q(qq)
              qlbl = f"q = {qq:.2f}"
              tc=colorlist[iq]
              ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.tight_layout()
        plt.show()


        fig = plt.figure(figsize=(10,8))
        ax=fig.add_subplot(111)
        plt.title(r"$\nabla \times \left(\eta \tilde j + \frac{1}{ne}\nabla \cdot \Pi_e\right)$" )
        plt.plot(self.rhon,np.sqrt(np.square(in_sum[0,3,:])+np.square(out_sum[0,3,:])),label='amp')
        plt.plot(self.rhon,in_sum[0,3,:],label='in phase')
        plt.plot(self.rhon,out_sum[0,3,:],label='90 out of phase')
        ax.axhline(0,ls='-',color='k')
        plt.xlabel(r'$\rho_N$',fontsize=fontsize)
        plt.ylabel(r'$|F_{m/n}|$ [Wb/s]',fontsize=fontsize)
        plt.ylim(-50,50)
        plt.legend()
        for iq,qq in enumerate(qlist):
            try:
              irho = self.get_rho_q(qq)
              qlbl = f"q = {qq:.2f}"
              tc=colorlist[iq]
              ax.axvline(irho,ls=':',color=tc, label=qlbl)
            except:
                print("q not found")
        plt.tight_layout()
        plt.show()

        return None
