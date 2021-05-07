#!/usr/bin/env python3
import eval_nimrod as eval
import fsa
import plot_nimrod as pn
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import h5py
import numpy as np
import pickle
import hcFields as hc
import nim_timer as timer
from scipy.interpolate import interp1d

class hcstepfsa:
    def __init__(self,dumpfile,nimrodin,lphi):
        self.dumpfile=dumpfile
        self.nimrodin=nimrodin
        self.fields=hc.hcfields(dumpfile,nimrodin)
        self.time=None
        self.step=None
        self.eval=None
        self.grid=None
        self.lphi=lphi
        self.fsa_power={}
        self.dvar_dict={}
        if lphi:
            self.set_phiplanes()

    def dump(self,file):
        pickle.dump(self.dumpfile,file)
        pickle.dump(self.nimrodin,file)
        pickle.dump(self.time,file)
        pickle.dump(self.step,file)
        pickle.dump(self.lphi,file)
        pickle.dump(self.nmodes,file)
        pickle.dump(self.nmax,file)
        pickle.dump(self.nphi,file)
        pickle.dump(self.phimax,file)
        pickle.dump(self.lphi,file)
        pickle.dump(self.fsa_power,file)
        pickle.dump(self.dvar_dict,file)


    def load(self,file):
        with open(file,'rb') as pickle_file:
            self.dumpfile=pickle.load(pickle_file)
            self.nimrodin=pickle.load(pickle_file)
            self.time=pickle.load(pickle_file)
            self.step=pickle.load(pickle_file)
            self.lphi=pickle.load(pickle_file)
            self.nmodes=pickle.load(pickle_file)
            self.nmax=pickle.load(pickle_file)
            self.nphi=pickle.load(pickle_file)
            self.phimax=pickle.load(pickle_file)
            self.lphi=pickle.load(pickle_file)
            self.fsa_power=pickle.load(pickle_file)
            self.dvar_dict=pickle.load(pickle_file)

    def interpolate_fsa(self,radial='rhon',npts=200, fsa = True):
        RadiusTitle = {
        'rhon':r"$\rho_N$",
        'phi':r"$\Phi$",
        'psin':r"$\psi_N$",
        'psi':r"$\psi$",
        }

        try:
            r_varable = self.dvar_dict[radial]
        except:
            print(f"Radial variable {radial} is not reconized")
            raise KeyError

        #use Radius Title to set default r_label
        self.r_label = RadiusTitle.get(radial,radial)
        print(self.r_label)
        rmin=np.min(r_varable)
        rmax=np.max(r_varable)

        self.r = np.linspace(rmin,rmax,1000,endpoint=True)
        self.dvar_interp = {}
        self.fsa_interp = {}
        if fsa:
            fac = 1.0
        else:
            fac = np.abs(self.dvar_dict['vprime'])

        for key, item in self.dvar_dict.items():
            self.dvar_interp[key] = interp1d(r_varable, item, kind='cubic')(self.r)
        for key, item in self.fsa_power.items():
            self.fsa_interp[key] = interp1d(r_varable, item*fac, kind='cubic')(self.r)

        self.r_of_q = interp1d(self.dvar_dict['q'], r_varable,
                               kind='cubic', fill_value="extrapolate")

    def get_r_of_q(self,qvalue):
        try:
            return self.r_of_q(qvalue)
        except:
            print(f"The safety factor {qvalue} is not it the domain")
            raise

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
        self.fields.grid = eval.EvalGrid(rzp3d)
        self.fields.grid.set_debug(debug)
        self.fields.grid.set_3d_symm()

    def set_2dgrid(self,start,stop,npts,debug=0):
        '''sets up a 2d grid, using non to determine the number of phi planes
           based on
        '''
        p1 = np.array([start[0], start[1], 0.0])
        p2 = np.array([stop[0], stop[1], 0.0])
        p3 = np.array([start[0], start[1], self.phimax])
        rzp2d = pn.PlotNimrod.grid_2d_gen(p1, p2, p3, npts,self.nphi)
        self.fields.grid=eval.EvalGrid(rzp2d)
        self.fields.grid.set_debug(debug)

    @staticmethod
    def calc_nmodes(lphi):
        nmodes=int(2**lphi/3)+1
        return nmodes

    @staticmethod
    def calc_nplanes(lphi):
        nplanes = 2**lphi
        return nplanes

    def set_phiplanes(self):
        self.nmodes = self.calc_nmodes(self.lphi)
        self.nmax = self.nmodes-1
        self.nphi = self.calc_nplanes(self.lphi)
        self.phimax = np.pi*2*(self.nphi-1)/self.nphi


    def set_fsagrid(self,r,z):
        rzp = np.array([np.broadcast_to(r,self.nphi),
                np.broadcast_to(z,self.nphi),
                np.linspace(0,self.phimax,self.nphi)])
        return rzp

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

    def find_comp_boundary(self,inpoint,outpoint,tol=1e-8,debug=False):
        index=0
        maxit=100
        ntest=self.fields.eval.eval_field('n',outpoint,dmode=0,eq=1)
        if ntest ==ntest:
            '''return if outpoint is in domain'''
            return outpoint
        fst=np.copy(inpoint)
        lst=np.copy(outpoint)
        if debug:
            print(fst,inpoint)
            print(lst,outpoint)
        dist=np.linalg.norm(lst-fst,2)
        if dist<tol:
            return fst
        while True:
            tst=(fst+lst)/2.0
            if debug:
                print(tst)
            ntest=self.fields.eval.eval_field('n',tst,dmode=0,eq=1)
            if ntest ==ntest:
                fst=tst
            else:
                lst=tst
            dist=np.linalg.norm(lst-fst,2)
            if debug:
                print(dist)
            if dist<tol:
                break
            index+=1
            if index>maxit:
                print("no convergence")
                break
        if debug:
            print(fst)
            print(lst)
        return fst

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

    @timer.timer_func
    def power_fsa_int(self,rzc,y,dy,evalnimrod,fargs):
        addpert=fargs.get("addpert",True)
        grid=self.set_fsagrid(rzc[0],rzc[1])
        self.fields.clean_up_fsa()
        self.fields.powerFlux(grid)
        self.fields.advectPowerFlux(grid)

#        dy[0-3] are used
        fac = fargs['sfac']
        minidx = fargs['nmin']
        idx = 4

        for key, item in self.fields.powerFluxDict.items():
            dy[idx:idx+self.nmodes]=(fac*item[minidx:minidx+self.nmodes] + 1.0)*dy[2]
            idx += self.nmodes
        for key, item in self.fields.advectDict.items():
            dy[idx:idx+self.nmodes]=(fac*item[minidx:minidx+self.nmodes] + 1.0)*dy[2]
            idx += self.nmodes
        return dy


    @timer.timer_func
    def calculate_power_fsa(self,rzo=None,rzx=None,nsurf=100,eqflag=0,fargs={},**kwargs):
        self.fields.set_method("powerflux")

#        nsurf=100
        dpow=kwargs.get('dpow',0.5)
        nmax = kwargs.get('nmax',5)
        rzo=np.array([1.76821,-0.0188439,0.0])
        oPoint=self.find_pf_null(self.fields.eval, rzo, flag=0)
        rzx=np.array([1.27,-1.14,0.0])
        xPoint=self.find_pf_null(self.fields.eval, rzx, flag=0)
        #find rmax
        rzmax=np.copy(oPoint)
        rzmax[0]=3.0
        rzmax=self.find_comp_boundary(oPoint,rzmax)
        #self.set_2dgrid(oPoint,rzmax,npts)
        print('rzmax',rzmax)


        self.fsa=True
        fargs['rtol']=1.e-8
        fargs['sfac']=1.e-2
        fargs['nmin']=1
        fargs['nmax']=5
        self.nmodes=nmax-fargs['nmin']+1
        self.nmax=nmax

        #to check nmax < nmodes
      #

#skip q integration
#        dvar, yvar, contours = fsa.FSA(self.fields.eval, rzo, self.dummy_fsa, 1, \
#          nsurf=nsurf,depvar='eta',dpow=0.5,rzx=rzx,flag=eqflag,normalize=True, \
#          fargs=fargs)

#        iend=-1
#        while np.isnan(yvar[:,iend]).any():
#            iend -= 1
#        iend += yvar.shape[1]+1
#        #unevaluated interpoate
#        self.qpsi=interp1d(dvar[2,:iend], dvar[7,:iend], kind='cubic',fill_value="extrapolate")


        #call induction at opoint to get number of fields
        self.fields.clean_up_fsa()
        self.fields.powerFlux(grid=self.set_fsagrid(*rzo[0:2]))
        self.fields.advectPowerFlux(grid=self.set_fsagrid(*rzo[0:2]))

        neq=0
        neq += len(self.fields.powerFluxDict) * self.nmodes
        neq += len(self.fields.advectDict) * self.nmodes
        self.fields.clean_up_fsa()

        dvar,yvar,contours = fsa.FSA(self.fields.eval, rzo, self.power_fsa_int, neq, \
          nsurf=nsurf,depvar='eta', dpow=dpow, flag=eqflag, normalize=False,\
          **fargs)
        iend=-1
        while np.isnan(yvar[:,iend]).any():
            iend -= 1
        iend += yvar.shape[1]+1

        #dvars
        self.dvar_dict={}
        self.dvar_dict['psin']=dvar[0,:iend]
        self.dvar_dict['rhon']=dvar[1,:iend]
        self.dvar_dict['psi']=dvar[2,:iend]
        self.dvar_dict['phi']=dvar[3,:iend] #"toroidal flux"
        self.dvar_dict['vprime']=dvar[6,:iend]
        self.dvar_dict['q']=dvar[7,:iend]

        idx=0
        self.fsa_power = {}
        for key in self.fields.powerFluxDict:
            self.fsa_power[key] = \
                (yvar[idx:idx+self.nmodes,:iend] - dvar[6,:iend]) / \
                (fargs['sfac'] * dvar[6,:iend])
            idx += self.nmodes
        for key in self.fields.advectDict:
            self.fsa_power[key] = \
                (yvar[idx:idx+self.nmodes,:iend] - dvar[6,:iend]) / \
                (fargs['sfac'] * dvar[6,:iend])
            idx += self.nmodes


        rhomin=np.min(dvar[1,:iend])
        rhomax=np.max(dvar[1,:iend])
        self.rhon = np.linspace(rhomin,rhomax,200,endpoint=True)

    #    self.raw_bcmn=bcmn
    #    self.raw_bsmn=bsmn
    #    self.raw_bmn_amp=bmn_amp
    #    self.raw_bmn_phase =bmn_phase

    #    self.interpolate_fsa()

        # neq=1+nterm*self.nfour*(4*self.mmax+1)
    #    fig = plt.figure(figsize=(10,8))
    #    ax=fig.add_subplot(111)
    #    conf=plt.plot(self.rhon,self.q)
    #    plt.show()

    #    fig =plt.figure(figsize=(10,8))
    #    ax=fig.add_subplot(111)
    #    conf=plt.plot(self.rhon,self.bmn_amp[0,2,:])
    #    plt.show()

    #    fig =plt.figure(figsize=(10,8))
    #    ax=fig.add_subplot(111)
    #    conf=plt.plot(self.rhon,self.bmn_phase[0,2,:])
    #    plt.show()

        return None

    def default_plot(self, fig_size = [8,6],qlist=None):

        DEFAULT_q = [-1.2, -1.5, -2, -3, -4]
        if not qlist:
            qlist = DEFAULT_q

        lin_power=self.fsa_interp['vxbeq'] + \
                  self.fsa_interp['jxbeq'] + \
                  self.fsa_interp['ngpp']
        qua_power=self.fsa_interp['vxbn0'] + self.fsa_interp['jxbn0']
        non_power=self.fsa_interp['vxbp'] + self.fsa_interp['jxbp']
        ohm_power=self.fsa_interp['etajp']
        visc_power=self.fsa_interp['divpip']
        neoi_power=self.fsa_interp['divPii']
        neoe_power=self.fsa_interp['divPie']
        poyn_power=self.fsa_interp['poynting']
        adv_power=self.fsa_interp['rhovdveq'] + \
                  self.fsa_interp['rhovdvn0'] + \
                  self.fsa_interp['rhovdvp']

        try:
            poyndis_power = self.fsa_interp['poyndis']
            poynlin_power = self.fsa_interp['poynlin']
            poynqln_power = self.fsa_interp['poynqln']
            poynnon_power = self.fsa_interp['poynnon']
        except:
            pass
        total_power = np.zeros_like(adv_power)
        skip_list = ['poyndis', 'poynlin', 'poynqln', 'poynnon']
        for key, item in self.fsa_interp.items():
            if key not in skip_list:
                total_power+=item

        #total_power -= self.fsa_interp['divpip']
        diss_power = ohm_power + neoi_power + neoe_power+ visc_power

        ## DEBUG:
        #print(np.max(self.dvar_dict['phi']))
        #print(np.min(self.dvar_dict['phi']))
        #fig, ax = plt.subplots(figsize=fig_size)
        #ax.plot(self.dvar_dict['phi'],marker='+')
        #plt.show()
        ##return

        y_label = r"$\Omega_n$ [MW/m^3]"
        y_label = r"$\Omega_n$ [MW/m]"
        xlim = [0,0.85]
        max_idx = -1
        while np.max(self.r[:max_idx])>0.9:
            max_idx -= 1

        nmin = self.nmax+1 -self.nmodes
        nlist = range(nmin,self.nmax+1)
        for idx,nn in enumerate(nlist):
            fig, ax = plt.subplots(figsize=fig_size)

            ax.plot(self.r[:max_idx], total_power[idx,:max_idx]/10**6, alpha=0.7, label="Total",color='tab:brown')
            ax.plot(self.r[:max_idx], lin_power[idx,:max_idx]/10**6, alpha=0.7, label="Lin",color='tab:orange')
            ax.plot(self.r[:max_idx], qua_power[idx,:max_idx]/10**6, alpha=0.7, label="QL",color ='tab:green')
            ax.plot(self.r[:max_idx], non_power[idx,:max_idx]/10**6, alpha=0.7, label="NL",color='tab:red')
            ax.plot(self.r[:max_idx], diss_power[idx,:max_idx]/10**6, alpha=0.7, label="Diss",color='tab:purple')
            ax.plot(self.r[:max_idx], poyn_power[idx,:max_idx]/10**6, alpha=0.7, label="PF",color='tab:blue')


            title = f"n={nn} power at {self.time*1000:.2f}ms"

            for q in qlist:
                print(q)
                try:
                    rq = self.get_r_of_q(q)
                    print(rq)
                    ax.axvline(rq,ls=':')
                except:
                    print("could not find q")
                    pass
            plt.legend(ncol=2, loc='upper left', fontsize = 18, frameon=True, framealpha=0.8,handlelength=1)

            ax.set(xlabel=self.r_label, ylabel=y_label, title=title, xlim=xlim)

            plt.axhline(0,color='k')
            #ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
            #            useOffset=None, useLocale=None, useMathText=True)

            plt.tight_layout()
            plt.show()
### plot 2
        for idx,nn in enumerate(nlist):
            fig, ax = plt.subplots(figsize=fig_size)

            ax.plot(self.r[:max_idx], total_power[idx,:max_idx]/10**6, alpha=0.7, label="Total",color='tab:brown')
            #ax.plot(self.r[:max_idx], (lin_power[idx,:max_idx]+poynlin_power[idx,:max_idx])/10**6, alpha=0.7, label="Lin+LPF",color='tab:orange')
            #ax.plot(self.r[:max_idx], (qua_power[idx,:max_idx]+poynqln_power[idx,:max_idx])/10**6, alpha=0.7, label="QL+QPF",color ='tab:green')
            ax.plot(self.r[:max_idx], (non_power[idx,:max_idx])/10**6, alpha=0.7, label="NL",color='tab:red')
            #ax.plot(self.r[:max_idx], (diss_power[idx,:max_idx]+poyndis_power[idx,:max_idx])/10**6, alpha=0.7, label="Diss+DPF",color='tab:purple')
            #ax.plot(self.r[:max_idx], poyn_power[idx,:max_idx]/10**6, alpha=0.7, label="PF",color='tab:blue')


            title = f"n={nn} power at {self.time*1000:.2f}ms"
            qlist2 = [-1.2,-4/3, -1.5,-5/3, -2, -3, -4]
            c2list = ['tab:blue','tab:red','tab:blue','tab:red','tab:blue','tab:blue','tab:blue']
            idx =0
            for q in qlist2:
                print(q)
                try:
                    rq = self.get_r_of_q(q)
                    print(rq)
                    ax.axvline(rq,ls=':',color = c2list[idx])
                    idx+=1
                except:
                    print("could not find q")
                    pass
            plt.legend(ncol=1, loc='upper left', fontsize = 18, frameon=True, framealpha=0.8,handlelength=1)
            ax.set(xlabel=self.r_label, ylabel=y_label, title=title, xlim=xlim)


            #ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
            #            useOffset=None, useLocale=None, useMathText=True)
            plt.axhline(0,color='k')
            plt.tight_layout()
            plt.show()
#plot 2b
        for idx,nn in enumerate(nlist):
            fig, ax = plt.subplots(figsize=fig_size)

            ax.plot(self.r[:max_idx], total_power[idx,:max_idx]/10**6, alpha=0.7, label="Total",color='tab:brown')
            ax.plot(self.r[:max_idx], (lin_power[idx,:max_idx]+poyn_power[idx,:max_idx])/10**6, alpha=0.7, label="Lin+PF",color='tab:orange')
            ax.plot(self.r[:max_idx], (qua_power[idx,:max_idx])/10**6, alpha=0.7, label="QL",color ='tab:green')
            ax.plot(self.r[:max_idx], (non_power[idx,:max_idx])/10**6, alpha=0.7, label="NL",color='tab:red')
            ax.plot(self.r[:max_idx], (diss_power[idx,:max_idx])/10**6, alpha=0.7, label="Diss",color='tab:purple')
#            ax.plot(self.r[:max_idx], poyn_power[idx,:max_idx]/10**6, alpha=0.7, label="PF",color='tab:blue')


            title = f"n={nn} power at {self.time*1000:.2f}ms"
            qlist2 = [-1.2,-4/3, -1.5,-5/3, -2, -3, -4]
            c2list = ['tab:blue','tab:red','tab:blue','tab:red','tab:blue','tab:blue','tab:blue']
            idx =0
            for q in qlist2:
                print(q)
                try:
                    rq = self.get_r_of_q(q)
                    print(rq)
                    ax.axvline(rq,ls=':',color = c2list[idx])
                    idx+=1
                except:
                    print("could not find q")
                    pass
            plt.legend(ncol=1, loc='upper left', fontsize = 18, frameon=True, framealpha=0.8,handlelength=1)

            ax.set(xlabel=self.r_label, ylabel=y_label, title=title, xlim=xlim)


            #ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
            #            useOffset=None, useLocale=None, useMathText=True)
            plt.axhline(0,color='k')
            plt.tight_layout()
            plt.show()
#plot 3
        for idx,nn in enumerate(nlist):
            fig, ax = plt.subplots(figsize=fig_size)

            ax.plot(self.r[:max_idx], (lin_power[idx,:max_idx]+poynlin_power[idx,:max_idx])/10**6, alpha=0.7, label="lin+poy")
            ax.plot(self.r[:max_idx], qua_power[idx,:max_idx]/10**6, alpha=0.7, label="qlin")
            ax.plot(self.r[:max_idx], (non_power[idx,:max_idx])/10**6, alpha=0.7, label="non")
            ax.plot(self.r[:max_idx], (non_power[idx,:max_idx]+poynnon_power[idx,:max_idx])/10**6, alpha=0.7, label="non+poy")
            ax.plot(self.r[:max_idx], diss_power[idx,:max_idx]/10**6, alpha=0.7, label="diss")
            ax.plot(self.r[:max_idx], total_power[idx,:max_idx]/10**6, alpha=0.7, label="tot")


            title = f"n={nn} power at {self.time*1000:.2f}ms"

            for q in qlist:
                print(q)
                try:
                    rq = self.get_r_of_q(q)
                    print(rq)
                    ax.axvline(rq,ls=':')
                except:
                    print("could not find q")
                    pass
            ax.legend(loc='best',ncol=2)

            ax.set(xlabel=self.r_label, ylabel=y_label, xlim=xlim, title=title)


            #ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
            #            useOffset=None, useLocale=None, useMathText=True)

            plt.tight_layout()
            plt.show()

        for idx,nn in enumerate(nlist):
            fig, ax = plt.subplots(figsize=fig_size)

            ax.plot(self.r[:max_idx], visc_power[idx,:max_idx]/10**6, alpha=0.7, label="visc")
            ax.plot(self.r[:max_idx], neoe_power[idx,:max_idx]/10**6, alpha=0.7, label="neoe")
            ax.plot(self.r[:max_idx], neoi_power[idx,:max_idx]/10**6, alpha=0.7, label="neoi")
            ax.plot(self.r[:max_idx], ohm_power[idx,:max_idx]/10**6, alpha=0.7, label="ohm")
            #ax.plot(self.r[:max_idx], adv_power[idx,:max_idx]/10**6, alpha=0.7, label="adv")


            title = f"n={nn} power at {self.time*1000:.2f}ms"

            for q in qlist:
                print(q)
                try:
                    rq = self.get_r_of_q(q)
                    print(rq)
                    ax.axvline(rq,ls=':')
                except:
                    print("could not find q")
                    pass
            ax.legend(loc='best',ncol=2)

            ax.set(xlabel=self.r_label, ylabel=y_label, xlim=xlim, title=title)


            #ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
            #            useOffset=None, useLocale=None, useMathText=True)

            plt.tight_layout()
            plt.show()


        return
    def clean_up(self):
        ''' Clean up fields to reduce memory foot print'''
        self.fields.clean_up()
        self.grid=None
