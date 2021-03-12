#!/usr/bin/env python3

import f90nml
import eval_nimrod as eval
import field_class as fc
import numpy as np
from scipy.fft import fft, ifft, rfft
import plot_nimrod as pn
import matplotlib.pyplot as plt
import time
import nim_timer as timer

class hcfields:
    ''' hc fields is a class for reading/storing nimrod data on a mesh
        it also calculates the different terms in nimrods field advances
        this seperates the evaulation of these fields from fsa interation
        or plotting on various meshs
    '''

    def __init__(self,dumpfile,nimrodin):
        self.dumpfile=dumpfile
        self.nimrodin=nimrodin
        self.eval=None
        self.grid=None
        self.fielddict={}
        self.ndict={}
        self.edict={}
        self.energyDict={}
        self.powerFluxDict={}
        self.evalb_timer=0.0
        self.fft_timer=0.0
        self.ndmode=0
        self.neq=False
        self.ne0=False
        self.npert=False
        self.pdmode=1
        self.peq=False
        self.p0=False
        self.ppert=False
        self.vdmode=0
        self.veq=False
        self.ve0=False
        self.vpert=False
        self.bdmode=0
        self.beq=False
        self.be0=False
        self.bpert=False
        self.jdmode=1
        self.jeq=False
        self.je0=False
        self.jpert=False
        self.diff_dmode=0

        if self.nimrodin is not None:
          self.nml=f90nml.read(self.nimrodin)
          self.grid_nml=self.nml.get('grid_input',{})
          self.equil_nml=self.nml.get('equil_input',{})
          self.const_nml=self.nml.get('const_input',{})
          self.physics_nml=self.nml.get('physics_input',{})
          self.closure_nml=self.nml.get('closure_input',{})
          self.solver_nml=self.nml.get('solver_input',{})
          self.output_nml=self.nml.get('output_input',{})
          self.set_physical_constants()
          self.set_evalnimrod()
        return None

    def set_physical_constants(self):
        '''
        Read namelists and set physical constants
        '''
        self.mu0=self.const_nml.get('mu0_input',np.pi*4.0e-7)
        self.me=self.const_nml.get('me_input',9.1093898e-31)
        self.mi=self.const_nml.get('mi_input',3.3435860e-27)
        self.zeff=self.const_nml.get('zeff_input',1.0)
        self.mtot=self.me + self.mi/self.zeff
        self.qe=self.const_nml.get('chrg_input',1.60217733e-19)
        self.gamma=self.const_nml.get('gam_input',5./3.)
        self.kblz=self.const_nml.get('kblz_input',1.60217733e-19)
        self.clight=self.const_nml.get('c_input',2.99792458e8)
        return None

    def set_evalnimrod(self):
        '''
        Set Eval Nimrod
        '''
        if self.eval is None:
            self.eval=eval.EvalNimrod(self.dumpfile,fieldlist='nvbtjpd')
        return None

    def get_gridrzp(self,grid):
        '''
        returns grid and rzp, and grid for a given grid
        if grid=None, try self.grid
        elif grid is an EvalGrid insatnace
        else np grid
        '''
        if grid is None:
          if self.grid is None:
            print("ntm_fields_grid is not set")
            raise ValueError
          else:
            grid=self.grid
            rzp=self.grid.rzp
        elif isinstance(grid,eval.EvalGrid):
          rzp=grid.rzp
        else:
          rzp=grid
        return grid,rzp

    @timer.timer_func
    def fft(self,pert,axis=-1,type=None):
        ''' NIMROD stores it's field data as f(phi) = sum_{-n}^n f_n exp(inphi)
            This implies that the normalization 1/N should be done in the transform
            from physical space to fourier space
            This is the one option that scipi.fft does not support, but why?
        '''
        #fpert = fft(pert.data,axis=axis,norm=None)/pert.data.shape[axis]
        #try rfft to save space (computation)
        fpert = rfft(pert.data,axis=axis,norm=None)/pert.data.shape[axis]
        if type=='s':
            fpert = fpert[0,...]
        elif type=='v':
            fpert = fpert[0:3,...]
        return fpert

    def set_method(self,method):
        if method == "induction":
            self.ndmode=1
            self.neq=True
            self.npert=False
            self.vdmode=1
            self.veq=True
            self.vpert=True
            self.beq=True
            self.bpert=True
            self.jeq=True
            self.jpert=True
            self.peq=True
            self.p0=True
            self.ppert=True
            self.diff_dmode=1
        elif method == 'powerflux':
            self.ndmode=1 #check
            self.neq=True
            self.ne0=True
            self.npert=False
            self.vdmode=2
            self.veq=True
            self.ve0=True
            self.vpert=True
            self.bdmode=1
            self.beq=True
            self.be0=True
            self.bpert=True
            self.jdmode=1
            self.jeq=True
            self.je0=True
            self.jpert=True
            self.pdmode=1
            self.peq=True
            self.p0=True
            self.ppert=True
            self.diff_dmode=1

    @timer.timer_func
    def eval_symm(self,fname,rzp,dmode,eq):
        if eq not in [1,3]:
            print("eval_symm only works for eq=1 or 3")
            raise ValueError
        newrzp=rzp[:,0]
        field=self.eval.eval_field(fname,newrzp,dmode=dmode,eq=eq)
        field=field.reshape(*field.shape,1)
        field=np.broadcast_to(field,(field.shape[0],rzp.shape[1]))
        return field

    @timer.timer_func
    def eval_n(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if self.neq or self.ne0:
            if len(rzp.shape)==2:
                field=self.eval_symm('n',rzp,dmode=self.ndmode,eq=1)
                self.fielddict['neq']=fc.Scalar(field,rzp,self.ndmode,True)
            else:
                field=self.eval.eval_field('n', grid, dmode=self.ndmode, eq=1)
                self.fielddict['neq']=fc.Scalar(field,rzp,self.ndmode,True)
        if self.ne0:
            if len(rzp.shape)==2:
                field2=self.eval_symm('n',rzp,dmode=self.ndmode,eq=3)
                field2=field2-field #sutract eq
                self.fielddict['ne0']=fc.Scalar(field2,rzp,self.ndmode,True)
            else:
                field2=self.eval.eval_field('n', grid, dmode=self.ndmode, eq=3)
                field2=field2-field #subtract eq
                self.fielddict['ne0']=fc.Scalar(field2,rzp,self.ndmode,True)
        if self.npert:
            field=self.eval.eval_field('n', grid, dmode=self.ndmode, eq=0)
            if self.ne0:
                field=field-field2 #remove n=0 perturbation
            self.fielddict['npert']=fc.Scalar(field,rzp,self.ndmode,True)
            if fft:
                self.fielddict['nfour']=self.fft(field,type='s')
        return None


    @timer.timer_func
    def eval_v(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if self.veq or self.ve0:
            if len(rzp.shape)==2:
                field=self.eval_symm('v',rzp,dmode=self.vdmode,eq=1)
                self.fielddict['veq']=fc.Vector(field,rzp,self.vdmode,True)
            else:
                field=self.eval.eval_field('v', grid, dmode=self.vdmode,eq=1)
                self.fielddict['veq']=fc.Vector(field,rzp,self.vdmode,True)
        if self.ve0:
            if len(rzp.shape)==2:
                field2=self.eval_symm('v',rzp,dmode=self.vdmode,eq=3)
                field2=field2-field #sutract eq
                self.fielddict['ve0']=fc.Vector(field2,rzp,self.vdmode,True)
            else:
                field2=self.eval.eval_field('v', grid, dmode=self.vdmode, eq=3)
                field2=field2-field #subtract eq
                self.fielddict['ve0']=fc.Vector(field2,rzp,self.vdmode,True)
        if self.vpert:
            field=self.eval.eval_field('v', grid, dmode=self.vdmode, eq=0)
            if self.ve0:
                field=field-field2 #subtract n=0 perturbation
            self.fielddict['vpert']=fc.Vector(field,rzp,self.vdmode,True)
            if fft:
                self.fielddict['vfour']=self.fft(field,type='v')
        return None

    @timer.timer_func
    def eval_b(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if self.beq or self.be0:
            if len(rzp.shape)==2:
                field=self.eval_symm('b',rzp,dmode=self.bdmode,eq=1)
                self.fielddict['beq']=fc.Vector(field,rzp,self.bdmode,True)
            else:
                field=self.eval.eval_field('b',grid,dmode=self.bdmode,eq=1)
                self.fielddict['beq']=fc.Vector(field,rzp,self.bdmode,True)
        if self.be0:
            if len(rzp.shape)==2:
                field2=self.eval_symm('b',rzp,dmode=self.bdmode,eq=3)
                field2=field2-field #sutract eq
                self.fielddict['be0']=fc.Vector(field2,rzp,self.bdmode,True)
            else:
                field2=self.eval.eval_field('b', grid, dmode=self.bdmode, eq=3)
                field2=field2-field #subtract eq
                self.fielddict['be0']=fc.Vector(field2,rzp,self.bdmode,True)
        if self.bpert:
            field=self.eval.eval_field('b',grid,dmode=self.bdmode, eq=0)
            if self.be0:
                field=field-field2 #subtract n=0 perturbation
            self.fielddict['bpert']=fc.Vector(field,rzp,self.bdmode,True)
            if fft:
              self.fielddict['bfour']=self.fft(field,type='v')
        return None

    @timer.timer_func
    def eval_fsa_beq2(self,grid=None):
        grid,rzp=self.get_gridrzp(grid)
        self.dump_fsa_beq2 =self.output_nml.get('dump_fsa_beq2',False)
        if self.dump_fsa_beq2:
            if len(rzp.shape)==2:
                field=self.eval_symm('fsa_beq2', grid, dmode=1, eq=1)
                self.fielddict['fsa_beq2']=fc.Scalar(field,rzp,1,True)
            else:
                field=self.eval.eval_field('fsa_beq2', grid, dmode=1, eq=1)
                self.fielddict['fsa_beq2']=fc.Scalar(field,rzp,1,True)
        else:
            self.fielddict['fsa_beq2']=0.0
        return None

    @timer.timer_func
    def eval_j(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if self.jeq or self.je0:
            if len(rzp.shape)==2:
                field=self.eval_symm('j',rzp,dmode=self.jdmode,eq=1)
                self.fielddict['jeq']=fc.Vector(field,rzp,self.jdmode,True)
            else:
                field=self.eval.eval_field('j', grid, dmode=self.jdmode, eq=1)
                self.fielddict['jeq']=fc.Vector(field,rzp,self.jdmode,True)
        if self.je0:
            if len(rzp.shape)==2:
                field2=self.eval_symm('j',rzp,dmode=self.jdmode,eq=3)
                field2=field2-field #sutract eq
                self.fielddict['je0']=fc.Vector(field2,rzp,self.jdmode,True)
            else:
                field2=self.eval.eval_field('j', grid, dmode=self.jdmode, eq=3)
                field2=field2-field #subtract eq
                self.fielddict['je0']=fc.Vector(field2,rzp,self.jdmode,True)
        if self.jpert:
            field=self.eval.eval_field('j', grid, dmode=self.jdmode, eq=0)
            if self.je0:
                field=field-field2
            self.fielddict['jpert']=fc.Vector(field,rzp,self.jdmode,True)
            if fft:
                self.fielddict['jfour']=self.fft(field,type='v')
        return None

    @timer.timer_func
    def eval_neo_mask(self,grid=None):
        grid,rzp=self.get_gridrzp(grid)
        r0=self.closure_nml.get('neo_axis_r',0)
        z0=self.closure_nml.get('neo_axis_z',0)
        rbump=self.closure_nml.get('neo_bump_r0',1.0)
        shape=list(rzp.shape)
        shape[0]=4
        fval=np.zeros(shape)
        fval[0]=1.0
        r2=(np.power(rzp[0]-r0,2)+np.power(rzp[1]-z0,2))/rbump**2
        if isinstance(r2,np.float):
            if r2<1.0:
                #check
                bump=np.exp(1-1./(1.-r2))
                dbumpdr2=-r2/(1-r2)**3
                dr2dx=2*(rzp[0]-r0)/rbump**2
                dr2dz=2*(rzp[1]-r0)/rbump**2
                d2r2dxx=2/rbump**2
                d2r2dyy=2/rbump**2
                d2r2dxy=0.0
                #print(type(bump),type(dbumpdr2),type(dr2dx))
                #print(dbumpdr2.shape)
                fval[0]=1.0-bump
                fval[1]=-bump*dbumpdr2*dr2dx
                fval[2]=-bump*dbumpdr2*dr2dz
        else:
            result=np.where(r2<1.0)
            for indicies in zip(*result):
                #check
                bump=np.exp(1-1./(1.-r2[indicies]))
                dbumpdr2=-r2[indicies]/(1-r2[indicies])**3
                dr2dx=2*(rzp[(0,)+indicies]-r0)/rbump**2
                dr2dz=2*(rzp[(1,)+indicies]-r0)/rbump**2
                d2r2dxx=2/rbump**2
                d2r2dyy=2/rbump**2
                d2r2dxy=0.0
                #print(type(bump),type(dbumpdr2),type(dr2dx))
                #print(dbumpdr2.shape)
                fval[(0,)+indicies]=1.0-bump
                fval[(1,)+indicies]=-bump*dbumpdr2*dr2dx
                fval[(2,)+indicies]=-bump*dbumpdr2*dr2dz
        return fc.Scalar(fval,rzp,1,True)

    @timer.timer_func
    def eval_p(self,grid=None,fft=False):
          grid,rzp=self.get_gridrzp(grid)
          if self.peq or self.p0:
              if len(rzp.shape)==2:
                  field=self.eval_symm('p',rzp,dmode=self.pdmode,eq=1)
                  self.fielddict['peq']=fc.Scalar(field,rzp,self.pdmode,True)
              else:
                  field=self.eval.eval_field('p',grid,dmode=self.pdmode,eq=1)
                  self.fielddict['peq']=fc.Scalar(field,rzp,self.pdmode,True)
          if self.p0:
              if len(rzp.shape)==2:
                  field2=self.eval_symm('p',rzp,dmode=self.pdmode,eq=3)
                  field2=field2-field #sutract eq
                  self.fielddict['p0']=fc.Scalar(field2,rzp,self.pdmode,True)
              else:
                  field2=self.eval.eval_field('p', grid, dmode=self.pdmode, eq=3)
                  field2=field2-field #subtract eq
                  self.fielddict['p0']=fc.Scalar(field2,rzp,self.pdmode,True)
          if self.ppert:
              field=self.eval.eval_field('p',grid,dmode=self.pdmode,eq=0)
              if self.p0:
                  field=field-field2 #remove n=0 perturbation
              self.fielddict['ppert']=fc.Scalar(field,rzp,self.pdmode,True)
              if fft:
                  self.fielddict['pfour']=self.fft(field,type='s')
          return None


    @timer.timer_func
    def eval_diff(self,grid=None):
        ''' Get the diff shape scalars from eval nimrod
            Some extralogic is needed to pull apart the different
            diff shape scalar fields from the output of eval_field
            note when ds_nqty>1
            elecd_diff is 0
            iso_visc_diff is 1
        '''
        diff_dmode=self.diff_dmode
        grid,rzp=self.get_gridrzp(grid)
        self.ds_nqty = self.equil_nml.get('ds_nqty',1)
        if len(rzp.shape)==2:
            field=self.eval_symm('d',rzp,dmode=diff_dmode,eq=1)
        else:
            field=self.eval.eval_field('d',grid,dmode=diff_dmode, eq=1)
        diff_shape=[]
        ishape=list(field.shape)
        ishape[0]=ishape[0]//self.ds_nqty
        ifield= np.zeros(ishape)
        for ii in range(self.ds_nqty):
            ifield[0] = field[ii]
            if diff_dmode>0:
                start=self.ds_nqty+ii*3
                end=start+3
                ifield[1:4]=field[start:end]
            if diff_dmode>1:
                start=self.ds_nqty*4+ii*6
                end=start+6
                ifield[4:10]=field[start:end]
            diff_shape.append(fc.Scalar(ifield,rzp,diff_dmode,True))
        self.fielddict['diff_shape']=diff_shape
        return None

    @timer.timer_func
    def energyDensity(self,grid=None):
        grid,rzp=self.get_gridrzp(grid)
        if 'veq' not in self.fielddict:
            self.eval_v(grid,fft=True)
        if 'beq' not in self.fielddict:
            self.eval_b(grid,fft=True)
        if 'neq' not in self.fielddict:
            self.eval_n(grid,fft=False)

        bTempP=self.fielddict['bfour']
        fac=1./(self.mu0)
        bpert=fac*np.sum(bTempP[0:3]*np.conjugate(bTempP[0:3]),axis=0,where=~np.isnan(bTempP[0:3]))
        bpert[...,0]=bpert[...,0]/2.0
        self.energyDict['bpert']=bpert
        fac=fac/2
        bTempE=self.fielddict['beq'].data[...,0]
        beq=fac*np.sum(bTempE[0:3]*np.conjugate(bTempE[0:3]),axis=0,where=~np.isnan(bTempE[0:3]))
        self.energyDict['beq']=beq
        bTempZ=bTempE+self.fielddict['be0'].data[...,0]
        be0=fac*np.sum(bTempZ[0:3]*np.conjugate(bTempZ[0:3]),axis=0,where=~np.isnan(bTempZ[0:3]))
        self.energyDict['be0']=be0

        vTempP=self.fielddict['vfour']
        nfour=vTempP.shape[-1]
        vTempE=self.fielddict['veq'].data[...,:nfour]
        vTempZ=vTempE[0:3,...,0]+self.fielddict['vfour'][0:3,...,0]
        nTempE=self.fielddict['neq'].data[...,:nfour]

        fac=self.mi
        vpert=fac*np.sum(vTempP[0:3]*np.conjugate(vTempP[0:3]),axis=0,where=~np.isnan(bTempP[0:3]))
        vpert=np.multiply(nTempE,vpert,out=vpert,where=~np.isnan(nTempE))
        vpert[...,0]=vpert[...,0]/2.0
        self.energyDict['vpert']=vpert

        veq=fac*np.sum(vTempE[0:3,...,0]*np.conjugate(vTempE[0:3,...,0]),axis=0,where=~np.isnan(vTempE[0:3,...,0]))
        veq=np.multiply(nTempE[...,0],veq,out=veq,where=~np.isnan(nTempE[...,0]))
        veq=veq/2.0
        self.energyDict['veq']=veq

        return

    @timer.timer_func
    def calculate_viscositiy(self,grid=None):
        ''' Calculate the divergence of the ion viscous stress tensor '''
        grid,rzp=self.get_gridrzp(grid)
        if 'veq' not in self.fielddict:
            self.eval_v(grid,fft=True)
        iso_visc = self.physics_nml.get('iso_visc',0.0)
        ndens = self.equil_nml.get('ndens',1.e20)
        mtot = self.mtot
        if iso_visc > 0:
            visc_fac = -1.0 * mtot * ndens * iso_visc * \
                       self.fielddict['diff_shape'][1]
            pi_tensor = self.fielddict['vpert'].grad() + \
                       self.fielddict['vpert'].grad().transpose() - \
                       2./3. * fc.eye(rzp) * self.fielddict['vpert'].div()
            pi_tensor *= visc_fac
            return pi_tensor.div()
        else:
            return None

    @timer.timer_func
    def calculate_neo_div_pi(self,grid=None):
        'calculates div Pi_i and 1/ne div Pi_e'
        grid,rzp=self.get_gridrzp(grid)
        if 'fsa_beq2' not in self.fielddict:
            self.eval_fsa_beq2(grid)
        mu_i=self.closure_nml.get('mu_i',0.0)
        mu_e=self.closure_nml.get('mu_e',0.0)
        beq=self.fielddict['beq']
        neq=self.fielddict['neq']
        fsa_beq2=self.fielddict['fsa_beq2']
        jpert=self.fielddict['jpert']
        vpert=self.fielddict['vpert']
        ephi=[0,0,1]
        etor=fc.basis_vector('p',rzp,torgeom=True)
        bep=beq-beq.dot(etor)*etor

        neo_mask=self.eval_neo_mask(grid)
        work_vec = neo_mask * fsa_beq2 /(beq.dot(bep)**2+1.0e-8) * bep

        #ecoef=self.mi*mu_i*neo_mask
        ifac = self.mi * mu_i
        efac =-self.me * mu_e / (self.qe**2)
        div_pi_i = ifac * neq * vpert.dot(bep) * work_vec
        div_pi_e = efac / neq * jpert.dot(bep) * work_vec
        self.fielddict['ndivPiipert'] = div_pi_i
        self.fielddict['ndivPiepert'] = div_pi_e


    @timer.timer_func
    def calculate_E(self,grid=None):
        ''' Calculate the electric field and the Curl of E for the Poynting Flux
        '''
        grid,rzp=self.get_gridrzp(grid)

        v_vec = self.fielddict['veq'] + self.fielddict['ve0'] + \
                self.fielddict['vpert']
        b_vec = self.fielddict['beq'] + self.fielddict['be0'] + \
                self.fielddict['bpert']

        j_vec = self.fielddict['jpert'] #   ONLY N!= IS needed

        div_pi_term = self.fielddict['ndivPiepert']

        Efield = - v_vec.cross(b_vec) \
                 + self.mu0 * self.elecd * j_vec - div_pi_term

        self.fielddict['epert'] = Efield

        return None

    @timer.timer_func
    def dotPlusCc(self,field1,field2):
        #first check dims
        if field1.shape!=field2.shape:
            print("array dimensions dont match")
            raise ValueError
        fshape=field1.shape[1:]
        result=np.zeros(fshape,dtype=complex)
        for ii in range(fshape[-1]):
            result[...,ii]=np.sum(field1[...,ii]*np.conj(field2[...,ii]) ,axis=0)
            if ii != 0:
                result[...,ii]+=np.conj(result[...,ii])
        return np.nan_to_num(np.real(result))

    @timer.timer_func
    def powerFlux(self,grid=None):
        grid,rzp=self.get_gridrzp(grid)
        if 'veq' not in self.fielddict:
            self.eval_v(grid,fft=True)
        if 'beq' not in self.fielddict:
            self.eval_b(grid,fft=True)
        if 'neq' not in self.fielddict:
            self.eval_n(grid,fft=False)
        if 'jeq' not in self.fielddict:
            self.eval_j(grid,fft=True)
        if 'peq' not in self.fielddict:
            self.eval_p(grid,fft=False)
        if 'diff_shape' not in self.fielddict:
            self.eval_diff(grid)

        ohmslaw=self.physics_nml.get('ohms','mhd')
        neoe_flag = self.closure_nml.get('neoe_flag',None)
        mu_e=self.closure_nml.get('mu_e',0.0)
        try:
            self.elecd=self.physics_nml.get('elecd',0.0)
            if isinstance(self.elecd,(np.ndarray,list)):
                self.elecd=self.elecd[0]
        except:
            print("Can't read elecd from nimrod.in")
            raise KeyError
        self.calculate_neo_div_pi()
        self.calculate_E()

        jeq_cross_bpert = self.fielddict['jeq'].cross(self.fielddict['bpert'],dmod=0)
        jpert_cross_beq = self.fielddict['jpert'].cross(self.fielddict['beq'],dmod=0)
        j0_cross_bpert = self.fielddict['je0'].cross(self.fielddict['bpert'],dmod=0)
        jpert_cross_b0 = self.fielddict['jpert'].cross(self.fielddict['be0'],dmod=0)
        jpert_cross_bpert =self.fielddict['jpert'].cross(self.fielddict['bpert'],dmod=0)

        j_cross_b_eq = self.fft(jeq_cross_bpert+jpert_cross_beq)
        j_cross_b_n0 = self.fft(j0_cross_bpert+jpert_cross_b0)
        j_cross_b_pert = self.fft(jpert_cross_bpert)

        veq_cross_bpert = self.fielddict['veq'].cross(self.fielddict['bpert'],dmod=0)
        vpert_cross_beq = self.fielddict['vpert'].cross(self.fielddict['beq'],dmod=0)
        v0_cross_bpert = self.fielddict['ve0'].cross(self.fielddict['bpert'],dmod=0)
        vpert_cross_b0 = self.fielddict['vpert'].cross(self.fielddict['be0'],dmod=0)
        vpert_cross_bpert = self.fielddict['vpert'].cross(self.fielddict['bpert'],dmod=0)

        v_cross_b_eq = self.fft(veq_cross_bpert+vpert_cross_beq)
        v_cross_b_n0 = self.fft(v0_cross_bpert+vpert_cross_b0)
        v_cross_b_pert = self.fft(vpert_cross_bpert)

        fac=1.0 #fac should be one, keep as 2 for convergence testing
        self.powerFluxDict['vxbeq'] = fac * \
            self.dotPlusCc(v_cross_b_eq,self.fielddict['jfour'][0:3])
        self.powerFluxDict['vxbn0'] = fac * \
            self.dotPlusCc(v_cross_b_n0,self.fielddict['jfour'][0:3])
        self.powerFluxDict['vxbp'] = fac * \
            self.dotPlusCc(v_cross_b_pert,self.fielddict['jfour'][0:3])

        fac=-1.0
        nfour=self.fielddict['jfour'].shape[-1]
        self.powerFluxDict['etajp'] = fac * \
            self.mu0 * self.elecd * \
            np.nan_to_num(self.fielddict['diff_shape'][0].data[...,:nfour]) * \
            self.dotPlusCc(self.fielddict['jfour'][0:3],
                           self.fielddict['jfour'][0:3])

        fac=1.0
        div_pie_pert = self.fft(self.fielddict['ndivPiepert']) #div Pi_e /ne
        self.powerFluxDict['divPie'] =fac * \
            self.dotPlusCc(div_pie_pert,self.fielddict['jfour'][0:3])

        fac=1.0 #fac should be one, keep as 2 for convergence testing
        self.powerFluxDict['jxbeq'] = fac * \
            self.dotPlusCc(j_cross_b_eq,self.fielddict['vfour'][0:3])
        self.powerFluxDict['jxbn0'] = fac * \
            self.dotPlusCc(j_cross_b_n0,self.fielddict['vfour'][0:3])
        self.powerFluxDict['jxbp'] = fac * \
            self.dotPlusCc(j_cross_b_pert,self.fielddict['vfour'][0:3])

        gradp_pert = self.fft(self.fielddict['ppert'].grad(dmod=0))
        fac = -1.0 #I think fac should be -1, keep as 2 for convergence testing
        self.powerFluxDict['ngpp'] = fac * \
            self.dotPlusCc(gradp_pert,self.fielddict['vfour'][0:3])

        div_pi = self.calculate_viscositiy()
        fac = - 1.0
        if div_pi is not None:
            div_pi_pert = self.fft(div_pi)
            self.powerFluxDict['divpip'] = fac * \
                self.dotPlusCc(div_pi_pert,self.fielddict['vfour'][0:3])

        fac= -1.0
        div_pii_pert = self.fft(self.fielddict['ndivPiipert'])
        self.powerFluxDict['divPii'] =fac * \
            self.dotPlusCc(div_pii_pert,self.fielddict['vfour'][0:3])

        eFour = self.fft(self.fielddict['epert'],type='v')
        curlEFour = self.fft(self.fielddict['epert'].curl(dmod=0),type='v')

        fac1 = -1.0/self.mu0
        fac2 = 1.0
        self.powerFluxDict['poynting'] = fac1 * \
            self.dotPlusCc(self.fielddict['bfour'][0:3],curlEFour[0:3]) + \
            fac2 * self.dotPlusCc(eFour[0:3],self.fielddict['jfour'][0:3])
