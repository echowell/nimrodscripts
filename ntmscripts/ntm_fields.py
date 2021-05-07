#!/usr/bin/env python3

import f90nml
import eval_nimrod as eval
import field_class as fc
import numpy as np
from scipy.fft import fft, ifft
import plot_nimrod as pn
import matplotlib.pyplot as plt
import time
import nim_timer as timer

class ntmfields:
    ''' ntm fields is a class for reading/storing nimrod data on a mesh
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
        self.evalb_timer=0.0
        self.fft_timer=0.0
        self.ndmode=0
        self.neq=False
        self.npert=False
        self.vdmode=0
        self.veq=False
        self.vpert=False
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
    def fft(self,pert,axis=-1):
        ''' NIMROD stores it's field data as f(phi) = sum_{-n}^n f_n exp(inphi)
            This implies that the normalization 1/N should be done in the transform
            from physical space to fourier space
            This is the one option that scipi.fft does not support, but why?
        '''
        fpert = fft(pert.data,axis=axis,norm=None)/pert.data.shape[axis]
        return fpert

    def set_method(self,method):
      if method == "induction":
        self.ndmode=1
        self.neq=True
        self.npert=False
        self.vdmode=1
        self.veq=True
        self.vpert=True
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
        if self.neq:
            if len(rzp.shape)==2:
                field=self.eval_symm('n',rzp,dmode=self.ndmode,eq=1)
                self.fielddict['neq']=fc.Scalar(field,rzp,self.ndmode,True)
            else:
                field=self.eval.eval_field('n', grid, dmode=self.ndmode, eq=1)
                self.fielddict['neq']=fc.Scalar(field,rzp,self.ndmode,True)
        if self.npert:
            field=self.eval.eval_field('n', grid, dmode=self.ndmode, eq=0)
            self.fielddict['npert']=fc.Scalar(field,rzp,self.ndmode,True)
        if fft:
            self.fielddict['nfour']=self.fft(field)
        return None

    @timer.timer_func
    def eval_v(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if self.veq:
            if len(rzp.shape)==2:
                field=self.eval_symm('v',rzp,dmode=self.vdmode,eq=1)
                self.fielddict['veq']=fc.Vector(field,rzp,self.vdmode,True)
            else:
                field=self.eval.eval_field('v', grid, dmode=self.vdmode,eq=1)
                self.fielddict['veq']=fc.Vector(field,rzp,self.vdmode,True)
        if self.vpert:
            field=self.eval.eval_field('v', grid, dmode=self.vdmode, eq=0)
            self.fielddict['vpert']=fc.Vector(field,rzp,self.vdmode,True)
        if fft:
            self.fielddict['vfour']=self.fft(field)
        return None

    @timer.timer_func
    def eval_b(self,grid=None,fft=False):
        grid,rzp=self.get_gridrzp(grid)
        if len(rzp.shape)==2:
            field=self.eval_symm('b',rzp,dmode=1,eq=1)
            self.fielddict['beq']=fc.Vector(field,rzp,1,True)
            field=self.eval_symm('b',rzp,dmode=1,eq=3)
            self.fielddict['b0']=fc.Vector(field,rzp,1,True)
        else:
            field=self.eval.eval_field('b',grid,dmode=1,eq=1)
            self.fielddict['beq']=fc.Vector(field,rzp,1,True)
            field=self.eval.eval_field('b',grid,dmode=1,eq=3)
            self.fielddict['b0']=fc.Vector(field,rzp,1,True)

        field=self.eval.eval_field('b',grid,dmode=1, eq=0)
        self.fielddict['bpert']=fc.Vector(field,rzp,1,True)
        if fft:
          self.fielddict['bfour']=self.fft(field)
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
        if len(rzp.shape)==2:
            field=self.eval_symm('j',rzp,dmode=1,eq=1)
            self.fielddict['jeq']=fc.Vector(field,rzp,1,True)
        else:
            field=self.eval.eval_field('j', grid, dmode=1, eq=1)
            self.fielddict['jeq']=fc.Vector(field,rzp,1,True)
        field=self.eval.eval_field('j', grid, dmode=1, eq=0)
        self.fielddict['jpert']=fc.Vector(field,rzp,1,True)
        if fft:
            self.fielddict['jfour']=self.fft(field)
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
    def eval_p(self,grid=None):
          grid,rzp=self.get_gridrzp(grid)
          if len(rzp.shape)==2:
              field=self.eval_symm('p',rzp,dmode=2,eq=1)
              self.fielddict['peq']=fc.Scalar(field,rzp,2,True)
          else:
              field=self.eval.eval_field('p',grid,dmode=2,eq=1)
              self.fielddict['peq']=fc.Scalar(field,rzp,2,True)
          field=self.eval.eval_field('p',grid,dmode=2,eq=0)
          self.fielddict['ppert']=fc.Scalar(field,rzp,2,True)
          self.fielddict['pfour']=self.fft(field)
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

    #begin coding of field equations
    def continuity(self,grid=None):
        if not 'neq' in self.fielddict:
          self.eval_n(grid)
        if not 'veq' in self.fielddict:
          self.eval_v(grid)

        nddiff = self.physics_nml.get('nd_diff',0.0)
        veq=self.fielddict['veq']
        vpert=self.fielddict['vpert']
        neq=self.fielddict['neq']
        npert=self.fielddict['npert']

        advecteq=veq.dot(neq.grad())
        advectlin=veq.dot(npert.grad())+vpert.dot(neq.grad())
        advectnon=vpert.dot(npert.grad())
        compeq=veq.div()*neq
        #todo add linear and nonlinear compresson
        nddiffusion=(nddiff*npert.grad()).div()
        self.ndict['advecteq']=advecteq
        self.ndict['advectlin']=advectlin
        self.ndict['advectnon']=advectnon
        self.ndict['dompeq']=compeq
        self.ndict['lindiff']=nddiffusion
        return None

    def momentum(self,grid=None):
      #todo
      pass

    def temperature(self,grid=None):
        #todo
        pass

    @timer.timer_func
    def ohms(self,grid=None):
        grid,rzp=self.get_gridrzp(grid)
        self.edict={}
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

        if 'veq' not in self.fielddict:
            self.eval_v(grid)
        if 'beq' not in self.fielddict:
            self.eval_b(grid)
        if 'jeq' not in self.fielddict:
            self.eval_j(grid)
        if 'diff_shape' not in self.fielddict:
            self.eval_diff(grid)
        veq=self.fielddict['veq']
        vpert=self.fielddict['vpert']
        beq=self.fielddict['beq']
        bpert=self.fielddict['bpert']
        jeq=self.fielddict['jeq']
        jpert=self.fielddict['jpert']
        diff_shape=self.fielddict['diff_shape']

        do_eq=False
        if do_eq:
            veqbeq=-veq.cross(beq)
            etajeq=self.mu0*self.elecd*diff_shape[0]*jeq
            self.edict['veqbeq']=veqbeq
            self.edict['etajeq']=etajeq

        vpbeq=-vpert.cross(beq)
        veqbp=-veq.cross(bpert)
        vpbp=-vpert.cross(bpert)
        etajpert=self.mu0*self.elecd*diff_shape[0]*jpert
        self.edict['vpbeq']=vpbeq
        self.edict['veqbp']=veqbp
        self.edict['vpbp']=vpbp
        self.edict['etajpert']=etajpert

        if ohmslaw in ['mhd&hall','2fl']:
        #TODO
            print("Hall and two fluid Ohms law are not yet supported")
            raise ValueError
            if 'neq' not in self.fielddict:
                self.eval_n(grid)
            if 'peq' not in self.fielddict:
                self.eval_p(grid)
            neq=self.fielddict['neq']
            npert=self.fielddict['npert']



        if neoe_flag in ['gianakon']:
            if 'neq' not in self.fielddict:
                self.eval_n(grid)
            if 'fsa_beq2' not in self.fielddict:
                self.eval_fsa_beq2(grid)
            neq=self.fielddict['neq']
            fsa_beq2=self.fielddict['fsa_beq2']

            ephi=[0,0,1]
            etor=fc.basis_vector('p',rzp,torgeom=True)
            bep=beq-beq.dot(etor)*etor
            neo_mask=self.eval_neo_mask(grid)
            coef=self.me*mu_e/(self.qe**2)*neo_mask
            coef1=self.me*mu_e/(self.qe)*neo_mask
            #coef=neo_mask*mu_e

            divpie=coef/neq*fsa_beq2*(jpert.dot(bep))/(beq.dot(bep)**2+1.0e-8) * bep
            self.edict['divpie']=divpie

            divpie1=coef1*fsa_beq2*(jpert.dot(bep))/(beq.dot(bep)**2+1.0e-8) * bep
            self.edict['divpie1']=divpie1

            divpie0=coef/neq*fsa_beq2*(bep.mag())/(beq.dot(bep)**2+1.0e-8) * bep
            self.edict['divpie0']=divpie0

            divpieb=coef1*fsa_beq2*(bep.mag())/(beq.dot(bep)**2+1.0e-8) * bep
            self.edict['divpieb']=divpieb
            
            #pn.PlotNimrod.plot_scalar_plane(rzp[:,:,:,0], fsa_beq2.data[:,:,0])
            #pn.PlotNimrod.plot_scalar_plane(rzp[:,:,:,0], divpie.data[0,:,:,0])
            #pn.PlotNimrod.plot_scalar_plane(rzp[:,:,:,0], divpie.data[1,:,:,0])
        return None

    @timer.timer_func
    def induction(self,grid=None):
        self.ohms(grid)
        self.dbdtdict={}
        for key, value in self.edict.items():
            self.dbdtdict[key]=-value.curl()
            '''
            fig=plt.figure(figsize=(6,8))
            ax=fig.add_subplot(111)
            ax.set_aspect('equal')
            #    levels = numpy.linspace(fmin, fmax, nlvls)
            levels = 50
            cf = ax.contourf(self.grid.rzp[0,:,:,0], self.grid.rzp[1,:,:,0], self.dbdtdict[key].data[0,:,:,0], levels)
            plt.colorbar(cf)
            plt.title(key)
            plt.tight_layout()
            plt.show()
            fig=plt.figure(figsize=(6,8))
            ax=fig.add_subplot(111)
            ax.set_aspect('equal')
            #    levels = numpy.linspace(fmin, fmax, nlvls)
            levels = 50
            cf = ax.contourf(self.grid.rzp[0,:,:,0], self.grid.rzp[1,:,:,0], self.dbdtdict[key].data[1,:,:,0], levels)
            plt.colorbar(cf)
            plt.title(key)
            plt.tight_layout()
            plt.show()
            thisfft=self.fft(self.dbdtdict[key])
            print(self.dbdtdict[key].data.shape)
            print(thisfft.shape)
            print(thisfft[0,50,:,1])
            print(thisfft[0,50,:,-1])
            print(type(thisfft.real),type(self.dbdtdict[key].data))
            print(self.dbdtdict[key].data[0,50,:,0])

            print(np.nanmax(thisfft.real[0,50,:,0]))
            fig=plt.figure(figsize=(6,8))
            ax=fig.add_subplot(111)
            ax.set_aspect('equal')
            #    levels = numpy.linspace(fmin, fmax, nlvls)
            levels = 50
            cf = ax.contourf(self.grid.rzp[0,:,:,0], self.grid.rzp[1,:,:,0], thisfft.real[0,:,:,1], levels)
            plt.colorbar(cf)
            plt.title(f"fft {key}")
            plt.tight_layout()
            plt.show()
            '''
        return None
