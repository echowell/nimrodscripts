#!/usr/bin/env python3

from plot_nimrod import PlotNimrod as pn
import f90nml
import numpy as np
from eval_nimrod import *
from field_class import *
from fsa import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,splev,UnivariateSpline
import os
import argparse

def readRaw(rawFile):
  datadict={}
  with open(rawFile, 'r') as f:
    while True:
      dataname = f.readline()
      if not dataname:
        break
      dataname = dataname.split()
      datalen = int(dataname[0])
      dataname = dataname[2]
      npdata=np.zeros((3,datalen),np.float)
      for line in range(datalen):
        datastr=f.readline().split()
        npdata[0][line]=float(datastr[0])
        npdata[1][line]=float(datastr[1])
        npdata[2][line]=float(datastr[2])
      datadict[dataname]=npdata
  return datadict

def basefsa(rzc, y, dy, eval_nimrod, fdict):
    '''
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    Set neq to number of outputs in FSA call
    and fill dy[4:4+neq]
    '''
    isurf=fdict.get('isurf')
    bextrema=fdict['bextrema']
    bigr=fdict['bigr']
    n = eval_nimrod.eval_field('n', rzc, dmode=0, eq=2)
    dy[4] = n[0]*dy[2] # ne
    dy[5] = dy[4] # nd #### n is not a vector for me
    dy[6] = 0.
    ti = eval_nimrod.eval_field('ti', rzc, dmode=0, eq=2)
    dy[7] = ti[0]*dy[2] # ti
    te = eval_nimrod.eval_field('te', rzc, dmode=0, eq=2)
    dy[8] = te[0]*dy[2] # te
    bf = eval_nimrod.eval_field('b', rzc, dmode=1, eq=2)
    B = Vector(bf, rzc, torgeom=True, dmod=1)
    bsq = B.dot(B,dmod=0).data
    dy[9] = bsq*dy[2] # B**2
    dy[10] = (B.hat(dmod=0).dot(grad(B.mag())).data)**2*dy[2] # (b.grad(|B|))**2
    dy[11] = rzc[0]*bf[2]*dy[2] # R B_Phi
    bmag = np.sqrt(bsq)
    bextrema[0,isurf] = min(bextrema[0,isurf], bmag)
    bextrema[1,isurf] = max(bextrema[1,isurf], bmag)
    bigr[isurf] = max(bigr[isurf], rzc[0])
    vf = eval_nimrod.eval_field('v', rzc, dmode=0, eq=2)
    dy[12] = vf[2]*dy[2]/rzc[0] # omega
    dy[13] = (vf[0]*bf[0]+vf[1]*bf[1])*dy[2]/np.sqrt(bf[0]*bf[0]+bf[1]*bf[1]) # Kpol
    return dy

def trapfsa(rzc, y, dy, eval_nimrod,fdict):
    '''
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    Set neq to number of outputs in FSA call
    and fill dy[4:4+neq]
    '''
    bf = eval_nimrod.eval_field('b', rzc, dmode=0, eq=2)
    B = Vector(bf, rzc, torgeom=True, dmod=0)
    bmag = B.mag().data
    bave=fdict['bave']
    nlam=fdict['nlam']
    lam=fdict['lam']
    rzo=fdict['rzo_copy']
    dy[4:4+nlam] = np.sqrt(1.0 - lam[:]*bmag/bave)
    dy[4+nlam:4+2*nlam] = dy[2]/(dy[4:4+nlam])
    dy[4:4+nlam]*=dy[2] #note dy[4:4+nlam] is used in above calc
    dy[4+2*nlam] = dy[2]*bmag/bave
    dy[4+2*nlam+1] = dy[2]*np.sqrt((rzc[0]-rzo[0])**2+(rzc[1]-rzo[1])**2)
    return dy

# load dump file
def neoclassical_calculator(dumpfile):
    dumpfile = 'dumpgll.00000.h5'
    nml = f90nml.read('nimrod.in')
    gmt = nml['grid_input']['geom']
    if gmt == 'tor':
        gmt=True
    else:
        gmt=False
    eval_nimrod = EvalNimrod(dumpfile, fieldlist='nvptbj') #, (e not in dumpfile coord='xyz')


    rzo = find_pf_null(eval_nimrod, [1.7, -0.2, 0])

    md=3.3435860e-27
    me=9.1093837015e-31
    #%zc=6
    #%mc=1.9944235e-26
    echrg=1.609e-19
    kboltz=echrg
    kb=1.609e-19
    eps0=8.85418782e-12
    #%
    nsurf = 150 # FSA surfaces
    bextrema = np.zeros([2,nsurf])
    bextrema[0,:] =  np.inf # min
    bextrema[1,:] = -np.inf # max
    bigr = np.zeros([nsurf])
    bigr[:] = -np.inf # max
    #%

################################################################################
# Calculate basic fsa quantitites
################################################################################
    dpow=1.0
    fsafilename = 'fsa.npz'
    if os.path.exists(fsafilename):
        fsaDict = np.load(fsafilename)
        dvar = fsaDict['arr_0']
        yvars = fsaDict['arr_1']
        contours = fsaDict['arr_2']
        bextrema = fsaDict['arr_3']
        bigr = fsaDict['arr_4']
    else:
#        dvar, yvars, contours = FSA(eval_nimrod, rzo, basefsa, 10, nsurf=nsurf, \
#                                    depvar='eta', dpow=dpow, rzx=[1.3, -1.14, 0],
#                                    bextrema=bextrema, bigr=bigr)
        dvar, yvars, contours = FSA(eval_nimrod, rzo, basefsa, 10, nsurf=nsurf, \
                                    depvar='eta', dpow=dpow,
                                    bextrema=bextrema, bigr=bigr)
        fsaArr = [dvar, yvars, contours, bextrema, bigr]
        np.savez(fsafilename,*fsaArr)

    # Determine where the FSA failed
    iend=-1
    while np.isnan(yvars[:,iend]).any():
        iend -= 1
    iend += yvars.shape[1]+1

################################################################################
# Calcualte trapped and passing fractions
################################################################################
    trapfilename = 'trap.npz'
    if os.path.exists(trapfilename):
        trapDict = np.load(trapfilename)
        f_pass = trapDict['arr_0']
        eps = trapDict['arr_1']
    else:
        # Arrays for passing/trapped fractions
        # ind 0 - Helander and Sigmar Eq. 11.24
        #   f_t = 1 - (3/4)*int_0^(Bave/Bmax) dlambda lambda/<SQRT(1-lambda*B/Bave)>
        # ind 1 - Lin-Liu and Miller (1995)
        #   f_tl = 1 - (3/4)*int_0^(Bave/Bmax) dlambda lambda*<1/SQRT(1-lambda*B/Bave)>
        # ind 2 - f_tu from Lin-Liu and Miller (1995)
        #   f_tu = 1 - (3/4)*int_0^(Bave/Bmax) dlambda lambda/SQRT(1-lambda*<B/Bave>)
        # int 3 - f_te from inverse aspect ratio as in B5
        #   f_c ~ 1 - 1.46*sqrt(eps) + 0.46*eps*sqrt(eps)
        f_pass = np.zeros([4,iend])
        eps = np.zeros([iend])
        # integrate from 0 to bmax
        for ii in range(iend):
            nlam = 100
            lam, weights = np.polynomial.legendre.leggauss(nlam)
            bave = np.sqrt(yvars[5,ii]) # sqrt(<B^2>)
            lam += 1
            lam *= bave/(2.0*bextrema[1,ii])
            weights *= bave/(2.0*bextrema[1,ii])
            rzp = [contours[0,0,ii], contours[1,0,ii], 0]
            intgr, contour = FSA(eval_nimrod, rzo, trapfsa, 2*nlam+2, nsurf=1,
                                 depvar='eta', rzp=rzp, bave=bave, nlam=nlam,
                                 lam=lam, rzo_copy=rzo)
            f_pass[0,ii] = 0.75*np.sum(weights*lam/intgr[0:nlam]) #Callen Eqn B5
            f_pass[1,ii] = 0.75*np.sum(weights*lam*intgr[nlam:2*nlam])
            f_pass[2,ii] = 0.75*np.sum(weights*lam/np.sqrt(1.0-lam*intgr[2*nlam]))
            eps[ii] = intgr[2*nlam+1]/rzo[0]
            f_pass[3,ii] = 1 + (-1.46 + 0.46*eps[ii])*np.sqrt(eps[ii])
            print(ii,dvar[1,ii],f_pass[:,ii])
        trapArr = [f_pass,eps]
        np.savez(trapfilename,*trapArr)
    f_trap = 1.0 - f_pass[:,:]

################################################################################
# Plot fsa quantities
################################################################################
    ne = yvars[0,:iend]
    nd = yvars[1,:iend]
    ti = yvars[3,:iend]
    te = yvars[4,:iend]
    fsabsq = yvars[5,:iend]
    fsabdgrBsq = yvars[6,:iend]
    rbphi = yvars[7,:iend]
    omega = yvars[8,:iend]
    kpol = yvars[9,:iend]
    rhon = dvar[1,:iend]
    psi = dvar[2,:iend]
    psix = dvar[2,-1]
    q = np.fabs(dvar[7,:iend])
    bigr = bigr[:iend]
    #%

    rhoofq=interp1d(q,rhon)
    rhoq2=rhoofq(2)
    rhoq3=rhoofq(3)
    rhoq4=rhoofq(4)
    print(rhoq2,rhoq3,rhoq4)

    #pn.plot_scalar_line(None, q, flabel=r'q',
    #                    xvar=rhon, xlabel=r'$\rho_N$', ylabel='',legend_loc='upper left')

    #%
    # Plot trapped fraction
    fig_size = [12,6.75]
    fig,ax = plt.subplots(figsize=fig_size)
    ax.axvline(rhoq2, ls=':')
    ax.axvline(rhoq3, ls=':')
    ax.axvline(rhoq4, ls=':')
    pn.plot_scalar_line(None, f_trap[0,:], flabel=r'$f_t$',
                        f2=f_trap[3,:], f2label=r'$f_{t}$ approx',
                        xvar=rhon, xlabel=r'$\rho_N$', ylabel=r'$f_t$',
                        style='varied',legend_loc='upper left',ax=ax)

    ft=interp1d(rhon,f_trap[0,:])
    ft_approx=interp1d(rhon,f_trap[3,:])
    print(f"q={2} f_t = {ft(rhoq2)} f_t approx = {ft_approx(rhoq2)}")
    print(f"q={3} f_t = {ft(rhoq3)} f_t approx = {ft_approx(rhoq3)}")
    print(f"q={4} f_t = {ft(rhoq4)} f_t approx = {ft_approx(rhoq4)}")

    # Plot fsa quants
    #pn.plot_scalar_line(None, fsabsq, flabel=r'\langle B^2 \rangle',
    #                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper left')
    #pn.plot_scalar_line(None, fsabdgrBsq, flabel=r'\langle(\mathbf{b}\cdot\nabla B)^2\rangle',
    #                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper right')
    fsa_approx = eps**2/(2*rzo[0]**2*q**2)
    #pn.plot_scalar_line(None, fsabdgrBsq/fsabsq,
    #                    flabel=r'\langle(\mathbf{b}\cdot\nabla B)^2\rangle/\langle B^2 \rangle',
    #                    f2=fsa_approx, f2label='\epsilon^2/(2 R_0^2 q^2)',
    #                    xvar=rhon, xlabel=r'\rho_N', ylabel='m^{-2}',legend_loc='upper right')

    fig_size = [12,6.75]
    fig,ax = plt.subplots(figsize=fig_size)
    ax.axvline(rhoq2, ls=':')
    ax.axvline(rhoq3, ls=':')
    ax.axvline(rhoq4, ls=':')
    pn.plot_scalar_line(None, eps, flabel=r'$\epsilon$',
                        xvar=rhon, xlabel=r'$\rho_N$', ylabel='',
                        legend_loc='upper right',ax=ax)

    eps_int=interp1d(rhon,f_trap[0,:])
    print(f"q={2} $\epsilon$ = {eps_int(rhoq2)}")
    temp = 1.46 * np.sqrt(.242)-0.46*np.power(.242,1.5)
    print(f"using epsilon b f_t = {temp} at q=2")
    #%# computed quantities
    #%### nustar ###
    #%

    for species in ['ion','electron']:
        lnLambda=24.-np.log(np.sqrt(ne/10**6)/te)
        sqrt2 = np.sqrt(2)
        if species == 'ion':
            z = 0.0
            ms=md
            ts=ti
            ns=nd
            nu_coef = 4.0 * np.sqrt(np.pi) * echrg**4 * lnLambda \
                    /((4.0 * np.pi * eps0)**2 * kb**1.5 * 3.0 * np.sqrt(md))
            eta_coef = 5/(sqrt2*12.0)*(205.0/(48.0*sqrt2))/(89/48) #z_star =0
            nustauss = 1.0/sqrt2
        elif species == 'electron':
            z = 1.0
            ms=me
            ts=te
            ns=ne
            nu_coef = 4.0 * np.sqrt(2*np.pi) * echrg**4 * lnLambda \
                    /((4.0 * np.pi * eps0)**2 * kb**1.5 * 3.0 * np.sqrt(me))

            eta_coef = 5/12.0 * \
            (17./4.*z**2+205*z/(48*sqrt2))/ \
            (2*z**2+301*z/(48*sqrt2)+89/48)
            nustauss = z
        else:
            print("species unknown")
            raise ValueError
        nu_s=nu_coef * ns/np.sqrt(ts**3)
        vt_s = np.sqrt(ts*kb/ms)
        lambda_s = vt_s/nu_s
        eta00_s = ms*ns*nu_s*lambda_s**2 # A17

        D = 1.2 * (2* np.power(z,2) + 301/(48.*sqrt2)*z+89./48.)




    # plot parallel viscosity as a diffusivity
    #    pn.plot_scalar_line(None, eta00_s/(md*nd), flabel=r'\eta^s_{00}/(m_d n_d)',
    #                    xvar=rhon, xlabel=r'\rho_N', ylabel='m^2/s', legend_loc='upper right')

        omegat_s = vt_s/(rzo[0]*q) # B10
        nustar_s = f_trap[0,:]/(2.92*f_pass[0,:]) * nu_s*omegat_s/vt_s**2 * fsabsq/fsabdgrBsq # B11
        nustar_s_aprx = nu_s*rzo[0]*q/(eps**(1.5)*vt_s) #
        #%


        # plot nustar
        if species=='ion':
            flabel=r'$\nu_{*i}$'
        else:
            flabel=r'$\nu_{*e}$'
        fig,ax = plt.subplots(figsize=fig_size)
        ax.axvline(rhoq2, ls=':')
        ax.axvline(rhoq3, ls=':')
        ax.axvline(rhoq4, ls=':')
        ax.axhline(y=1.0, ls='--',color='k')
        ax.set_ylim([0,10])
        pn.plot_scalar_line(None, nustar_s, flabel=flabel,
                            f2=eps**(-1.5), f2label=r'$\epsilon^{-3/2}$',
                            xvar=rhon, xlabel=r'$\rho_N$', ylabel='',
                            style='varied',legend_loc='upper right',ax=ax)

        nut=interp1d(rhon,nustar_s)
        nut_approx=interp1d(rhon,nustar_s_aprx)
        print(f"q={2} nu_* = {nut(rhoq2)} nu_* approx = {nut_approx(rhoq2)}")
        print(f"q={3} nu_* = {nut(rhoq3)} nu_* approx = {nut_approx(rhoq3)}")
        print(f"q={4} nu_* = {nut(rhoq4)} nu_* approx = {nut_approx(rhoq4)}")

    ### NC ion poloidal flow ###



        tauss = nustauss/nu_s

        K00B_s = (z + sqrt2-np.log(1+sqrt2))/nustauss

        K01B_s = (z + 1.0/sqrt2)/nustauss # isn't this just 1?
        K11B_s = (2.0*z + 9.0/(4.0*sqrt2))/nustauss
        K00P_s = np.sqrt(np.pi)
        K01P_s = 3 * np.sqrt(np.pi)
        K11P_s = 12 * np.sqrt(np.pi)
        K00PS_s = (17.0*z/4.0 + 205.0/(sqrt2*48.0))/D
        K01PS_s = (7.0/2.0)*(23.0*z/4.0 + 241.0/(sqrt2*48.0))/D
        K11PS_s = (49.0/4.0)*(33.0*z/4.0 + 325.0/(sqrt2*48.0))/D

        K00tot_s = K00B_s/(1.0 + np.sqrt(nustar_s) + 2.92*nustar_s*K00B_s/K00P_s) \
              /(1.0 + 2.0*K00P_s/(3.0*omegat_s*tauss*K00PS_s))
        K01tot_s = K01B_s/(1.0 + np.sqrt(nustar_s) + 2.92*nustar_s*K01B_s/K01P_s) \
              /(1.0 + 2.0*K01P_s/(3.0*omegat_s*tauss*K01PS_s))
        K11tot_s = K11B_s/(1.0 + np.sqrt(nustar_s) + 2.92*nustar_s*K11B_s/K11P_s) \
              /(1.0 + 2.0*K11P_s/(3.0*omegat_s*tauss*K11PS_s))

        muiratio = 5.0/2.0 - K01tot_s/K00tot_s
#        mu00_s = K00tot_s * nu_s * tauss * f_trap[0,:]/f_pass[0,:] #TODO

        mu00_s = K00tot_s * nu_s * f_trap[0,:]/f_pass[0,:] #TODO remove tauss
        mu00_s_norm = K00tot_s* f_trap[0,:]/f_pass[0,:]


        if species=='ion':
            flabel=r'$\mu_{i}$'
        else:
            flabel=r'$\mu_{e}$'
        fig_size = [12,6.75]
        fig,ax = plt.subplots(figsize=fig_size)
        ax.axvline(rhoq2, ls=':')
        ax.axvline(rhoq3, ls=':')
        ax.axvline(rhoq4, ls=':')
        pn.plot_scalar_line(None, mu00_s, flabel=flabel,
                        xvar=rhon, xlabel=r'$\rho_N$', ylabel=r'Damping Frequency [1/s]',
                        legend_loc='upper left',ax=ax)


        mus=interp1d(rhon,mu00_s)
        mus_norm=interp1d(rhon,mu00_s_norm)
        print(f"q={2} mu_s = {mus(rhoq2)} mu_s/nu approx = {mus_norm(rhoq2)}")
        print(f"q={3} mu_s = {mus(rhoq3)} mu_s/nu approx = {mus_norm(rhoq3)}")
        print(f"q={4} mu_s = {mus(rhoq4)} mu_s/nu approx = {mus_norm(rhoq4)}")

#        mu01_s = 5.0/2.0*K00tot_s - K01tot_s
#        mu11_s = K11tot_s - 5.0*K01tot_s + 25.0/4.0*K00tot_s
#        nu11_s = sqrt2 # lots of assumptions here







    # Setup FSA quants spline as function of rhon
    #zvars = np.zeros((2, )+yvars.shape) # psin,q + yvars
    #zvars[0,:] = dvar[0,:] # psin
    #zvars[1,:] = dvar[7,:] # q
    #zvars[2:,:] = yvars[:,:]
    #splrho = interp1d(rhon, zvars, kind='cubic')

    #pn.plot_scalar_line(None, yvars[0,:], flabel=r'p', xvar=dvar[1,:], ylabel=r'pressure (Pa)', xlabel=r'\rho_N', legend_loc='upper right')
    #pn.plot_scalar_line(None, yvars[1,:], flabel=r'n', xvar=dvar[1,:], ylabel=r'density (m^{-3})', xlabel=r'\rho_N', legend_loc='upper right')

    ##2D plots -- forces
    #xres=250; yres=250 # need y=1000
    #grid = pn.grid_2d_gen([1.14, -1.19, 0], [1.14, 1.02, 0], [2.31, -1.19, 0], xres, yres)
    #
    #maskfilename = 'mask'+str(xres)+'x'+str(yres)+'.npy'
    #if os.path.exists(maskfilename):
    #    mask = np.load(maskfilename)
    #else:
    #    mask = np.zeros((1,grid.shape[1],grid.shape[2]))
    #    for ix in range(grid.shape[1]):
    #        for iy in range(grid.shape[2]):
    #            if contourContains(Rsep,Zsep,grid[0,ix,iy],grid[1,ix,iy]):
    #                mask[0,ix,iy] = 1.0
    #            else:
    #                mask[0,ix,iy] = np.nan
    #    np.save(maskfilename,mask)
    #
    ##fval = eval_nimrod.eval_field('ti', grid, dmode=0, eq=2)
    ##ti = Scalar(fval, grid, torgeom=gmt, dmod=0)
    ##fval = eval_nimrod.eval_field('n', grid, dmode=1, eq=2)
    ##ne = Scalar(fval, grid, torgeom=gmt, dmod=1, nqty=0)
    ##nd = Scalar(fval, grid, torgeom=gmt, dmod=1, nqty=1)
    ##nc = Scalar(fval, grid, torgeom=gmt, dmod=1, nqty=2)
    #fval = eval_nimrod.eval_field('v', grid, dmode=2, eq=2)
    #v = Vector(fval, grid, torgeom=gmt, dmod=2)
    #fval = eval_nimrod.eval_field('b', grid, dmode=1, eq=2)
    #B = Vector(fval, grid, torgeom=gmt, dmod=1)
    #
    #divv = 2.0/3.0 * div(v)
    #pn.plot_scalar_plane(grid, mask[0]*divv.data, ax=setBdryPlot())
    #BdcvcB = B.dot( curl( v.cross(B) ) ) / B.dot(B)
    #pn.plot_scalar_plane(grid, mask[0]*BdcvcB.data, ax=setBdryPlot())
    #PiParFast = BdcvcB + divv
    #pn.plot_scalar_plane(grid, mask[0]*PiParFast.data, ax=setBdryPlot(), fmin=-10000., fmax=10000.)

    # 1D force-balance test
    grid = pn.grid_1d_gen([1.76821, -0.0188439, 0], [2.7, -0.0188439, 0], 1000)
    fval = eval_nimrod.eval_field('p', grid, dmode=2, eq=1)
    p = Scalar(fval, grid, torgeom=gmt, dmod=2)
    fval = eval_nimrod.eval_field('b', grid, dmode=1, eq=1)
    B = Vector(fval, grid, torgeom=gmt, dmod=1)
    fval = eval_nimrod.eval_field('j', grid, dmode=1, eq=1)
    j = Vector(fval, grid, torgeom=gmt, dmod=1)
    fval = eval_nimrod.eval_field('v', grid, dmode=1, eq=1)
    v = Vector(fval, grid, torgeom=gmt, dmod=1)
    fval = eval_nimrod.eval_field('n', grid, dmode=1, eq=1)
    nd = Scalar(fval, grid, torgeom=gmt, dmod=1, nqty=1)
    fval = eval_nimrod.eval_field('ti', grid, dmode=1, eq=1)
    ti = Scalar(fval, grid, torgeom=gmt, dmod=1)
    #nc = Scalar(fval, grid, torgeom=gmt, dmod=1, nqty=2)
    # computed quantities
    jxb=j.cross(B)
    grp=p.grad()
    mc=1.9944235e-26
    #rhovdgrdv = (md*nd+mc*nc)*v.dot(grad(v))
    rhovdgrdv = (md*nd)*v.dot(grad(v))

    force=jxb-grp
#    pn.plot_vector_line(grid, force.data, flabel=r'Force', f2=jxb.data,
#                        f2label=r'$J\times B$',
#                        f3=grp.data, f3label=r'$\nabla p$',
#                        f4=rhovdgrdv.data, f4label=r'$v\cdot \nabla v$',
#                        comp='perp', style='varied',
#                        legend_loc='lower left', ylabel=r'Force density N/m^3')

#    pn.plot_scalar_line(grid, ti.data, flabel=r'T_i', style='varied',
#                        legend_loc='lower left', ylabel=r'Temperature eV')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Perform mu calculations from UW-CPTC-09-6R'
    )
    parser.add_argument('file',help='dumpfile',default='dumpgll.00000.h5')
    args = vars(parser.parse_args())
    neoclassical_calculator(dumpfile=args['file'])
