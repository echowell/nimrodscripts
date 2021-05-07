#!/usr/bin/env python3

from plot_nimrod import PlotNimrod as pn
import f90nml
from eval_nimrod import *
from field_class import *
from fsa import *
#from particle import *
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d,splev,UnivariateSpline
import os
#import pfileReader

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
      npdata=numpy.zeros((3,datalen),numpy.float)
      for line in range(datalen):
        datastr=f.readline().split()
        npdata[0][line]=float(datastr[0])
        npdata[1][line]=float(datastr[1])
        npdata[2][line]=float(datastr[2])
      datadict[dataname]=npdata
  return datadict

# load dump file
dumpfile = 'dumpgll.00000.h5'
nml = f90nml.read('nimrod.in')
gmt = nml['grid_input']['geom']
if gmt == 'tor':
    gmt=True
else:
    gmt=False
eval_nimrod = EvalNimrod(dumpfile, fieldlist='nvptbj') #, (e not in dumpfile coord='xyz')
print("After eval nimrod")

#%# read sol.grn and prepare plots
#%Rsol = numpy.loadtxt('sol.grn', skiprows=1, max_rows=1)
#%Zsol = numpy.loadtxt('sol.grn', skiprows=3, max_rows=1)
#%Rsep = numpy.loadtxt('sol.grn', skiprows=5, max_rows=1)
#%Zsep = numpy.loadtxt('sol.grn', skiprows=7, max_rows=1)
#%Rwall = numpy.loadtxt('sol.grn', skiprows=9, max_rows=1)
#%Zwall = numpy.loadtxt('sol.grn', skiprows=11, max_rows=1)
#%Rwall = numpy.append(Rwall,Rwall[0])
#%Zwall = numpy.append(Zwall,Zwall[0])
rzo = find_pf_null(eval_nimrod, [1.7, -0.2, 0])
#%def setBdryPlot():
#%    fig, ax = pl.subplots()
#%    fig.set_size_inches(8,12)
#%    ax.plot(Rsol[:],Zsol[:],color='k')
#%    ax.plot(Rsep[:],Zsep[:],color='k')
#%    ax.plot(Rwall[:],Zwall[:],color='k')
#%    ax.set(ylabel="Z (m)", xlabel="R (m)")
#%    ax.ticklabel_format(axis='both', style='sci', scilimits=(10**3,10**-3),
#%                        useOffset=None, useLocale=None, useMathText=True)
#%    pl.tight_layout()
#%    # find o-point
#%    ax.plot(rzo[0],rzo[1],'o')
#%    return ax
#%def contourContains(rcon,zcon,rr,zz):
#%    intersections = 0
#%    for ii in range(rcon.size):
#%        if (zcon[ii]-zz)*(zcon[ii-1]-zz)<0. or zcon[ii]-zz==0.:
#%            if rcon[ii]>rr and rcon[ii-1]>rr:
#%                intersections += 1
#%            elif rcon[ii]<rr and rcon[ii-1]<rr:
#%                continue
#%            else:
#%                if rcon[ii]>rcon[ii-1]:
#%                    if (rr-rcon[ii-1])/(rcon[ii]-rcon[ii-1]) \
#%                            +(zz-zcon[ii])/(zcon[ii-1]-zcon[ii]) < 1.e-8:
#%                        intersections += 1
#%                else:
#%                    if (rr-rcon[ii])/(rcon[ii-1]-rcon[ii]) \
#%                            +(zz-zcon[ii-1])/(zcon[ii]-zcon[ii-1]) < 1.e-8:
#%                        intersections += 1
#%    if intersections%2 == 0:
#%        return False
#%    else:
#%        return True
#%

md=3.3435860e-27
#%zc=6
#%mc=1.9944235e-26
echrg=1.609e-19
kboltz=echrg
kb=1.609e-19
eps0=8.85418782e-12
#%
nsurf = 150 # FSA surfaces
bextrema = numpy.zeros([2,nsurf])
bextrema[0,:] =  numpy.inf # min
bextrema[1,:] = -numpy.inf # max
bigr = numpy.zeros([nsurf])
bigr[:] = -numpy.inf # max
#%
# Do FSA
def basefsa(rzc, y, dy, eval_nimrod, fdict):
    '''
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    Set neq to number of outputs in FSA call
    and fill dy[4:4+neq]
    '''
    isurf=fdict.get('isurf')
    n = eval_nimrod.eval_field('n', rzc, dmode=0, eq=2)
    dy[4] = n[0]*dy[2] # ne
    dy[5] = dy[4] # nd #### n is not a vector for me
    if eval_nimrod.ndnq>1:
      dy[6] = n[2]*dy[2] # nc
    else:
      dy[6] = 0.
    ti = eval_nimrod.eval_field('ti', rzc, dmode=0, eq=2)
    dy[7] = ti[0]*dy[2] # ti
    te = eval_nimrod.eval_field('te', rzc, dmode=0, eq=2)
    dy[8] = te[0]*dy[2] # te
    bf = eval_nimrod.eval_field('b', rzc, dmode=1, eq=2)
    B = Vector(bf, rzc, torgeom=gmt, dmod=1)
    bsq = B.dot(B,dmod=0).data
    dy[9] = bsq*dy[2] # B**2
    dy[10] = (B.hat(dmod=0).dot(grad(B.mag())).data)**2*dy[2] # (b.grad(|B|))**2
    dy[11] = rzc[0]*bf[2]*dy[2] # R B_Phi
    bmag = numpy.sqrt(bsq)
    bextrema[0,isurf] = min(bextrema[0,isurf], bmag)
    bextrema[1,isurf] = max(bextrema[1,isurf], bmag)
    bigr[isurf] = max(bigr[isurf], rzc[0])
    vf = eval_nimrod.eval_field('v', rzc, dmode=0, eq=2)
    dy[12] = vf[2]*dy[2]/rzc[0] # omega
    dy[13] = (vf[0]*bf[0]+vf[1]*bf[1])*dy[2]/numpy.sqrt(bf[0]*bf[0]+bf[1]*bf[1]) # Kpol
    return dy
#%
dpow=1.0
fsafilename = 'fsa.npz'
if os.path.exists(fsafilename):
    fsaDict = numpy.load(fsafilename)
    dvar = fsaDict['arr_0']
    yvars = fsaDict['arr_1']
    contours = fsaDict['arr_2']
    bextrema = fsaDict['arr_3']
    bigr = fsaDict['arr_4']
else:
    dvar, yvars, contours = FSA(eval_nimrod, rzo, basefsa, 10, nsurf=nsurf, \
                                depvar='eta', dpow=dpow, rzx=[1.3, -1.14, 0])
    fsaArr = [dvar, yvars, contours, bextrema, bigr]
    numpy.savez(fsafilename,*fsaArr)

# Determine where the FSA failed
iend=-1
while numpy.isnan(yvars[:,iend]).any():
    iend -= 1
iend += yvars.shape[1]+1
#%
# Do FSA for f_passing
def trapfsa(rzc, y, dy, eval_nimrod, isurf):
    '''
    Flux surface averge quantities (f/bdgrth where y[2]=1/bdgrth)
    Set neq to number of outputs in FSA call
    and fill dy[4:4+neq]
    '''
    bf = eval_nimrod.eval_field('b', rzc, dmode=0, eq=2)
    B = Vector(bf, rzc, torgeom=gmt, dmod=0)
    bmag = B.mag().data
    dy[4:4+nlam] = numpy.sqrt(1.0 - lam[:]*bmag/bave)
    dy[4+nlam:4+2*nlam] = dy[2]/(dy[4:4+nlam])
    dy[4:4+nlam] *= dy[2]
    dy[4+2*nlam] = dy[2]*bmag/bave
    dy[4+2*nlam+1] = dy[2]*numpy.sqrt((rzc[0]-rzo[0])**2+(rzc[1]-rzo[1])**2)
    return dy
#%
trapfilename = 'trap.npz'
if os.path.exists(trapfilename):
    trapDict = numpy.load(trapfilename)
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
    f_pass = numpy.zeros([4,iend])
    eps = numpy.zeros([iend])
    # integrate from 0 to bmax
    for ii in range(iend):
        nlam = 100
        lam, weights = numpy.polynomial.legendre.leggauss(nlam)
        bave = numpy.sqrt(yvars[5,ii]) # sqrt(<B^2>)
        lam += 1
        lam *= bave/(2.0*bextrema[1,ii])
        weights *= bave/(2.0*bextrema[1,ii])
        rzp = [contours[0,0,ii], contours[1,0,ii], 0]
        intgr, contour = FSA(eval_nimrod, rzo, trapfsa, 2*nlam+2, nsurf=1, depvar='eta', rzp=rzp)
        f_pass[0,ii] = 0.75*numpy.sum(weights*lam/intgr[0:nlam])
        f_pass[1,ii] = 0.75*numpy.sum(weights*lam*intgr[nlam:2*nlam])
        f_pass[2,ii] = 0.75*numpy.sum(weights*lam/numpy.sqrt(1.0-lam*intgr[2*nlam]))
        eps[ii] = intgr[2*nlam+1]/rzo[0]
        f_pass[3,ii] = 1 + (-1.46 + 0.46*eps[ii])*numpy.sqrt(eps[ii])
        print(ii,dvar[1,ii],f_pass[:,ii])
    trapArr = [f_pass,eps]
    numpy.savez(trapfilename,*trapArr)
f_trap = 1.0 - f_pass[:,:]

#%#import IPython; IPython.embed()

# set fields as views to the arrays
ne = yvars[0,:iend]
nd = yvars[1,:iend]
#nc = yvars[2,:iend]
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
q = numpy.fabs(dvar[7,:iend])
bigr = bigr[:iend]
#%
pn.plot_scalar_line(None, psi, flabel=r'\psi',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper left')
#%
# Plot trapped fraction
pn.plot_scalar_line(None, f_trap[0,:], flabel=r'f_t',
                    f2=f_trap[1,:], f2label=r'f_{tl}',
                    f3=f_trap[2,:], f3label=r'f_{tu}',
                    f4=f_trap[3,:], f4label=r'f_{t\epsilon}',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',
                    style='varied',legend_loc='upper left')

# Plot fsa quants
pn.plot_scalar_line(None, fsabsq, flabel=r'\langle B^2 \rangle',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper left')
pn.plot_scalar_line(None, fsabdgrBsq, flabel=r'\langle(\mathbf{b}\cdot\nabla B)^2\rangle',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper right')
fsa_approx = eps**2/(2*rzo[0]**2*q**2)
pn.plot_scalar_line(None, fsabdgrBsq/fsabsq,
                    flabel=r'\langle(\mathbf{b}\cdot\nabla B)^2\rangle/\langle B^2 \rangle',
                    f2=fsa_approx, f2label='\epsilon^2/(2 R_0^2 q^2)',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='m^{-2}',legend_loc='upper right')

#%# computed quantities
#%### nustar ###
#%
nc=0
zc=0

lnLambda=24.-numpy.log(numpy.sqrt(ne/10**6)/te)
sqrt2 = numpy.sqrt(2)
coef = 4.0 * numpy.sqrt(numpy.pi) * echrg**4 * lnLambda \
        /((4.0 * numpy.pi * eps0)**2 * kb**1.5 * 3.0 * numpy.sqrt(md))
nui = coef*(nd+sqrt2*nc*zc**2)*lnLambda / numpy.sqrt(ti**3) # A14
vti = numpy.sqrt(ti*kb/md)
lambdai = vti/nui     # A15
zstar = nc*zc**2 / nd # A13
coef = 5*(1+sqrt2*zstar)/(sqrt2*12.0)*((17.0*zstar/4.0)+205.0/(48.0*sqrt2)) \
         /(2*zstar**2 + 301*zstar/(48.0*sqrt2)+89/48)
etai00 = md*nd*nui*lambdai**2 # A17

# plot parallel viscosity as a diffusivity
pn.plot_scalar_line(None, etai00/(md*nd), flabel=r'\eta^i_{00}/(m_d n_d)',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='m^2/s', legend_loc='upper right')

omegati = vti/(rzo[0]*q) # B10
nustari = f_trap[0,:]/(2.92*f_pass[0,:]) * nui*omegati/vti**2 * fsabsq/fsabdgrBsq # B11
nustari_aprx = nui*rzo[0]*q/(eps**(1.5)*vti) #
#%
pn.plot_scalar_line(None, eps, flabel=r'\epsilon',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='upper right')

# plot nustar
pn.plot_scalar_line(None, nustari, flabel=r'\nu_*',
                    f2=nustari_aprx, f2label=r'R_0 q \nu_i / \epsilon^{3/2} v_{Ti}',
                    f3=eps**(-1.5), f3label=r'\epsilon^{-3/2}',
                    xvar=rhon, xlabel=r'\rho_N', ylabel='',
                    style='varied',legend_loc='upper right')
#%
#%### IOL quants ###
#%
#%# Compute the electric potential
#%def int_enorm(delta, rzst, norm):
#%    '''
#%    Integrand for E normal to compute electric potential.
#%    '''
#%    rz = rzst + norm*delta
#%    e = eval_nimrod.eval_field('e', rz, dmode=0, eq=2)
#%    return numpy.dot(e[0:3],norm)
#%norm = numpy.array([contours[0,0,iend]-contours[0,0,0],
#%                    contours[1,0,iend]-contours[1,0,0],0.0])
#%norm = norm/numpy.linalg.norm(norm)
#%pot = numpy.zeros([iend])
#%delta0 = 0.0
#%rzst = rzo
#%for ii in range(iend):
#%    rzp = numpy.array([contours[0,0,ii],contours[1,0,ii],0.0])
#%    delta1 = numpy.linalg.norm(rzst-rzp)
#%    pot[ii], err = quad(int_enorm, delta0, delta1, args=(rzst, norm))
#%    #print(rzp,delta1,pot[ii],err)
#%    delta0 = delta1
#%pot = numpy.cumsum(pot)
#%# plot the potential
#%#pn.plot_scalar_line(None, pot, flabel=r'\Phi',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='V',legend_loc='upper left')
#%# Make a spline vs psi
#%splpot = UnivariateSpline(numpy.fabs(psi), pot, k=4, s=0)
#%d2pot = splpot.derivative(n=2)(numpy.fabs(psi))
#%#pn.plot_scalar_line(None, d2pot, flabel=r'\Phi^{\prime \prime}',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='V/T^2m^4',legend_loc='upper left')
#%
#%fval = eval_nimrod.eval_field('b', rzo, dmode=0, eq=2)
#%B0 = Vector(fval, rzo, torgeom=gmt, dmod=0)
#%omega0 = echrg*B0.mag().data/md
#%S = 1.0 + (rbphi/omega0)**2 * (echrg*d2pot/md)
#%#pn.plot_scalar_line(None, S, flabel=r'S',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='',legend_loc='lower left')
#%omegai = echrg*bextrema[0,:iend]/md
#%term1 = nustari**0.25
#%term2 = omegai*numpy.sqrt(S)/(rbphi*vti*numpy.sqrt(2.0*eps))*numpy.fabs(psi-psix)
#%#pn.plot_scalar_line(None, term1, flabel=r'\nu_*^{1/4}',
#%#                    f2=term2, f2label=r'\Omega_i \sqrt{|S|} |\psi - \psi_x| / F v_{Ti} \sqrt{2\epsilon}',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='',
#%#                    legend_loc='upper right')
#%dndtorb = -2.25*nd*nui/numpy.sqrt(numpy.pi*2.0*S*eps) \
#%            *numpy.exp(-(term1+term2)**2.0)
#%S0 = 1.0
#%term02 = omegai*numpy.sqrt(S0)/(rbphi*vti*numpy.sqrt(2.0*eps))*numpy.fabs(psi-psix)
#%dndtorb0 = -2.25*nd*nui/numpy.sqrt(numpy.pi*2.0*S0*eps) \
#%            *numpy.exp(-(term1+term2)**2.0)
#%#pn.plot_scalar_line(None, dndtorb, flabel=r'\langle\partial n / \partial t\rangle_{orb}',
#%#                    f2=dndtorb0, f2label=r'\langle\partial n / \partial t\rangle_{orb} E=0',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='m^{-3} s^{-1}',
#%#                    legend_loc='lower left')
#%
#%dr=0.01 # 1cm estimated
#%jorbBtor = echrg*dr*dndtorb*rbphi/bigr # estimate, be more careful with BPhi
#%#pn.plot_scalar_line(None, jorbBtor, flabel=r'e \Gamma_{orb} B_{\Phi}',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel='N m^{-3}',
#%#                    legend_loc='lower left')
#%

### NC poloidal flow ###

D = 6.0/5.0 * (2.0*zstar**2+301.0*zstar/(sqrt2*48.0)+89.0/48.0)
nuitauii = 1.0/sqrt2 + zstar
tauii = nuitauii/nui
K00iB = (zstar + sqrt2-numpy.log(1+sqrt2))/nuitauii
K01iB = (zstar + 1.0/sqrt2)/nuitauii # isn't this just 1?
K11iB = (2.0*zstar + 9.0/(4.0*sqrt2))/nuitauii
K00iP = numpy.sqrt(numpy.pi)
K01iP = 3 * numpy.sqrt(numpy.pi)
K11iP = 12 * numpy.sqrt(numpy.pi)
K00iPS = (17.0*zstar/4.0 + 205.0/(sqrt2*48.0))/D
K01iPS = (7.0/2.0)*(23.0*zstar/4.0 + 241.0/(sqrt2*48.0))/D
K11iPS = (49.0/4.0)*(33.0*zstar/4.0 + 325.0/(sqrt2*48.0))/D
K00itot = K00iB/(1.0 + numpy.sqrt(nustari) + 2.92*nustari*K00iB/K00iP) \
          /(1.0 + 2.0*K00iP/(3.0*omegati*tauii*K00iPS))
K01itot = K01iB/(1.0 + numpy.sqrt(nustari) + 2.92*nustari*K01iB/K01iP) \
          /(1.0 + 2.0*K01iP/(3.0*omegati*tauii*K01iPS))
K11itot = K11iB/(1.0 + numpy.sqrt(nustari) + 2.92*nustari*K11iB/K11iP) \
          /(1.0 + 2.0*K11iP/(3.0*omegati*tauii*K11iPS))
muiratio = 5.0/2.0 - K01itot/K00itot
mui00 = K00itot
mui01 = 5.0/2.0*K00itot - K01itot
mui11 = K11itot - 5.0*K01itot + 25.0/4.0*K00itot
nui11 = sqrt2 # lots of assumptions here
csharp = 1.0/(1.0+mui11-mui01*muiratio)/nui11
splti = UnivariateSpline(rhon, ti, k=4, s=0)
splpsi = UnivariateSpline(rhon, psi, k=4, s=0)
dtidpsi = splti.derivative(n=1)(rhon)/splpsi.derivative(n=1)(rhon)
pn.plot_scalar_line(None, dtidpsi, flabel=r'\partial T_i / \partial \psi',
                    xvar=rhon, xlabel=r'\rho_N', ylabel=r'eV/Tm^2',
                    legend_loc='upper left')
U0itheta = -csharp*muiratio*rbphi/(echrg*fsabsq)*kboltz*dtidpsi

#%slicefileName='../../haskey/slice_164988.3525'
#%rfile = readRaw(slicefileName)
#%fitfileName='../../haskey/fit_164988.3525'
#%ffile = pfileReader.pfileReader(fitfileName)
#%vpolcfit   = ffile.getPfileData('V_pol_12C6_FIT(m/s)')
#%vpolcdat   = rfile['V_pol_12C6_SLICE(m/s)']
#%splvpfit = UnivariateSpline(vpolcfit[0], vpolcfit[1], k=4, s=0)
#%v0ipol = numpy.zeros(iend)
#%Uitheta = numpy.zeros(iend)
#%fsaNCforce = numpy.zeros(iend)
#%for ii in range(iend):
#%    rz = numpy.array([bigr[ii], 0, 0])
#%    bf = eval_nimrod.eval_field('b', rz, dmode=1, eq=2)
#%    v0ipol[ii] = bf[1]*U0itheta[ii]
#%    B = Vector(bf, rzo, torgeom=gmt, dmod=1)
#%    bgrB = B.hat(dmod=0).dot(grad(B.mag()))
#%    Uitheta[ii] = splvpfit(rhon[ii])/bf[1]
#%    fsaNCforce[ii] = md*nd[ii]*mui00[ii]*(fsabsq[ii]/fsabdgrBsq[ii])*(bgrB.data)**2*(Uitheta[ii]-U0itheta[ii])
#%#pn.plot_scalar_line(None, U0itheta, flabel=r'U^0_{i\theta}',
#%#                    #f2=Uitheta, f2label=r'U_{i\theta}\;C6\;CER\;fit',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel=r'm/sT',
#%#                    legend_loc='upper left')
#%#fig_size = [12,6.75]
#%#fig, ax = pl.subplots(figsize=fig_size)
#%#ax.errorbar(vpolcdat[0],vpolcdat[1],yerr=vpolcdat[2],fmt='x',color='b', label='C6 CER measurements')
#%#ax.plot(vpolcfit[0],vpolcfit[1],color='r',label='C6 CER fit')
#%#pn.plot_scalar_line(None, v0ipol, flabel=r'U^0_{i\theta} B_p\;(Z=0)',
#%#                    xvar=rhon, xlabel=r'\rho_N', ylabel=r'm/s',
#%#                    legend_loc='lower left', ax=ax)
#%pn.plot_scalar_line(None, fsaNCforce,
#%                    flabel=r'm_i n_i \mu_i \langle B^2_0 \rangle (U_{i\theta} - U_{i\theta}^0)',
#%                    xvar=rhon, xlabel=r'\rho_N', ylabel=r'N m^{-3}',
#%                    legend_loc='upper right')
#%


# Setup FSA quants spline as function of rhon
#zvars = numpy.zeros((2, )+yvars.shape) # psin,q + yvars
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
#    mask = numpy.load(maskfilename)
#else:
#    mask = numpy.zeros((1,grid.shape[1],grid.shape[2]))
#    for ix in range(grid.shape[1]):
#        for iy in range(grid.shape[2]):
#            if contourContains(Rsep,Zsep,grid[0,ix,iy],grid[1,ix,iy]):
#                mask[0,ix,iy] = 1.0
#            else:
#                mask[0,ix,iy] = numpy.nan
#    numpy.save(maskfilename,mask)
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
md=3.3435860e-27
mc=1.9944235e-26
#rhovdgrdv = (md*nd+mc*nc)*v.dot(grad(v))
rhovdgrdv = (md*nd)*v.dot(grad(v))

force=jxb-grp
pn.plot_vector_line(grid, force.data, flabel=r'Force', f2=jxb.data,
                    f2label=r'$J\times B$',
                    f3=grp.data, f3label=r'$\nabla p$',
                    f4=rhovdgrdv.data, f4label=r'$v\cdot \nabla v$',
                    comp='perp', style='varied',
                    legend_loc='lower left', ylabel=r'Force density N/m^3')

pn.plot_scalar_line(grid, ti.data, flabel=r'T_i', style='varied',
                    legend_loc='lower left', ylabel=r'Temperature eV')
