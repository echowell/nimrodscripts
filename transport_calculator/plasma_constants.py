#!/usr/local/bin/python3

import numpy as np
me = 9.1093898e-31
md = 3.3435860e-27
qe = 1.60217733e-19
mu0 = 4 * np.pi * 1e-7
gamma = 5./3.
kboltz = 1.60217733e-19
clight = 2.99792458e8
eps0 = 8.854187817e-12
 
def coulombLog(te, ne):
#  logLam = 18.4 - 1.15 * np.log10(ne)+2.3*np.log10(te)
# Jim Callen uses 17
  logLam = 24 - np.log(np.sqrt(ne*1.e-6)/te)
  return logLam

def nu_e(ne,te,zeff,logLam):
  nu = np.sqrt(2*np.pi)* ne * zeff * qe**4 * logLam
  nu *= 1.0/(12.0*(np.pi * eps0)**2 * np.sqrt(me) * np.power(kboltz*te,1.5))
  return nu

def nu_i(ni,ti,mi,zstar,logLam):
  nu = np.sqrt(np.pi)* ni* (1.0+np.sqrt(2.0)*zstar)  * qe**4 * logLam
  nu *= 1.0/(12.0*(np.pi * eps0)**2 * np.sqrt(mi) * np.power(kboltz*ti,1.5))
  return nu

def vte(te):
  return np.sqrt(2*kboltz*te/me)

def vts(ts,ms):
  return np.sqrt(2*kboltz*ts/ms)

def zStar(ne,ni,zeff):
  return (ne/ni)*zeff-1.0

def circulatingFraction(eps):
  return 1-1.46*np.sqrt(eps) + 0.46 * np.power(eps,1.5)

def nuBanana(nu,eps):
  return np.sqrt(eps)*nu

def nuStar(R0,q0,epsilon,lam):
  return R0*q0/(np.power(epsilon,1.5)*lam)

def neoMue(nue,nueStar,zeff,eps):
  numer = 1.46 * np.sqrt(eps) * (1.+0.533/zeff)*nue
  denomFac = (2.4*zeff**2+5.32*zeff+2.225)/(4.25*zeff**2+3.02*zeff)
  denom = (1.+np.sqrt(nueStar)+1.65*(1.+0.533/zeff)*nueStar)
  denom *= (1+1.18*denomFac*np.power(eps,1.5)*nueStar)
  return numer/denom

def neoMui(nu,nuStar,zStar,eps):
  numer = 1.46 * np.sqrt(eps) * (zStar+0.533)/(zStar+0.707)*nu
  denomFac = (2.4*zStar**2+5.32*zStar+2.225)/((zStar+0.707)*(4.25*zStar+3.02))
  denom = (1.+np.sqrt(nuStar)+1.65*(zStar+0.533)/(zStar+0.707)*nuStar)
  denom *= (1+1.18*denomFac*np.power(eps,1.5)*nuStar)
  return numer/denom