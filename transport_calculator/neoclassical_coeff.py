#!/usr/local/bin/python3

import numpy as np
import plasma_constants as pc


### set background plasma parameters
ne = 6.82e19 #m^-3
te = 1504. #eV
ti = 1754. #eV
zeff = 3.25
ni = 3.49e19 #m^-3 rough guess

epsilon = 0.25 # magnetic inverse aspect ratio
R0 = 1.696
q0 = 2.0

#calcularte classical values for electrons and ions with impurities
#here I assume the the Coulomb Log is the same for both electron and 
#ions
logLam = pc.coulombLog(te,ne)
nu_e = pc.nu_e(ne,te,zeff,logLam)
vte = pc.vte(te)
lambda_e = vte/nu_e
zStar=pc.zStar(ne,ni,zeff)
nu_i = pc.nu_i(ni,ti,pc.md,zStar,logLam)
vti = pc.vts(ti,pc.md) 
lambda_i = vti/nu_i

fc = pc.circulatingFraction(epsilon)#circuling fraction
nueBanana = pc.nuBanana(nu_e,epsilon)
nuiBanana = pc.nuBanana(nu_i,epsilon)

nueStar = pc.nuStar(R0,q0,epsilon,lambda_e)
nuiStar = pc.nuStar(R0,q0,epsilon,lambda_i)
mu_e = pc.neoMue(nu_e,nueStar,zeff,epsilon)
mu_i = pc.neoMue(nu_i,nuiStar,zStar,epsilon)
print(logLam)
print(nu_e)
print(nu_i)
print (lambda_e)
print (lambda_i)
print(nu_i/nu_e)
print(fc)
print(nueStar)
print(nuiStar)
print(mu_e)
print(mu_i)