#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

homeDir = os.environ['HOME']
scratchPath= "/SCRATCH/174446_novac_debug/eq28/"
scratchPath="/SCRATCH/166439/03300_fgnimeq_q104_reorder_temp/"
fileName = "eqfsa.h5"

fullFile = homeDir + scratchPath + fileName

fontSize=16
eqNames = ['q',r'$\lambda$','beq2','nd','ti','te','pi','pe']
eqNames += [r'$\omega$','kpol']
psiName = [r'$\psi$']
fc=h5py.File(fullFile, 'r')

print(fc['prof'].shape)
data=np.zeros(fc['prof'].shape)
for ii in range(fc['prof'].shape[0]):
  for jj in range(fc['prof'].shape[1]):
    data[ii,jj] =fc['prof'][ii,jj]

psi=np.zeros(fc['psi'].shape)
for ii in range(fc['psi'].shape[0]):
  psi[ii]=fc['psi'][ii]

fig=plt.figure(figsize=(6,8))

#fig.subplots_adjust(left=0.5)
#ax.yaxis.labelpad=35
ax = fig.add_subplot(311)
ax.plot(psi,abs(data[0,:]))
ax.hlines(1.0,0,1,linestyle='dotted')
ax.axvline(0.637122,ls='-',c='k')
plt.xlim(0.0,1.0)
plt.ylabel(eqNames[0],fontsize=fontSize, rotation=90)
plt.title("Saftey Factor",fontsize=fontSize)
ax = fig.add_subplot(312)
ax.plot(psi,abs(data[1,:]))
plt.xlim(0.0,1.0)
ax.axvline(0.637122,ls='-',c='k')
plt.ylabel(eqNames[1]+r' $(m^{-1})$',fontsize=fontSize, rotation=90)
plt.title("Parallel Current",fontsize=fontSize)
ax = fig.add_subplot(313)
ax.plot(psi,abs(data[8,:]/1000))
plt.xlim(0.0,1.0)
ax.axvline(0.637122,ls='-',c='k')
plt.xlabel(psiName[0],fontsize=fontSize)
plt.ylabel(eqNames[8]+r' (krad/s)',fontsize=fontSize, rotation=90)
plt.title("Rotation Profile",fontsize=fontSize)
ax.hlines(0.0,0,1,linestyle='solid')
plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(4,4))

#fig.subplots_adjust(left=0.5)
#ax.yaxis.labelpad=35
ax = fig.add_subplot(111)
ax.plot(psi,data[8,:]/1000)
plt.xlim(0.0,1.05)
#ax.axvline(0.637122,ls='-',c='k')
plt.xlabel(psiName[0],fontsize=fontSize)
plt.ylabel(eqNames[8]+r' (kiloR/s)',fontsize=fontSize, rotation=90)
plt.title("Rotation Profile",fontsize=fontSize)
ax.hlines(0.0,0,1.05,linestyle='solid')
plt.tight_layout()
plt.show()
print(eqNames)
