#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
import h5py


homeDir = os.environ['HOME']
scratchDir = homeDir + '/SCRATCH'
fileName = scratchDir+'/174446_novac_fl/eq26/1908013/05000/surfmn.05000.h5'
#fileName = scratchDir+'/174446_novac_debug/vac_eq27_rmp/surfmn.00050.h5'
#fileName = scratchDir+'/166439//03300_vac_eq/complexconj_rmp_vac/surfmn.00005.h5'


profNames = ['Vprime','q']
with h5py.File(fileName,'r') as fc:
  for aname, avalue in fc.attrs.items():
    print(aname,avalue)
  mrGrid = fc['surfmnGrid'][:]
  bmn = fc['Bmn001'][:]
  rho = fc['rho'][:]
  profs = fc['prof'][:]
  print(fc.keys())
  

print(mrGrid.shape)
print(bmn.shape)
plt.set_cmap('nipy_spectral')
plt.contourf(mrGrid[0,:,:],mrGrid[1,:,:],bmn,levels=300)
plt.colorbar()
plt.show()

plt.plot(rho,profs[0,:],label=profNames[0])
plt.legend()
plt.show()

plt.plot(rho,profs[1,:],label=profNames[1])
plt.legend()
plt.show()

#plt.plot(mrGrid[1,:,1],bmn[:,0], label='m=-5')
plt.plot(mrGrid[1,:,1],bmn[:,1], label='m=-4')
plt.plot(mrGrid[1,:,1],bmn[:,2], label='m=-3')
plt.plot(mrGrid[1,:,1],bmn[:,3], label='m=-2')
plt.plot(mrGrid[1,:,1],bmn[:,4], label='m=-1')
#plt.plot(mrGrid[1,:,1],bmn[:,5], label='m=0')
plt.plot(mrGrid[1,:,1],bmn[:,6], label='m=1')
plt.plot(mrGrid[1,:,1],bmn[:,7], label='m=2')
plt.plot(mrGrid[1,:,1],bmn[:,8], label='m=3')
plt.plot(mrGrid[1,:,1],bmn[:,9], label='m=4')
#plt.plot(mrGrid[1,:,1],bmn[:,10], label='m=5')
plt.axvline(0.537,ls=':',c='k')
plt.axvline(0.736,ls=':',c='k')
plt.axvline(0.848,ls=':',c='k')
plt.legend(loc=0)
plt.show()