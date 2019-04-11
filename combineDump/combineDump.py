#!/usr/local/bin/python3
#
# Input files:
#     firstDump - hdf5 dump file with equilibrium and first nk perturbation
#     secondDump - hdf5 dump file with last nk perturbation
# Ouput file:
#     finalDump -hdf5 dump file that has combines the perturbed fourier nmodes
#     from each dump file
#     I assume that the is no overlap in keff
import h5py
from itertools import product
import numpy as np

veList = [ 'bq' , 'diff', 'jq', 'nq', 'peq', 'prq', 'rz', 'vq', 'psi_eq']
vvList = [ 'imbe', 'imve', 'rebe', 'reve']

vsList = [ 'imconc', 'imnd', 'impe', 'impr', 'imte', 'imti', \
           'reconc', 'rend', 'repe', 'repr', 'rete', 'reti' ]


basePath = '/home/research/ehowell/SCRATCH/166439/03300_q104_reorder_combine/vac/'

firstDump = basePath+'n0-5/dumpglln04.h5'
secondDump = basePath+'n5/dumpglln5.h5'
finalDump = basePath+'n0-5/dumpglln05.h5'

newStep=0
newTime=0.0

f1 = h5py.File(firstDump, 'r') #fe
f2 = h5py.File(secondDump, 'r') #fp
fc = h5py.File(finalDump, 'w') #f0

# reset time and step
f1.copy(f1['dumpTime'], fc)
fc['dumpTime'].attrs.modify('vsStep',newStep)
fc['dumpTime'].attrs.modify('vsTime',newTime)

nk1=f1['keff'].size
nk2=f2['keff'].size
nkc = nk1 + nk2
newKeff = np.zeros(nkc)
for ii in range(nk1):
    newKeff[ii]=f1['keff'][ii]
for ii in range(nk2):
    newKeff[nk1+ii]=f2['keff'][ii]
fc.create_dataset('keff', data=newKeff)

fc.copy(f1['seams'], fc)
# copy file attriubtes and update nmodes
for aname, avalue in f1.attrs.items():
    fc.attrs[aname] = avalue 
fc.attrs['nmodes'] = nkc

fc.create_group('rblocks')
# rblocks has no attributes
#for aname, avalue in f1['rblocks'].attrs.items():
#    print(aname,avalue)
#  fc['rblocks'].attrs[aname] = avalue 
#fc['rblocks'].attrs['nmodes'] = nkc

#loop over rblocks in list
for re in f1['rblocks'].keys():
    print('Processing rblock ' + re)
    g1 = f1['rblocks/'+re]
    g2 = f2['rblocks/'+re]
    gc = fc.create_group('rblocks/'+re)
    for aname, avalue in g1.attrs.items():
        gc.attrs[aname] = avalue
    gc.attrs['nfour'] = nkc

    for d1key, d1value in g1.items():
# copy eq fieds from first dumpfile
        if d1key.startswith(tuple(veList)):
            gc.create_dataset(d1key, data=d1value)
            for aname, avalue in g1[d1key].attrs.items():
                gc[d1key].attrs[aname] = avalue
            continue
        d2value=g2[d1key][:]
        if(d1key.startswith(tuple(vsList))): #scalar field
            dcvalue=np.zeros([d1value.shape[0],d1value.shape[1],nkc])
            for (iv,jv) in product(range(d1value.shape[0]),range(d1value.shape[1])):
                dcvalue[iv,jv,0:nk1-1]=d1value[iv][jv][0:nk1-1]
                dcvalue[iv,jv,nk1:nkc-1]=d2value[iv][jv][nk1:nkc-1]
        else: #vector field
            dcvalue=np.zeros([d1value.shape[0],d1value.shape[1],3*nkc])
            for (iv,jv) in product(range(d1value.shape[0]),range(d1value.shape[1])):
                for nn in range(nk1):
                    dcvalue[iv,jv,3*nn]=d1value[iv][jv][3*nn]
                    dcvalue[iv,jv,3*nn+1]=d1value[iv][jv][3*nn+1]
                    dcvalue[iv,jv,3*nn+2]=d1value[iv][jv][3*nn+2]
                for nn in range(nk1,nkc):
                    dcvalue[iv,jv,3*nn]=d2value[iv][jv][3*(nn-nk1)]
                    dcvalue[iv,jv,3*nn+1]=d2value[iv][jv][3*(nn-nk1)+1]
                    dcvalue[iv,jv,3*nn+2]=d2value[iv][jv][3*(nn-nk1)+2]
        gc.create_dataset(d1key, data=dcvalue)
        for aname, avalue in g1[d1key].attrs.items():
            gc[d1key].attrs[aname] = avalue
        print(d1key, d1value.shape, d2value.shape, dcvalue.shape)
#        print(d1key,d1value.shape,d2value.shape)
#        for (iv,jv) in product(range(len(d1value)),range(len(d1value[0]))):
#            if d1key.startswith(tuple(vsList)):
#                print(d1key)
'''
  for de in ge.iteritems() :
    dse = de[1][:]
    if de[0].startswith('psi_eq') : 
      continue
    dsp = gp[de[0]][:]
    for (iv,jv) in product(range(len(dse)),range(len(dse[0]))) :
      if de[0].startswith(tuple(vsList)) :
	dse[iv][jv] = 0. * len(dse[iv][jv])
	for n in nList :
          dse[iv][jv][n] = dsp[iv][jv][nList[n]] 
      elif de[0].startswith(tuple(vvList)) :
	#if de[0]=='reve0001' :
	#  print de[0], iv, jv, len(dse[iv][jv]), dse[iv][jv][3*0:3*0+2], dsp[iv][jv][3*nList[0]:3*nList[0]+2]
	dse[iv][jv] = 0. * len(dse[iv][jv])
	for n in nList :
	  dse[iv][jv][3*n]   = dsp[iv][jv][3*nList[n]]
	  dse[iv][jv][3*n+1] = dsp[iv][jv][3*nList[n]+1]
	  dse[iv][jv][3*n+2] = dsp[iv][jv][3*nList[n]+2]
      
    g0.create_dataset(de[0], data=dse)
    for aname, avalue in ge[de[0]].attrs.items():
      g0[de[0]].attrs[aname] = avalue
'''
     