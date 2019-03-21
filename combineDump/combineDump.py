#!/usr/bin/python 
#
# Input files:
#     dumpeq.h5 - hdf5 dump file with equilibrium
#     dumppert.h5 - hdf5 dump file with perturbations
# Ouput file:
#     dumpgll.00000.h5 -hdf5 dumpe fle that will use the 
#     equilibrium profiles form dumpeq.h5 and perturbations
#     dumppert.h5. The dumpgll.00000.h5 will have the same
#     keff and the same zperiod as dumpeq.h5.
#
# nList is a dictionary to set a correspondance between
# records for different mode numbers in dumpeq.h5 and 
# dumppert.h5. The key values are for the records of 
# different mode numbers in dumpeq.h5. The values are
# for the records in dumppert.h5 that will be used to
# set the perturbations in dumpgll.00000.h5.
#
import h5py
from itertools import product
import numpy as np

vvList = [ 'imbe', 'imve', 'rebe', 'reve' ]

vsList = [ 'imconc', 'imnd', 'impe', 'impr', 'imte', 'imti', \
           'reconc', 'rend', 'repe', 'repr', 'rete', 'reti' ]

nList = {0: 0, 1: 11, 2: 22, 3 : 33 }

fe = h5py.File('dumpeq.h5', 'r')
fp = h5py.File('dumppert.h5', 'r')
f0 = h5py.File('dumpgll.00000.h5', 'w')

fe.copy(fe['dumpTime'], f0)
f0.create_dataset('keff', data=fe['keff'][:])
fe.copy(fe['seams'], f0)
for aname, avalue in fe.attrs.items():
 f0.attrs[aname] = avalue 

f0.create_group('rblocks')
for aname, avalue in fe['rblocks'].attrs.items():
  f0['rblocks'].attrs[aname] = avalue 
for re in fe['rblocks'].iteritems() :
  ge = fe['rblocks/'+re[0]]
  gp = fp['rblocks/'+re[0]]
  g0 = f0.create_group('rblocks/'+re[0])
  for aname, avalue in ge.attrs.items():
    g0.attrs[aname] = avalue
  print('Processing rblock ' + re[0])
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
 
     