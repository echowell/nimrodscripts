#!/usr/local/bin/python3
#
# Input files:
#     firstDump - hdf5 dump file with equilibrium and first nk perturbation
#     secondDump - hdf5 dump file with last nk perturbation
# Ouput file:
#     finalDump -hdf5 dump file that has combines the perturbed fourier nmodes
#     from each dump file
#     I assume that the is no overlap in keff
import os
import h5py
from itertools import product
import numpy as np
import sys
homeDir = os.environ['HOME']

#begin user inputs

basePath = homeDir +'/SCRATCH/166439/03300_q104_reorder_combine/vac/'

dumpList = ['n0/dumpglln0.h5','n1/dumpglln1.h5','n2/dumpglln2.h5','n3/dumpglln3.h5','n4/dumpglln4.h5','n5/dumpglln5.h5']

finalDump = basePath+'n0-5/test.h5'
newStep=0
newTime=0.0

#end of inputs
veList = [ 'bq' , 'diff', 'jq', 'nq', 'peq', 'prq', 'rz', 'vq', 'psi_eq']
vvList = [ 'imbe', 'imve', 'rebe', 'reve']

vsList = [ 'imconc', 'imnd', 'impe', 'impr', 'imte', 'imti', \
           'reconc', 'rend', 'repe', 'repr', 'rete', 'reti' ]

fileList = []
nkList = []
kStart=[]
nmodes=0

fc = h5py.File(finalDump, 'w') 

for id, dumpFile in enumerate(dumpList):
    fileList.append(h5py.File(basePath+dumpFile, 'r'))
    kStart.append(nmodes)
    thisK=fileList[id]['keff'].size
    nmodes+=thisK
    nkList.append(thisK)


# reset time and step
fileList[0].copy(fileList[0]['dumpTime'], fc)
fc['dumpTime'].attrs.modify('vsStep',newStep)
fc['dumpTime'].attrs.modify('vsTime',newTime)

newKeff = np.zeros(nmodes)
for id in range(len(fileList)):
    for ii in range(nkList[id]):
        newKeff[kStart[id]+ii]=fileList[id]['keff'][ii]
fc.create_dataset('keff', data=newKeff)
fc.copy(fileList[0]['seams'],fc)
for aname, avalue in fileList[0].attrs.items():
    fc.attrs[aname] = avalue 
fc.attrs['nmodes'] = nmodes

fc.create_group('rblocks')
for aname, avalue in fileList[0]['rblocks'].attrs.items():
    fc['rblocks'].attrs[aname] = avalue 
for re in fileList[0]['rblocks'].keys():
    print('Processing rblock '+re)
    sys.stdout.flush()
    gList=[]
    for id in range(len(fileList)):
        gList.append(fileList[id]['rblocks/'+re])
    gc = fc.create_group('rblocks/'+re)
    for aname, avalue in gList[0].attrs.items():
        gc.attrs[aname] = avalue
    gc.attrs['nfour'] = nmodes

    for d1key, d1value in gList[0].items():
        print('Processing '+d1key)
        sys.stdout.flush()
# copy eq fieds from first dumpfile
        if d1key.startswith(tuple(veList)):
            gc.create_dataset(d1key, data=d1value)
            for aname, avalue in gList[0][d1key].attrs.items():
                gc[d1key].attrs[aname] = avalue
            continue
        if not(d1key.startswith(tuple(vsList+vvList))):
            print("Unreconized key: "+d1key)
            continue
        for (iv,jv) in product(range(d1value.shape[0]),range(d1value.shape[1])):
            if(d1key.startswith(tuple(vsList))): #scalar field
                dcvalue=np.zeros([d1value.shape[0],d1value.shape[1],nmodes])
                for id in range(len(fileList)):
                    dvalue=gList[id][d1key][:]
                    for nn in range(nkList[id]):
                        dcvalue[iv,jv,nn+kStart[id]]=dvalue[iv][jv][nn]        
            else: #vector field
                dcvalue=np.zeros([d1value.shape[0],d1value.shape[1],3*nmodes])
                for id in range(len(fileList)):
                    dvalue=gList[id][d1key][:]
                    for nn in range(nkList[id]):
                        dcvalue[iv,jv,3*(nn+kStart[id])]=dvalue[iv][jv][3*nn]
                        dcvalue[iv,jv,3*(nn+kStart[id])+1]=dvalue[iv][jv][3*nn+1]
                        dcvalue[iv,jv,3*(nn+kStart[id])+2]=dvalue[iv][jv][3*nn+2]
        gc.create_dataset(d1key, data=dcvalue)
        for aname, avalue in gList[0][d1key].attrs.items():
            gc[d1key].attrs[aname] = avalue

