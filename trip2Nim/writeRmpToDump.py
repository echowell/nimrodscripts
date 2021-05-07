#!/usr/bin/env python3

# This script writes rmp data stored in brmpn##.dat file to a brmp fields in
# the NIMROD dumpfile.
#
# Input files:
#     dumpFile - nimrod hdf5 dump file, this file will be modified
#     brmp##.dat - dat files storing the rmp fields at the nimrod boudnary
# Ouput file:
#     dumpFile -hdf5 dump file

#to do list...
#script inputs -dump file, brmpfile directory, flag to add rmp to be, optional n-list

#read brmp files
#check for keff capibility
#
import h5py
import numpy as np

def getBlockKey(bl):
    '''Return block key string corresponding to the block number bl'''
    if not (isinstance(bl,(int,np.integer)) ):
        raise TypeError
    elif (bl>=10000 or bl<1):
        print("bl must be in [1,9999]")
        raise ValueError
    return f'{bl:04d}'

def getSeamKey(bl):
    '''Return seam key string corresponding to the block number bl'''
    if not (isinstance(bl,(int,np.integer)) ):
        raise TypeError
    elif (bl>=10000 or bl<0):
        print("bl must be in [0,9999]")
        raise ValueError
    return f'{bl:04d}'

def getBeKeys(bl):
    '''Return Be keys strings corresponding to the block number bl'''
    if not (isinstance(bl,(int,np.integer)) ):
        raise TypeError
    elif (bl>=10000 or bl<1):
        print("bl must be in [1,9999]")
        raise ValueError
    reKey='rebe'+f'{bl:04d}'
    imKey='imbe'+f'{bl:04d}'
    return reKey,imKey

def getBrmpKeys(bl):
    '''Return Brmp keys strings corresponding to the block number bl'''
    if not (isinstance(bl,(int,np.integer)) ):
        raise TypeError
    elif (bl>=10000 or bl<1):
        print("bl must be in [1,9999]")
        raise ValueError
    reKey='rebrmp'+f'{bl:04d}'
    imKey='imbrmp'+f'{bl:04d}'
    return reKey,imKey

def getRzedgeKey(bl):
    '''Return rzedge key string corresponding to the block number bl'''
    if not (isinstance(bl,(int,np.integer)) ):
        raise TypeError
    elif (bl>=10000 or bl<1):
        print("bl must be in [1,9999]")
        raise ValueError
    reKey='rzedge'+f'{bl:04d}'
    return reKey

def readBrmp(nn,fdir='',prefix='brmpn',suffix='.dat'):
    '''reads brmp file and returns an brmp as a np array'''
    if not (isinstance(nn,(int,np.integer)) ):
        raise TypeError
    elif (nn<0):
        print("nn must be zero or positive")
        raise ValueError
    brmpFile=fdir+prefix+f'{nn:02d}'+suffix
    try:
        tempData=np.genfromtxt(brmpFile,delimiter=',')
        #second index of tempData should have 6 values
        #Each corresponding or Re or Im of each vector component of Bn
        if not (tempData.shape[1]==6):
            print(f'Data in {brmpFile} has wrong dimensions')
            raise ValueError
        brmp = np.zeros([tempData.shape[0],3],dtype=np.complex128)
        for ii in range(tempData.shape[0]):
            brmp[ii,0]=tempData[ii,0]+tempData[ii,1]*1.0j
            brmp[ii,1]=tempData[ii,2]+tempData[ii,3]*1.0j
            brmp[ii,2]=tempData[ii,4]+tempData[ii,5]*1.0j
    except OSError:
        print(f'{brmpFile} not found')
        raise
    return brmp


def writeBrmp(dumpFile,nList=[],hardFail=False):
    '''Main driver function to write Brmp to dump file'''
    with h5py.File(dumpFile, 'r+') as h5DumpFile:
        #get list of nmodes
        nmodeDict = {}
        fillList=False
        if not nList:
            fillList=True
        for ii, nn in enumerate(h5DumpFile['keff']):
            nmodeDict[int(nn)]=ii
            if fillList:
                nList.append(int(nn))
        #analyze mesh need to check if mesh is a polar mesh and find nybl
        polarMesh=True
        nybl=0
        seam0=h5DumpFile['seams'][getSeamKey(0)]
        edgeBlock=set()
        for vertex in seam0['vertex']:
            edgeBlock.add(vertex[0])
        nybl=len(edgeBlock)
        #search for corners, if found we know mesh is not polar
        for excorner in seam0['excorner']:
            if (excorner):
                polarMesh=False
                break
        if not polarMesh:
            raise Exception("Mesh in dumpfile is not a polar mesh.")
        #loop over nList and modify Brmp
        for nn in nList:
            if not nn in nmodeDict:
                print(f'Warning {nn} is not a Fourier mode in {dumpFile}')
                if hardFail:
                    raise ValueError
                else:
                    continue
            try:
                brmp=readBrmp(nn)
                nIndex=nmodeDict[nn]
                print(nIndex)
            except OSError:
                if hardFail:
                    raise
                else:
                    print("Continuing with next foruier mode")
            except:
                raise
        print(nmodeDict)
        print(nList)

if __name__ == "__main__":
    dumpFile = "/home/ehowell/SCRATCH/166439/03300_fgnimeq_q104_reorder_normp/dumpgll.00000.h5"
    writeBrmp(dumpFile,nList=[15],hardFail=True)




#exit()
dumpFile = "/home/ehowell/SCRATCH/166439/03300_fgnimeq_q104_reorder_test/dumpgll.00000.h5"
h5DumpFile = h5py.File(dumpFile, 'r') #h5instance of dumpfile r+ for read/write

#file structure
# h5dump.keys = dumpTime, keff, rblocks, seams
for key in h5DumpFile.attrs.items():
    print(key)
for key in h5DumpFile['dumpTime'].attrs:
    print(key)
print(h5DumpFile['dumpTime'].attrs.get('vsTime'))
print(h5DumpFile['dumpTime'])
#get list of nmodes
nmodeList = []
for ii in h5DumpFile['keff']:
    nmodeList.append(int(ii))

#get list of external vertices from seam0
#seam zero has 3 keys, np, vertex, and excorner
#vertex has a list of block and vertex id
#vertex indicate how many vertex share an exterior vertex
#excorner is a flag to indicate vertex is a corner
seam0=h5DumpFile['seams'][getSeamKey(0)]
h5np=list(seam0['np'])
print(h5np)
edgeBlock=set()
for key in seam0:
    print(key)
for vertex in seam0['vertex']:
    edgeBlock.add(vertex[0])
print(edgeBlock)
print(len(edgeBlock))
polar=True
for excorner in seam0['excorner']:
    if (excorner):
        polar=False
        break
    print(excorner)
print(polar)
seam1=h5DumpFile['seams'][getSeamKey(32)]
for key in seam1:
    print(key)
for ii in seam1['intxy']:
    print(ii)
exit()
vertex=list(seam0['vertex'])
print(vertex)
print(vertex[0][0])
sum=0
for block, node in vertex:
    sum+=1
    print(block,node)
    print(h5DumpFile['rblocks'][getBlockKey(block)][getBrmpKeys(block)[0]][node])
    print(h5DumpFile['rblocks'][getBlockKey(block)][getBrmpKeys(block)[1]][node])
print(sum)
exit()

for key in h5DumpFile['seams']['0032']:
    print(key)
for ii in h5DumpFile['seams']['0032']['intxy']:
    print(ii)
for key in h5DumpFile['rblocks']['0032']:
    print(key)
print(h5DumpFile['rblocks']['0032']['rz0032'])
print(getBrmpKeys(32))
print(h5DumpFile['rblocks']['0032'][getRzedgeKey(32)])
for ii in h5DumpFile['rblocks']['0032'][getBrmpKeys(32)[0]]:
    print(ii)
exit()
for key in  h5DumpFile.keys():
    print(key)

for keff in h5DumpFile['keff']:
    print(keff)

for rblock in h5DumpFile['rblocks']:
    print(rblock)

for seams in h5DumpFile['seams']:
    print(seams)

seam0=h5DumpFile['seams']['0000']
for key in seam0:
    print(key)

sum=0
for ii in seam0['np']:
    print(ii)
    sum+=ii

print(sum)
for ii in seam0['vertex']:
    print(ii)
