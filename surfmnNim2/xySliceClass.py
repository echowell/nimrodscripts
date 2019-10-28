#!/usr/local/bin/python3

import struct
import numpy as np
import matplotlib.pyplot as plt

class xyClass:
  ''' Base clase for reading ax xy file'''
  def readxySlice(self):
    ''' Read xy_slice.bin '''
    npx=self.mx*self.pd+1
    npy=self.my*self.pd+1
    tf=np.zeros(47,dtype='f',order='F')
    self.ix=np.zeros([npx,npy],dtype='f',order='F')
    self.iy=np.zeros([npx,npy],dtype='f',order='F')
    self.R=np.zeros([npx,npy],dtype='f',order='F')
    self.Z=np.zeros([npx,npy],dtype='f',order='F')
    self.B0R=np.zeros([npx,npy],dtype='f',order='F')
    self.B0Z=np.zeros([npx,npy],dtype='f',order='F')
    self.B0T=np.zeros([npx,npy],dtype='f',order='F')
    self.J0R=np.zeros([npx,npy],dtype='f',order='F')
    self.J0Z=np.zeros([npx,npy],dtype='f',order='F')
    self.J0T=np.zeros([npx,npy],dtype='f',order='F')  
    self.V0R=np.zeros([npx,npy],dtype='f',order='F')
    self.V0Z=np.zeros([npx,npy],dtype='f',order='F')
    self.V0T=np.zeros([npx,npy],dtype='f',order='F')
    self.P0=np.zeros([npx,npy],dtype='f',order='F')
    self.PE0=np.zeros([npx,npy],dtype='f',order='F')
    self.n0=np.zeros([npx,npy],dtype='f',order='F')
    self.diff=np.zeros([npx,npy],dtype='f',order='F')
    self.BRr=np.zeros([npx,npy],dtype='f',order='F')
    self.BZr=np.zeros([npx,npy],dtype='f',order='F')
    self.BTr=np.zeros([npx,npy],dtype='f',order='F')
    self.BRi=np.zeros([npx,npy],dtype='f',order='F')
    self.BZi=np.zeros([npx,npy],dtype='f',order='F')
    self.BTi=np.zeros([npx,npy],dtype='f',order='F')
    self.JRr=np.zeros([npx,npy],dtype='f',order='F')
    self.JZr=np.zeros([npx,npy],dtype='f',order='F')
    self.JTr=np.zeros([npx,npy],dtype='f',order='F')
    self.JRi=np.zeros([npx,npy],dtype='f',order='F')
    self.JZi=np.zeros([npx,npy],dtype='f',order='F')
    self.JTi=np.zeros([npx,npy],dtype='f',order='F')
    self.VRr=np.zeros([npx,npy],dtype='f',order='F')
    self.VZr=np.zeros([npx,npy],dtype='f',order='F')
    self.VTr=np.zeros([npx,npy],dtype='f',order='F')
    self.VRi=np.zeros([npx,npy],dtype='f',order='F')
    self.VZi=np.zeros([npx,npy],dtype='f',order='F')
    self.VTi=np.zeros([npx,npy],dtype='f',order='F')
    self.Pr=np.zeros([npx,npy],dtype='f',order='F')
    self.Pi=np.zeros([npx,npy],dtype='f',order='F')
    self.PEr=np.zeros([npx,npy],dtype='f',order='F')
    self.PEi=np.zeros([npx,npy],dtype='f',order='F')
    self.Nr=np.zeros([npx,npy],dtype='f',order='F')
    self.Ni=np.zeros([npx,npy],dtype='f',order='F')
    self.Cr=np.zeros([npx,npy],dtype='f',order='F')
    self.Ci=np.zeros([npx,npy],dtype='f',order='F')
    self.TEr=np.zeros([npx,npy],dtype='f',order='F')
    self.TEi=np.zeros([npx,npy],dtype='f',order='F')
    self.TIr=np.zeros([npx,npy],dtype='f',order='F')
    self.TIi=np.zeros([npx,npy],dtype='f',order='F')
    self.N=np.zeros([npx,npy],dtype='f',order='F')

    with open(self.file,'rb') as thisFile:
      jj=0
      thisFile.seek(0)

      while jj < npy:
        ii=0
        while ii < npx:
          thisLine=thisFile.read(4)
          blah=struct.unpack(">l",thisLine)
          thisLine=thisFile.read(188)
          tf = struct.unpack(">"+47*'f', thisLine)
          if jj==0 and (ii==0 or ii==1):
            tf1=tf
          thisLine=thisFile.read(4)
          blah=struct.unpack(">l",thisLine)
          self.ix[ii,jj]=tf[0]
          self.iy[ii,jj]=tf[1]
          self.R[ii,jj]=tf[2]
          self.Z[ii,jj]=tf[3]
          self.B0R[ii,jj]=tf[4]
          self.B0Z[ii,jj]=tf[5]
          self.B0T[ii,jj]=tf[6]
          self.J0R[ii,jj]=tf[7]
          self.J0Z[ii,jj]=tf[8]
          self.J0T[ii,jj]=tf[9]
          self.V0R[ii,jj]=tf[10]
          self.V0Z[ii,jj]=tf[11]
          self.V0T[ii,jj]=tf[12]
          self.P0[ii,jj]=tf[13]
          self.PE0[ii,jj]=tf[14]
          self.n0[ii,jj]=tf[15]
          self.diff[ii,jj]=tf[16]
          self.BRr[ii,jj]=tf[17]
          self.BZr[ii,jj]=tf[18]
          self.BTr[ii,jj]=tf[19]
          self.BRi[ii,jj]=tf[20]
          self.BZi[ii,jj]=tf[21]
          self.BTi[ii,jj]=tf[22]
          self.JRr[ii,jj]=tf[23]
          self.JZr[ii,jj]=tf[24]
          self.JTr[ii,jj]=tf[25]
          self.JRi[ii,jj]=tf[26]
          self.JZi[ii,jj]=tf[27]
          self.JTi[ii,jj]=tf[28]
          self.VRr[ii,jj]=tf[29]
          self.VZr[ii,jj]=tf[30]
          self.VTr[ii,jj]=tf[31]
          self.VRi[ii,jj]=tf[32]
          self.VZi[ii,jj]=tf[33]
          self.VTi[ii,jj]=tf[34]
          self.Pr[ii,jj]=tf[35]
          self.Pi[ii,jj]=tf[36]
          self.PEr[ii,jj]=tf[37]
          self.PEi[ii,jj]=tf[38]
          self.Nr[ii,jj]=tf[39]
          self.Ni[ii,jj]=tf[40]
          self.Cr[ii,jj]=tf[41]
          self.Ci[ii,jj]=tf[42]
          self.TEr[ii,jj]=tf[43]
          self.TEi[ii,jj]=tf[44]
          self.TIr[ii,jj]=tf[45]
          self.TIi[ii,jj]=tf[46]
          ii=ii+1
        if (jj<(npy-1)):
          thisLine= thisFile.read(8)
        jj=jj+1

    for i in range(len(self.BRr[:,0])-1):
      for j in self.slicesy:
        self.N[i,j]=-(self.R[i,j]/self.BZr[i,j])*(self.BZr[i+1,j]-self.BZr[i,j])/(self.R[i+1,j]-self.R[i,j])
      if self.BZr[i,j]==0:
        self.N[i,j]=0
      if (self.R[i+1,j]-self.R[i,j])==0:
        self.N[i,j]=0
      self.N[-1,j]=-(self.R[-1,j]/self.BZr[-1,j])*(self.BZr[-1,j]-self.BZr[-2,j])/(self.R[-1,j]-self.R[-2,j])    

  def plotxySlice(self):
    ''' Plot fields in xy_slice.bin. Currently only eq fields are plotted'''
    eqRows=3
    eqCols=4
    eqFigSize=(12,10)
    eqFig, eqAxs = plt.subplots(eqRows,eqCols,figsize=eqFigSize)
    eqAxs[0,0].contourf(self.R,self.Z,self.B0R)
    eqAxs[0,0].set_aspect('equal')
    eqAxs[0,0].set_title(r'$B_r$')

    eqAxs[0,1].contourf(self.R,self.Z,self.B0Z)
    eqAxs[0,1].set_aspect('equal')
    eqAxs[0,1].set_title(r'$B_z$')

    eqAxs[0,2].contourf(self.R,self.Z,self.B0T)
    eqAxs[0,2].set_aspect('equal')
    eqAxs[0,2].set_title(r'$RB_\phi$')

    eqAxs[0,3].contourf(self.R,self.Z,self.P0)
    eqAxs[0,3].set_aspect('equal')
    eqAxs[0,3].set_title(r'$Pr$')

    eqAxs[1,0].contourf(self.R,self.Z,self.J0R)
    eqAxs[1,0].set_aspect('equal')
    eqAxs[1,0].set_title(r'$J_r$')

    eqAxs[1,1].contourf(self.R,self.Z,self.J0Z)
    eqAxs[1,1].set_aspect('equal')
    eqAxs[1,1].set_title(r'$J_z$')

    eqAxs[1,2].contourf(self.R,self.Z,self.J0T)
    eqAxs[1,2].set_aspect('equal')
    eqAxs[1,2].set_title(r'$J_\phi/R$')

    eqAxs[1,3].contourf(self.R,self.Z,self.PE0)
    eqAxs[1,3].set_aspect('equal')
    eqAxs[1,3].set_title(r'$Pr_e$')

    eqAxs[2,0].contourf(self.R,self.Z,self.V0R)
    eqAxs[2,0].set_aspect('equal')
    eqAxs[2,0].set_title(r'$V_r$')

    eqAxs[2,1].contourf(self.R,self.Z,self.V0Z)
    eqAxs[2,1].set_aspect('equal')
    eqAxs[2,1].set_title(r'$V_z$')

    eqAxs[2,2].contourf(self.R,self.Z,self.V0T)
    eqAxs[2,2].set_aspect('equal')
    eqAxs[2,2].set_title(r'$V_\phi$')

    eqAxs[2,3].contourf(self.R,self.Z,self.n0)
    eqAxs[2,3].set_aspect('equal')
    eqAxs[2,3].set_title(r'$nd$')

    plt.show()
  def __init__(self,xyFile,mx,my,pd,plot):
    self.file = xyFile
    self.mx = mx
    self.my = my
    self.pd = pd
    self.intSize=4
    self.floatSize=4
    self.slicesy=[0] #not sure what this does
    self.readxySlice()
    if plot:
      self.plotxySlice()
