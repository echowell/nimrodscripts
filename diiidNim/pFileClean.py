#!/usr/local/bin/python3
''' This script is useful for cleaning up the pfile'''

###################################
import os
import numpy as np
from scipy.interpolate import splev,splrep
import matplotlib.pyplot as plt
homeDir = os.environ['HOME']


def readField(thisFile,npts):
    field = np.zeros([npts,3])
    for ii in range(npts):
        thisLine = thisFile.readline()
        thisWord = thisLine.split()
        for jj in range(3):
            field[ii,jj]=float(thisWord[jj])
    return field

def writeField(writeFile, field):
    for ii in range(field.shape[0]):
        thisline = writeFile.write('{:-f} {:-f} {:-f}\n'.format(field[ii,0],field[ii,1],field[ii,2]))

class modField:
    field0 = 0.0
    dFeild0 = 0.0
    weldPsi = -1.0
    mode = 0
    quadA = 0.0
    quadB = 0.0
    quadC = 0.0
    name = ""
    def __init__ (self, name, field0, dField0, weldPsi, mode):
      self.name = name
      self.field0 = field0
      self.dField0 = dField0
      self.weldPsi = weldPsi
      self.mode = mode
    def quadFit(self,fWeld,dfWeld):
      ''' Assume that the function is quadratic near the axis'''
      self.quadC = self.field0
      if (self.mode==0): # fit derivative at 0
        self.quadB = self.dField0
        self.quadA = (fWeld-self.quadB * self.weldPsi - self.quadC)/(self.weldPsi**2)
      else: #fit derivative at xweld
        self.quadA = (dfWeld * self.weldPsi +self.quadC - fWeld)/self.weldPsi**2
        self.quadB = dfWeld - 2.0* self.quadA * self.weldPsi
    def quadEval(self, psi):
      return self.quadA * psi**2 + self.quadB * psi + self.quadC
    def smooth (self, field):
      ''' This function applies a bump function fit to smooth the field'''
      splineField = splrep(field[:,0],field[:,1],k=3)
      fieldAtWeld = splev(self.weldPsi,splineField)
      dFieldAtWeld = splev(self.weldPsi,splineField,der=1)
      self.quadFit(fieldAtWeld,dFieldAtWeld)
      tempField = np.zeros(field.shape)
      for ix, ipsi in enumerate(field[:,0]):
        tempField[ix,0] = ipsi
        if (ipsi < self.weldPsi):
          tempField[ix,1] = self.quadEval(ipsi)
        else:
          tempField[ix,1] = field[ix,1]

      newSplineField = splrep(tempField[:,0],tempField[:,1],k=3)
      for ix, ipsi in enumerate(tempField[:,0]):
        tempField[ix,2]= splev(ipsi,newSplineField,der=1)
# plot fields
      x2 = np.linspace(0, 1, 200)
      y2 = splev(x2, splineField)
      y3 = splev(x2, newSplineField)
      plt.plot(x2, y2, x2, y3)
      plt.show()
    
      return tempField


filePath = homeDir + "/SCRATCH/174446_debug2/eq1/"
pFileName = "p174446.3390.0_new_rot_fits"
pFileWrite = pFileName + ".smooth"


# Set up field 
omgeb = modField("omgeb(kRad/s)",21.6,0.0,0.1,1)
kpol = modField("kpol(km/s/T)",-7.5,0.0,0.1,1)

fixList = [omgeb,kpol]

print(filePath + pFileName)
thisFile = open(filePath+pFileName,"r")
writeFile = open(filePath+pFileWrite, "w")

print ("Reading file " + filePath+pFileName )
while True:
    thisLine = thisFile.readline()
    if len(thisLine)==1:
        break
    thisWord = thisLine.split()
    if len(thisWord)>4: break
    writeFile.write(thisLine)
    print ("Reading Field " + thisWord[2])
    thisField = readField(thisFile,int(thisWord[0]))
    for iFix in fixList:
        if iFix.name == thisWord[2]:
            print(iFix.name)
            thisField=iFix.smooth(thisField)
    print ("Writing Field " + thisWord[2])
    writeField(writeFile,thisField)


thisFile.close()
writeFile.close()
