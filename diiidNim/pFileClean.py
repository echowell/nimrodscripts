#!/usr/local/bin/python3
''' This script is useful for cleaning up the pfile'''

###################################
import os
import numpy as np

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
    fitPsi = -1.0
    fieldPsi = -1.0
    name = ""
    def __init__ (self, name, field0, fitPsi, fieldPsi):
        self.name = name
        self.field0 = field0
        self.fitPsi = fitPsi
        self.fieldPsi = fieldPsi

filePath = homeDir + "/SCRATCH/rmp_wkdir/166439_03300/"
pFileName = "p166439.03300"
pFileWrite = pFileName + ".smooth"

# Set up field 
omgeb = modField("omgeb(kRad/s)",0.0,0.162306,55955.2)

fixList = [omgeb]

print(filePath + pFileName)
thisFile = open(filePath+pFileName,"r")
writeFile = open(filePath+pFileWrite, "w")

print ("Reading file " + filePath+pFileName )
while True:
    thisLine = thisFile.readline()
    if len(thisLine)==1:
        break
    writeFile.write(thisLine)
    thisWord = thisLine.split()
    print ("Reading Field " + thisWord[2])
    thisField = readField(thisFile,int(thisWord[0]))
    for iFix in fixList:
        if iFix.name == thisWord[2]:
            print(iFix.name)
    print ("Writing Field " + thisWord[2])
    writeField(writeFile,thisField)



thisFile.close()
writeFile.close()
