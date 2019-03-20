#!/usr/local/bin/python3
''' Base class for storing and reading trip3D data '''
import numpy as np
class TripClass:
    probeRzFile = ""
    probeBFile = ""
    probeAFile = ""
    bData = np.zeros(1)
    aData = np.zeros(1)
    npoints = 0
    nphi = 0
    rr = np.zeros(1)
    zz = np.zeros(1)
    phi = np.zeros(1)
    brReal = np.zeros(1)
    bzReal = np.zeros(1)
    btReal = np.zeros(1)
    brPhase = np.zeros(1,dtype=np.complex_)
    bzPhase = np.zeros(1,dtype=np.complex_)
    shiftindex=0
    def __init__ (self,rzFile, aFile, bFile, shiftindex):
        ''' Initialize TripClass '''
        self.probeRzFile = rzFile
        self.probeBFile = bFile
        self.probeAFile = aFile
        self.shiftindex = shiftindex
    def readRz (self):
        ''' Read probeRz File as Save data '''
        #todo
    def readBFile (self):
        ''' Read probeBFile as Save Data '''
        #b data is stores phi, r, z, B_phi, B_R, B_Z, B_p, B_mag, PsiN_pol
        self.bData = np.loadtxt(self.probeBFile,comments='%')
    def reorderBFile(self):
        ''' Reorder the bfile to account for a shift in the data '''
        if (self.shiftindex==0):
            return
        else:
            # reorder data becuase jake changed the indexing in NIMROD
            tempData = self.bData
            for ip in range(self.nphi):
                startIndex= ip * self.npoints
                for ii in range(self.npoints):
                    if(ii>=self.shiftindex):
                        self.bData[ii-self.shiftindex+startIndex,:]=tempData[ii+startIndex,:]
                    else:
                        self.bData[self.npoints-self.shiftindex+ii+startIndex,:]=tempData[ii+1+startIndex,:]
#            self.bData = np.zeros(tempData.shape)
#            maxData=self.bData.shape[0]
#            for ii in range(maxData):                
    def readAFile (self):
        ''' Read probeAFile as Save Data '''
        #b data is stores phi, r, z, A_phi, A_R, A_Z, A_p, A_mag, PsiN_pol
        self.aData = np.loadtxt(self.probeAFile,comments='%')
    def findNumPoints (self, data):
        phi0 = data[0,0]
        for ii in range(data.shape[0]):
            if (data[ii,0]!=phi0):
                self.npoints =ii 
                break
        self.nphi = int(data.shape[0]/self.npoints)
    def flipPhi (self):
        for ii, iphi in enumerate(self.phi):
            print (ii, iphi)
            self.phi[ii] = 360.0 - iphi
            if (self.phi[ii]==360.0): self.phi[ii]=0.0
    def processBFile (self):
        self.readBFile()
        self.findNumPoints(self.bData)
        self.reorderBFile()
        self.phi = np.zeros(self.nphi)
        self.brReal = np.zeros([self.npoints,self.nphi])
        self.bzReal = np.zeros([self.npoints,self.nphi])
        self.btReal = np.zeros([self.npoints,self.nphi])
        self.rr = self.bData[:self.npoints,1]
        self.zz = self.bData[:self.npoints,2]
        for ii in range(self.nphi):
            startIndex = ii * self.npoints
            self.phi[ii] = self.bData[startIndex,0]
            self.brReal[:,ii] = self.bData[startIndex:startIndex+self.npoints,4]
            self.bzReal[:,ii] = self.bData[startIndex:startIndex+self.npoints,5]
            self.btReal[:,ii] = self.bData[startIndex:startIndex+self.npoints,3]
        # i don't know if I need to flip phi or change the sign of Bphi
        # self.flipPhi()
        # to be consistant with nimrod I should use the forward fft when going
        # from real space to fourier space, and I also need to devide by nphi
        # numpy does not have an option for normalizing the FFT by n,
        # it can only do sqrt(N) normalization or normalize the IFFT by N.
        self.brPhase = np.fft.fft(self.brReal,axis=1)/(float(self.nphi))
        self.bzPhase = np.fft.fft(self.bzReal,axis=1)/(float(self.nphi))
        self.btPhase = np.fft.fft(self.btReal,axis=1)/(float(self.nphi))
    def writeNimrodBext(self,path,baseFileName,fileExt):
        if (self.nphi % 2 == 0): #even
            maxnphi = int(self.nphi/2)
        else: #odd
            maxnphi = int((self.nphi+1)/2)
        for ii in range (maxnphi +1):
            
            if ii==maxnphi:
                fac=0.5
            else:
                fac=1.0
            print(ii, maxnphi, fac)
            tempFileName = path + baseFileName +"{0:0=2d}".format(ii)  + fileExt
            thisFile = open(tempFileName,'w')
            for jj in range(self.brPhase.shape[0]):
                thisLine = '{: 16.16e}'.format(fac*self.brPhase[jj,ii].real) + ", " 
                thisLine+= '{: 16.16e}'.format(fac*self.brPhase[jj,ii].imag) + ", "
                thisLine+= '{: 16.16e}'.format(fac*self.bzPhase[jj,ii].real) + ", " 
                thisLine+= '{: 16.16e}'.format(fac*self.bzPhase[jj,ii].imag) + ", "
                thisLine+= '{: 16.16e}'.format(fac*self.btPhase[jj,ii].real) + ", " 
                thisLine+= '{: 16.16e}'.format(fac*self.btPhase[jj,ii].imag) + "\n"
                thisFile.write(thisLine)
            thisFile.close()