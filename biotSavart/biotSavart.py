#!/usr/local/bin/python3
#
# Input files:
# Ouput file:


import os
import numpy as np
import sys
import coilClass as cc
import integrateBS as bs
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

homeDir = os.environ['HOME']
filePath = homeDir + '/SCRATCH/174446_novac_new_eq/surfmn_scan/case_15/'
filePath = homeDir+ '/SCRATCH/174446_novac_new_eq/192601_eq_valgrind/'
rzfile = filePath + 'nimrod_bdry_rz.txt'
baseFileName = "brmpn"
fileExt = ".dat"


phiPlanes = 4 #n=0->2
segmentsPerCoil = 50
baseCurrent = 1.0
nPert = 1


########### begin code #########
#nimrod RZ node locations
nodeRZ = np.loadtxt(rzfile,comments='%',delimiter=',', skiprows=1)
# allocate xyz and b vectors
print(nodeRZ.shape)

nodeXYZ = np.zeros([3,nodeRZ.shape[0],phiPlanes])
bRZPhi = np.zeros([3,nodeRZ.shape[0],phiPlanes])
bRPhase = np.zeros(1,dtype=np.complex_)
bZPhase = np.zeros(1,dtype=np.complex_)
bPhiPhase = np.zeros(1,dtype=np.complex_)
#convert node locations to xyz coordinates at multiple phi planes
for iPhi in range(phiPlanes):
    sinPhi = np.sin(iPhi*2.0*np.pi/phiPlanes)
    cosPhi = np.cos(iPhi*2.0*np.pi/phiPlanes)
    nodeXYZ[0,:,iPhi]=nodeRZ[:,0]*cosPhi
    nodeXYZ[1,:,iPhi]=nodeRZ[:,0]*sinPhi
    nodeXYZ[2,:,iPhi]=nodeRZ[:,1]

# set up c coils
# in the future I can move this to an function
coilList=[]
#for ii in range(6):
#    thisCurrent = baseCurrent * np.cos(np.pi*(ii)*nPert/3.0)
#    thisCoil = cc.coil(thisCurrent,segmentsPerCoil)
#    coilList.append(thisCoil)
#    coilList[ii].cCoil(ii+1)

# coil rmp 2 or 3
distance_coil = 1.2
coil_theta = [0, .35, .5, .65]
coil_tilt = [0.25, .35, .75, 1.1]
coil_tilt = [0.25, .40, .75, 1.1]
# coils for eq 28
coil_theta = [0, .15, .5, .85]
coil_tilt = [0.25, .65, .75, 0.85]
#delta theta
delta_theta=0.20
delta_tilt=.17

coil_theta = [0, delta_theta, .5, 1-delta_theta]
coil_tilt = [0.25, .5+delta_tilt, .75, 1-delta_tilt]
scope=False
#scope=True

rmag_axis=1.7682
zmag_axis=0.0
coil_r = 0.6
for ii in range(6):
  phi = np.pi * ii/3.0
  thisCurrent = baseCurrent * np.cos(phi*nPert)
  for jj in range(4): 
    theta = coil_theta[jj] * 2.0 * np.pi
#    theta = np.pi * jj/2.0
    if (jj%2 == 0):
      this_m_current = -thisCurrent 
    else:
      this_m_current= thisCurrent

    r0 = rmag_axis + distance_coil * np.cos(theta)
    z0 = zmag_axis + distance_coil * np.sin(theta)
    x0 = r0 * np.cos(phi)
    y0 = r0 * np.sin(phi)
    tx = 0.0 
#    ty = (jj+1) * np.pi/2.0
    ty = coil_tilt[jj] * 2.0 * np.pi
    tz = phi

    # temp set rotation angles to zero to debug
#    tx=0.0
#    ty=0.0
#    tz=0.0

    thisCoil = cc.coil(this_m_current,segmentsPerCoil)
    thisCoil.planarCoil(x0,y0,z0,coil_r,tx,ty,tz)
    coilList.append(thisCoil)

# plot coils
fig =plt.figure()
ax = fig.add_subplot(111,projection='3d')
for iCoil in coilList:
#  print(iCoil.xyz[2,:])
  iCoil.plot_coil(ax)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()

rd3d= [2.443303419949772, 2.445003360252620, 2.435206281190796,
       2.397229070763543, 2.349096435450708, 2.296697511366329,
       2.226202284287925, 2.132775562815938, 2.027583832063605, 
       1.910630360626907, 1.778445679008138, 1.635079922339590,
       1.484330005656411, 1.324367385828326, 1.162101456042379,
       1.027709774120046, 0.9531004737753517, 0.9291925921464999,
       0.9306128757867301, 0.9373734199854492, 0.9374716427460326,
       0.9309174125983126, 0.9309062023801610, 0.9330912388264760,
       0.9290759363555934, 0.9313924528732267, 0.9325192643910063, 
       0.9298801313664510, 0.9355275529633224, 0.9383072037111380, 
       0.9315511559488159, 0.9286559025958544, 0.9465601634341064, 
       1.000382927938164, 1.092264926894742, 1.215568147385899, 
       1.343869169462016, 1.443290163237878, 1.545839461727406, 
       1.659387736523215, 1.785869422676203, 1.919004656822144, 
       2.048477467998743, 2.170484269722559, 2.273692536052303, 
       2.345258986779728, 2.405126351903085, 2.441499638224660, 
       2.444671091311594, 2.441550619712000, 2.443303419949772 ]
       
zd3d =[0.1070010083064705, 0.2407488482372501, 0.3804810616106032, 
       0.5159387361983736, 0.6503946165564584, 0.7935848357649634, 
       0.9352181481113857, 1.059155770383247, 1.178561334034986,
       1.293970423413222, 1.384524187526246, 1.432284507191533,
       1.453445898862813, 1.463657741995547, 1.445934022958783,
       1.359242878727507, 1.195323634649968, 1.004001380805289, 
       0.8237611508094629, 0.6662857938612209, 0.5291230710159969,
       0.4026518905384536, 0.2771781125246763, 0.1533390886307173, 
       2.974376768179259E-02, -9.479066317849494E-02, -0.2219425844279092,
       -0.3553396294163588,-0.4901821828518301, -0.6351465365923717,
       -0.8027612626835731, -0.9874117999240944,-1.170409819654592,
       -1.321107993701108, -1.416381737427392, -1.457845801192058,
       -1.463452619441855, -1.457140387484942,-1.446593354617252,
       -1.427322714872727,-1.380546621923976, -1.286917884397917,
       -1.155789433396968, -1.012819480748024,-0.8474364195572380,
       -0.6613050759653780,-0.4935139420070444, -0.3285717485346793,
       -0.1676609207868461,-2.517907021385858E-02,0.1070010083064715 ]

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
coilList[0].plot_2D_coil(ax,'b')
coilList[1].plot_2D_coil(ax,'g')
coilList[2].plot_2D_coil(ax,'b')
coilList[3].plot_2D_coil(ax,'g')
plt.plot(rd3d,zd3d,label="Wall",c='k',ls=':')
ax.set_aspect(aspect=1.0)
ax.set_xlabel('R [m]')
ax.set_ylabel('Z [m]')
plt.title("Error Field Coils")
plt.tight_layout()
plt.show()
if (scope):
  sys.exit()

for iNode in range(nodeXYZ.shape[1]):
    print("Calculating node: " + str(iNode))
    for iPhi in range(nodeXYZ.shape[2]):
        print("Calculating plane: " + str(iPhi))
        sys.stdout.flush()
        bXYZ=np.zeros(3)
        for iCoil in coilList:
            bXYZ[:]+=bs.intCoil(iCoil,nodeXYZ[:,iNode,iPhi])
        phi = 2.0*np.pi*iPhi/phiPlanes
### transform to bRZPhi
# bRZPhi accounts for the negative in Bphi due to rzphi coordinates
        bRZPhi[0,iNode,iPhi] = bXYZ[0]*np.cos(phi)+bXYZ[1]*np.sin(phi)
        bRZPhi[1,iNode,iPhi] = bXYZ[2]
        bRZPhi[2,iNode,iPhi] = bXYZ[0]*np.sin(phi)-bXYZ[1]*np.cos(phi)

bRPhase=np.fft.fft(bRZPhi[0,:,:],axis=1)/(float(phiPlanes))
bZPhase=np.fft.fft(bRZPhi[1,:,:],axis=1)/(float(phiPlanes))
bPhiPhase=np.fft.fft(bRZPhi[2,:,:],axis=1)/(float(phiPlanes))

### write brmp files
if (phiPlanes % 2 == 0): #even
    maxnphi = int(phiPlanes/2)
else: #odd
    maxnphi = int((phiPlanes+1)/2)
for ii in range (maxnphi +1):
    if ii==maxnphi:
        fac=0.5
    else:
        fac=1.0
    print(ii, maxnphi, fac)
    tempFileName = filePath + baseFileName +"{0:0=2d}".format(ii)  + fileExt
    thisFile = open(tempFileName,'w')
    for jj in range(bRPhase.shape[0]):
        thisLine ='{: 16.16e}'.format(fac*bRPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bRPhase[jj,ii].imag) + ", "
        thisLine+='{: 16.16e}'.format(fac*bZPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bZPhase[jj,ii].imag) + ", "
        thisLine+='{: 16.16e}'.format(fac*bPhiPhase[jj,ii].real) + ", " 
        thisLine+='{: 16.16e}'.format(fac*bPhiPhase[jj,ii].imag) + "\n"
        thisFile.write(thisLine)
    thisFile.close()

#for iCoil in coilList:
#    bvec+=bs.intCoil(iCoil,np.asarray([0.0,0.0,0.0]))

