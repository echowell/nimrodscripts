#!/usr/local/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os

homeDir = os.environ['HOME']
relDir = "/SCRATCH/174446_novac_debug/nonlin1_eq26_bamp3_nimfl/4000/"
fileName = "nimfl0004000.dat"

fullFile = homeDir+relDir+fileName
data = np.loadtxt(fullFile)
sizes=0.005
sizes2=0.05

surfaces = []
n3m1 = [53,63,66,97,98,100,116,139,236,240,264,276,283,287]
n2m1 = [94,109,170,199,252]
print(data.shape)
last_surface = -1
n_surface = 0
for ii in range(data.shape[0]):
  if data[ii,3]== last_surface:
    n_surface = n_surface+1
  else:
    if n_surface > 0:
      this_surface = np.zeros([n_surface,2])
      this_surface = data[i_start:ii,0:2]
      surfaces.append(this_surface)
    i_start = ii
    n_surface = 1
    last_surface = data[ii,3]
this_surface = data[i_start:,0:2]
surfaces.append(this_surface)

print(len(surfaces))
lastl =0.0
nskip=20


fig=plt.figure(figsize=(4,6))
ax = fig.add_subplot(111)
for i_surface in range(0,len(surfaces),nskip):
  this_surface = surfaces[i_surface]
  ax.scatter(this_surface[:,0],this_surface[:,1],s=sizes,c='k')
for i_surface in n3m1:
  this_surface = surfaces[i_surface]
  ax.scatter(this_surface[:,0],this_surface[:,1],s=sizes2)
for i_surface in n2m1:
  this_surface = surfaces[i_surface]
  ax.scatter(this_surface[:,0],this_surface[:,1],s=sizes2)

  
ax.set_xlabel(r"$R$", size=16)
ax.set_ylabel(r"$Z$",size=16,rotation=0)
plt.axis('equal')
ax.set_xlim([1,2.5])
plt.tight_layout()
plt.show()
