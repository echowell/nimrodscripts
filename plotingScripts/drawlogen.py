#!/usr/local/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import math as m

data = np.loadtxt('/home/research/ehowell/SCRATCH/166439/03300_2_fl/19091702/lphi5_nolinear_fresh/logen.txt')

print (data.shape)
vals = data.shape
#print vals[0]
fmode =11
lenl=int(vals[0]/fmode)
time = np.zeros(int(vals[0]/fmode),dtype='float')
logme = np.zeros((fmode, int(vals[0]/fmode)),dtype='float')
logke = np.zeros((fmode, int(vals[0]/fmode)),dtype='float')
me = np.zeros((fmode, int(vals[0]/fmode)),dtype='float')
ke = np.zeros((fmode, int(vals[0]/fmode)),dtype='float')


for ii in range(0,lenl):
  time[ii] = data[ii*fmode,1]*1000000
  for jj in range(0,fmode):
    logme[jj,ii]=data[ii*fmode+jj,4]
    logke[jj,ii]=data[ii*fmode+jj,5]
    me[jj,ii]=m.pow(10,data[ii*fmode+jj,4])
    ke[jj,ii]=m.pow(10,data[ii*fmode+jj,5])
    



fig=plt.figure(figsize=(10,4))
#fig.subplots_adjust(left=0.5)
#ax.yaxis.labelpad=35
ax = fig.add_subplot(121)
#ax.plot(time,logme[0,:],ls='-.',label='n=0')
ax.plot(time,logme[1,:],ls='-',label='n=1')
ax.plot(time,logme[2,:],ls=':',label='n=2')
ax.plot(time,logme[3,:],ls=':',label='n=3')
ax.plot(time,logme[4,:],ls='-',label='n=4')
ax.plot(time,logme[5,:],ls='-.',label='n=5')
ax.plot(time,logme[6,:],ls='-',label='n=6')
ax.plot(time,logme[7,:],ls=':',label='n=7')
ax.plot(time,logme[8,:],ls=':',label='n=8')
ax.plot(time,logme[9,:],ls='-',label='n=9')
ax.plot(time,logme[10,:],ls='-.')
ax.set_ylabel(r"Log Energy", size=16)
ax.set_xlabel(r"Time($\mu$s)",size=16)
#plt.legend(ncol=2, loc=4)
#ax.set_xlim([0.0,0.5])
ax.set_ylim([-1,2])
#ax.yaxis.labelpad=35
plt.title(r'Magnetic Energy', size=16)
ax = fig.add_subplot(122)
#ax.plot(time,logke[0,:],ls='-.',label='n=0')
ax.plot(time,logke[1,:],ls='-',label='n=1')
ax.plot(time,logke[2,:],ls=':',label='n=2')
ax.plot(time,logke[3,:],ls=':',label='n=3')
ax.plot(time,logke[4,:],ls='-',label='n=4')
ax.plot(time,logke[5,:],ls='-.',label='n=5')
ax.plot(time,logke[6,:],ls='-',label='n=6')
ax.plot(time,logke[7,:],ls=':', label='n=7')
ax.plot(time,logke[8,:],ls=':', label='n=8')
ax.plot(time,logke[9,:],ls='-',label='n=9')
ax.plot(time,logke[10,:],ls='-.')
ax.set_ylabel(r"Log Energy", size=16)
ax.set_xlabel(r"Time($\mu$s)",size=16)
#ax.set_xlim([0.0,0.5])
ax.set_ylim([-3.0,1])
plt.legend(ncol=3, loc=4, columnspacing =1,prop={'size':12})
#ax.yaxis.labelpad=35
plt.title(r'Kinetic Energy', size=16)
plt.tight_layout()
#ax.legend(loc=2)
plt.show()

fig=plt.figure(figsize=(10,4))
#fig.subplots_adjust(left=0.5)
#ax.yaxis.labelpad=35
ax = fig.add_subplot(121)
ax.plot(time,me[0,:],ls='-.',label='n=0')
ax.set_ylabel(r"Energy", size=16)
ax.set_xlabel(r"Time($\mu$s)",size=16)
#plt.legend(ncol=2, loc=4)
#ax.set_xlim([0.0,0.5])
#ax.set_ylim([-1,2])
#ax.yaxis.labelpad=35
plt.title(r'Magnetic Energy', size=16)
ax = fig.add_subplot(122)
ax.plot(time[1:],ke[0,1:],ls='-.',label='n=0')
ax.set_ylabel(r"Energy", size=16)
ax.set_xlabel(r"Time($\mu$s)",size=16)
#ax.set_xlim([0.0,0.5])
#ax.set_ylim([4.14,4.16])
plt.legend(ncol=3, loc=1, columnspacing =1,prop={'size':12})
#ax.yaxis.labelpad=35
plt.title(r'Kinetic Energy', size=16)
plt.tight_layout()
#ax.legend(loc=2)
plt.show()
