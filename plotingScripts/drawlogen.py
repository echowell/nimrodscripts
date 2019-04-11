import numpy as np
import matplotlib.pyplot as plt
import math as m

data = np.loadtxt('logen.txt')

print data.shape
vals = data.shape
print vals[0]
fmode =43
lenl=vals[0]/fmode
time = np.zeros(vals[0]/fmode,dtype='float')
logme = np.zeros((fmode, vals[0]/fmode),dtype='float')
logke = np.zeros((fmode, vals[0]/fmode),dtype='float')


for ii in range(0,lenl):
  time[ii] = data[ii*fmode,1]*1000
  for jj in range(0,fmode):
    logme[jj,ii]=data[ii*fmode+jj,4]
    logke[jj,ii]=data[ii*fmode+jj,5]
    



fig=plt.figure(figsize=(10,5))
#fig.subplots_adjust(left=0.5)
#ax.yaxis.labelpad=35
ax = fig.add_subplot(121)
ax.plot(time,logme[0,:],ls='-.',label='n=0')
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
ax.plot(time,logme[11,:],ls='-')
ax.plot(time,logme[15,:],ls='-.')
ax.plot(time,logme[20,:],ls='-.')
ax.set_ylabel(r"Log Energy", size=16)
ax.set_xlabel(r"Time(ms)",size=16)
ax.annotate('Island Coalescence', xy=(2.4,-0.6),xytext=(0.1,-4), arrowprops=dict(width=2,headwidth=3,facecolor='black', shrink=0.10))
ax.annotate('Sawtooth Oscillations', xy=(3.2,-4.8),xytext=(2.5,-10), arrowprops=dict(width=2,headwidth=3,facecolor='black', shrink=0.10))
ax.annotate('', xy=(4.0,-4.8),xytext=(4.0,-9.0), arrowprops=dict(width=2,headwidth=3,facecolor='black', shrink=0.10))
ax.annotate('', xy=(4.75,-4.8),xytext=(4.5,-9.0), arrowprops=dict(width=2,headwidth=3,facecolor='black', shrink=0.10))
#plt.legend(ncol=2, loc=4)
#ax.set_xlim([0.0,0.5])
ax.set_ylim([-17,6])
#ax.yaxis.labelpad=35
plt.title(r'Magnetic Energy', size=16)
ax = fig.add_subplot(122)
ax.plot(time,logke[0,:],ls='-.',label='n=0')
ax.plot(time,logke[1,:],ls='-',label='n=1')
ax.plot(time,logke[2,:],ls=':',label='n=2')
ax.plot(time,logke[3,:],ls=':',label='n=3')
ax.plot(time,logke[4,:],ls='-',label='n=4')
ax.plot(time,logke[5,:],ls='-.',label='n=5')
ax.plot(time,logke[6,:],ls='-',label='n=6')
ax.plot(time,logke[7,:],ls=':',label='n=7')
ax.plot(time,logke[8,:],ls=':',label='n=8')
ax.plot(time,logke[9,:],ls='-',label='n=9')
ax.plot(time,logke[10,:],ls='-.',label='n=10')
ax.plot(time,logke[11,:],ls='-',label='n=11')
ax.plot(time,logke[15,:],ls='-.',label='n=15')
ax.plot(time,logke[20,:],ls='-.',label='n=20')
ax.set_ylabel(r"Log Energy", size=16)
ax.set_xlabel(r"Time(ms)",size=16)
#ax.set_xlim([0.0,0.5])
ax.set_ylim([-20.0,1])
plt.legend(ncol=2, loc=4, columnspacing =1,prop={'size':12})
#ax.yaxis.labelpad=35
plt.title(r'Kinetic Energy', size=16)
plt.tight_layout()
#ax.legend(loc=2)
plt.show()
