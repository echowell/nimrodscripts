#!/usr/bin/python

import struct
import numpy as np
#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from scipy import interpolate

# counts based on the fact that the end of the file has two integer zeros (8) and each time adds numfields*4+8 (beginning and ending codes for data)
def count_time(filename,numfields):
  f=open(filename, 'rb')
  f.seek(0,2)
  fsize=f.tell()
  f.close
  return (fsize-8)/(numfields*4+8)


largef=16
smallf=8
titlef=18
tickf=15


# file size = 4* (49*(mx*pd+1)+2)*(my*pd+1)

#plt.rcParams['text.usetex']=True

int_size =4
float_size = 4

# list of file names to be read
files=['fgprofs.bin']
# Create dictionaries for each file with pertinent information
# Not necessary for files with time data as the "1" dimension
npx={}
npy={}

# data in binary with number of fields in numfields
numfields=16
psifac={}
rho={}
frb={}
mu0_p={}
qg={}
mach={}
n_den={}
mu0_pe={}
omgeb={}
Kpol={}
omgp={}
hot_den={}
hot_temp={}
ti={}
te={}
residual={}

#mu=(m**-1)

for findex in range(len(files)):
  fieldcount=count_time(files[findex],numfields)
  npx[files[findex]]=fieldcount

  psifac[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  rho[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  frb[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  mu0_p[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  qg[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  mach[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  n_den[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  mu0_pe[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  omgeb[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  Kpol[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  omgp[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  hot_den[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  hot_temp[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  ti[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  te[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')
  residual[files[findex]]=np.zeros(npx[files[findex]], dtype = 'f',order = 'F')

  f=open(files[findex], 'rb')

  jj=0
  f.seek(0)

  while jj < npx[files[findex]] :
    #first bite is length of a string
    temp=f.read(4)
    blah=struct.unpack(">l",temp)
    temp=f.read(numfields*4)
    tf = struct.unpack(">"+numfields*'f', temp)
    #last byte is length of string read
    temp=f.read(4)
    psifac[files[findex]][jj]=tf[0]
    rho[files[findex]][jj]=tf[1]
    frb[files[findex]][jj]=tf[2]
    mu0_p[files[findex]][jj]=tf[3]
    qg[files[findex]][jj]=tf[4]
    mach[files[findex]][jj]=tf[5]
    n_den[files[findex]][jj]=tf[6]
    mu0_pe[files[findex]][jj]=tf[7]
    omgeb[files[findex]][jj]=tf[8]    
    Kpol[files[findex]][jj]=tf[9]
    omgp[files[findex]][jj]=tf[10]
    hot_den[files[findex]][jj]=tf[11]
    hot_temp[files[findex]][jj]=tf[12]
    ti[files[findex]][jj]=tf[13]
    te[files[findex]][jj]=tf[14]
    residual[files[findex]][jj]=tf[15]
    
    jj=jj+1

  f.close()

#for i in range(len(qg[files[0]])):
#    qg[files[0]][i]=qg[files[0]][i]

for i in range(len(qg[files[0]])):
    mid2=(qg[files[0]][i]+2.)*(qg[files[0]][i+1]+2.)
    if mid2<0:
        irho2=i
    mid3=(qg[files[0]][i]+3.)*(qg[files[0]][i+1]+3.)
    if mid3<0:
        irho3=i
    mid4=(qg[files[0]][i]+4.)*(qg[files[0]][i+1]+4.)
    if mid4<0:
        irho4=i
        break

fig1,ax1=plt.subplots(1,2,figsize=(16,5))
#plt.subplots_adjust(left=0.10, bottom=0.12, right=0.95, top=0.92, wspace=0.175)

plt.locator_params(axis='y',nbins=4)

mup=[]
for i in range(len(mu0_p[files[0]])):
    mup.append(100*mu0_p[files[0]][i])

#ax1.set_title(r'$f$ v. $\rho$',fontsize=titlef)
plt.setp(ax1[0].get_xticklabels(), fontsize=20)
plt.setp(ax1[0].get_yticklabels(), fontsize=20)
ax1[0].plot(rho[files[0]][:],mup,'b',label=r'$100\times\mu_0p$',lw=3,color='b')
ax1[0].plot(rho[files[0]][:],abs(qg[files[0]][:]),'b',label=r'$abs(q)$',lw=3,color='r')
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax1[0].axvline(x=rho[files[0]][irho2],lw=3,color='g',ls='dashed',label=r'$q=2$')
ax1[0].axvline(x=rho[files[0]][irho3],lw=3,color='m',ls='dashed',label=r'$q=3$')
ax1[0].axvline(x=rho[files[0]][irho4],lw=3,color='orange',ls='dashed',label=r'$q=4$')

#ax2=ax1.twinx()
#ax2.plot(rho[files[0]][:],omgeb[files[0]][:],'b',label=r'$\Omega_{E\times B}$',lw=3,color='g')
#ax1.plot(0,0,'b',label=r'$\Omega_{E\times B}$ v. $\rho$',lw=3,color='g')

ax1[0].axhline(y=1,lw=2,ls='dotted',color='k')


ax1[0].legend(fontsize=15,loc=2,ncol=2)

#ax2.yaxis.major.formatter._useMathText = True
#ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
#ax1[0].yaxis.offsetText.set_fontsize(16)



#ax.set_xlim(0.6,0.8)
ax1[0].set_ylim(0,5.5)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')

#plt.subplots_adjust(left=0.10, bottom=0.12, right=0.95, top=0.92, wspace=0.175)


#ax1.set_title(r'$f$ v. $\rho$',fontsize=titlef)
plt.setp(ax1[1].get_xticklabels(), fontsize=20)
plt.setp(ax1[1].get_yticklabels(), fontsize=20)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax1[1].axvline(x=rho[files[0]][irho2],lw=3,color='g',ls='dashed',label=r'$q=2$')
ax1[1].axvline(x=rho[files[0]][irho3],lw=3,color='m',ls='dashed',label=r'$q=3$')
ax1[1].axvline(x=rho[files[0]][irho4],lw=3,color='orange',ls='dashed',label=r'$q=4$')

ax1[1].plot(rho[files[0]][:],omgeb[files[0]][:],'b',label=r'$\Omega_{E\times B}$',lw=5,color='b')


ax1[1].axhline(y=0,lw=2,ls='dotted',color='k')

ax1[1].yaxis.major.formatter._useMathText = True
ax1[1].ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
ax1[1].yaxis.offsetText.set_fontsize(18)
ax1[1].tick_params(axis='y',labelsize=20)


#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax1[0].set_xlabel(r'$\rho$',fontsize=26)
ax1[1].set_xlabel(r'$\rho$',fontsize=26)
ax1[1].set_ylabel(r'$\Omega_{E\times B}$',fontsize=26)

#plt.savefig('fgprofs2.png',bbox_inches='tight')


plt.show()


fig1=plt.figure(figsize=(12,12))
#plt.subplots_adjust(left=0.10, bottom=0.12, right=0.95, top=0.92, wspace=0.175)
ax=fig1.add_subplot(331)
plt.title(r'$f$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],frb[files[0]][:],'b',label=r'$f$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$f$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(332)
plt.title(r'$\mu_0 p$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],mu0_p[files[0]][:],'b',label=r'$\mu_0p$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$mu_0p$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(333)
plt.title(r'$q$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],qg[files[0]][:],'b',label=r'$q$')
ax.axvline(x=rho[files[0]][irho2],color='g')
ax.axvline(x=rho[files[0]][irho3],color='m')
ax.axvline(x=rho[files[0]][irho4],color='orange')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$q$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(334)
plt.title(r'$Mach$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],mach[files[0]][:],'b',label=r'${\rm Mach}$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'${\rm Mach}$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(335)
plt.title(r'$\Omega_{E\times B}$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],omgeb[files[0]][:],'b',label=r'$\Omega_{E\times B}$ v. $\rho$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$\Omega_{E\times B}$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(336)
plt.title(r'$\Omega_{\nabla P}$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],omgp[files[0]][:],'b',label=r'${\rm Mach}$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$\Omega_{\nabla P}$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(337)
plt.title(r'$n$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],n_den[files[0]][:],'b',label=r'$q$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$n$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

ax=fig1.add_subplot(338)
plt.title(r'$ln(Residual)$ v. $\rho$',fontsize=titlef)
plt.setp(ax.get_xticklabels(), fontsize=tickf)
plt.setp(ax.get_yticklabels(), fontsize=tickf)
#sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
ax.plot(rho[files[0]][:],residual[files[0]][:],'b',label=r'$q$')
#ax.set_xlim(0.6,0.8)
#ax.set_ylim(-1,0.25)
#ax.legend(loc=2,prop={'size':8})
#plt.axvline(x=rs, color='k')
ax.set_ylabel(r'$ln(Residual)$',fontsize=largef)
ax.set_xlabel(r'$\rho$',fontsize=largef)

#plt.savefig('fgprofs1.png',bbox_inches='tight')

plt.show()
