#!/usr/bin/python

import struct
import numpy as np
#from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
from fgProfs import rho,qg,irho2,irho3,irho4


# xy have elecd=0.25
# file size = 4* (49*(mx*pd+1)+2)*(my*pd+1)

#plt.rcParams['text.usetex']=True

int_size = 4
float_size = 4
nimmx=144
nimmy=128
nimpd=5

plotfields=0
# list of file names to be read
files=['/home/research/ehowell/SCRATCH/166439/footpoint_03300_q104/lphi5/S7Pr1e2_surfmn/xy_slice01.bin']
files=['/home/research/ehowell/SCRATCH/166439/footpoint_03300_q104/lphi5/vac_surfmn/xy_slice01.bin']
files=['/home/research/ehowell/SCRATCH/166439/footpoint_03300_q104/lphi5/vac_surfmn/xy_slice01.bin']
#,'xy_slice17.bin','xy_slice260.bin']
#files=['xy_slice15000.bin']
# slicesy is between 0 and my*pd for all yslices
slicesy=[0]
color=cm.rainbow(np.linspace(0,1,len(slicesy)))
#I assume all files have the same  size
mx={files[0]:nimmx}#,files[1]:128,files[2]:128}
my={files[0]:nimmy}#,files[1]:72,files[2]:72}
pd={files[0]:nimpd}#,files[1]:5,files[2]:5}
npx={}
npy={}

tf={}
ix={}
iy={}
R={}
Z={}
B0R={}
B0Z={}
B0T={}
J0R={}
J0Z={}
J0T={}
V0R={}
V0Z={}
V0T={}
P0={}
PE0={}
n0={}
diff={}
BRr={}
BZr={}
BTr={}
BRi={}
BZi={}
BTi={}
JRr={}
JZr={}
JTr={}
JRi={}
JZi={}
JTi={}
VRr={}
VZr={}
VTr={}
VRi={}
VZi={}
VTi={}
Pr={}
Pi={}
PEr={}
PEi={}
Nr={}
Ni={}
Cr={}
Ci={}
TEr={}
TEi={}
TIr={}
TIi={}
N={}

#each field will be a 2x2 array of size(mx*pd+1,my*pd+1)
for findex in range(len(files)):
  npx[files[findex]]=mx[files[findex]]*pd[files[findex]]+1
  npy[files[findex]]=my[files[findex]]*pd[files[findex]]+1

  tf[files[findex]]=np.zeros(47,dtype='f',order = 'F')
  ix[files[findex]]=np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  iy[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  R[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Z[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  B0R[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  B0Z[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  B0T[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  J0R[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  J0Z[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  J0T[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  V0R[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  V0Z[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  V0T[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  P0[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  PE0[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  n0[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  diff[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BRr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BZr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BTr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BRi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BZi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  BTi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JRr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JZr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JTr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JRi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JZi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  JTi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VRr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VZr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VTr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VRi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VZi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  VTi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Pr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Pi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  PEr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  PEi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Nr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Ni[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Cr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  Ci[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  TEr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  TEi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  TIr[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  TIi[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')
  N[files[findex]]= np.zeros((npx[files[findex]],npy[files[findex]]), dtype = 'f',order = 'F')




  f=open(files[findex], 'rb')
  limity=my[files[findex]]*pd[files[findex]]+1
  limitx=mx[files[findex]]*pd[files[findex]]+1

  jj=0
  f.seek(0)

  while jj < limity :
    ii=0
    while ii < limitx :
      #first bite is length of a string
      temp=f.read(4)
      blah=struct.unpack(">l",temp)
#      if jj==0 and (ii==0 or ii==1):
#        print blah
      temp=f.read(188)
      tf[files[findex]] = struct.unpack(">"+47*'f', temp)
      if jj==0 and (ii==0 or ii==1):
        tf1=tf
#        print tf1,'\n'
      #last byte is length of string read
      temp=f.read(4)
      blah=struct.unpack(">l",temp)
#      if jj==0 and (ii==0 or ii==1):
#        print blah

      ix[files[findex]][ii,jj]=tf[files[findex]][0]
      iy[files[findex]][ii,jj]=tf[files[findex]][1]
      R[files[findex]][ii,jj]=tf[files[findex]][2]
      Z[files[findex]][ii,jj]=tf[files[findex]][3]
      B0R[files[findex]][ii,jj]=tf[files[findex]][4]
      B0Z[files[findex]][ii,jj]=tf[files[findex]][5]
      B0T[files[findex]][ii,jj]=tf[files[findex]][6]
      J0R[files[findex]][ii,jj]=tf[files[findex]][7]
      J0Z[files[findex]][ii,jj]=tf[files[findex]][8]
      J0T[files[findex]][ii,jj]=tf[files[findex]][9]
      V0R[files[findex]][ii,jj]=tf[files[findex]][10]
      V0Z[files[findex]][ii,jj]=tf[files[findex]][11]
      V0T[files[findex]][ii,jj]=tf[files[findex]][12]
      P0[files[findex]][ii,jj]=tf[files[findex]][13]
      PE0[files[findex]][ii,jj]=tf[files[findex]][14]
      n0[files[findex]][ii,jj]=tf[files[findex]][15]
      diff[files[findex]][ii,jj]=tf[files[findex]][16]
      BRr[files[findex]][ii,jj]=tf[files[findex]][17]
      BZr[files[findex]][ii,jj]=tf[files[findex]][18]
      BTr[files[findex]][ii,jj]=tf[files[findex]][19]
      BRi[files[findex]][ii,jj]=tf[files[findex]][20]
      BZi[files[findex]][ii,jj]=tf[files[findex]][21]
      BTi[files[findex]][ii,jj]=tf[files[findex]][22]
      JRr[files[findex]][ii,jj]=tf[files[findex]][23]
      JZr[files[findex]][ii,jj]=tf[files[findex]][24]
      JTr[files[findex]][ii,jj]=tf[files[findex]][25]
      JRi[files[findex]][ii,jj]=tf[files[findex]][26]
      JZi[files[findex]][ii,jj]=tf[files[findex]][27]
      JTi[files[findex]][ii,jj]=tf[files[findex]][28]
      VRr[files[findex]][ii,jj]=tf[files[findex]][29]
      VZr[files[findex]][ii,jj]=tf[files[findex]][30]
      VTr[files[findex]][ii,jj]=tf[files[findex]][31]
      VRi[files[findex]][ii,jj]=tf[files[findex]][32]
      VZi[files[findex]][ii,jj]=tf[files[findex]][33]
      VTi[files[findex]][ii,jj]=tf[files[findex]][34]
      Pr[files[findex]][ii,jj]=tf[files[findex]][35]
      Pi[files[findex]][ii,jj]=tf[files[findex]][36]
      PEr[files[findex]][ii,jj]=tf[files[findex]][37]
      PEi[files[findex]][ii,jj]=tf[files[findex]][38]
      Nr[files[findex]][ii,jj]=tf[files[findex]][39]
      Ni[files[findex]][ii,jj]=tf[files[findex]][40]
      Cr[files[findex]][ii,jj]=tf[files[findex]][41]
      Ci[files[findex]][ii,jj]=tf[files[findex]][42]
      TEr[files[findex]][ii,jj]=tf[files[findex]][43]
      TEi[files[findex]][ii,jj]=tf[files[findex]][44]
      TIr[files[findex]][ii,jj]=tf[files[findex]][45]
      TIi[files[findex]][ii,jj]=tf[files[findex]][46]
      ii=ii+1
    if (jj<(limity-1)):
      temp= f.read(8)
    jj=jj+1

  f.close()



for m in range(len(files)):
  for i in range(len(BRr[files[m]][:,0])-1):
    for j in slicesy:
      N[files[m]][i,j]=-(R[files[m]][i,j]/BZr[files[m]][i,j])*(BZr[files[m]][i+1,j]-BZr[files[m]][i,j])/(R[files[m]][i+1,j]-R[files[m]][i,j])
      if BZr[files[m]][i,j]==0:
#        print "BZr=0",BZr[files[m]][i,j]
        N[files[m]][i,j]=0
      if (R[files[m]][i+1,j]-R[files[m]][i,j])==0:
#        print "R=0",R[files[m]][i+1,j],R[files[m]][i,j]
        N[files[m]][i,j]=0
  N[files[m]][-1,j]=-(R[files[m]][-1,j]/BZr[files[m]][-1,j])*(BZr[files[m]][-1,j]-BZr[files[m]][-2,j])/(R[files[m]][-1,j]-R[files[m]][-2,j])


qcon1=np.zeros((npx[files[findex]],npy[files[findex]]))
#fgfile='fgprofs.bin'
fgfile='/home/research/ehowell/SCRATCH/166439/footpoint_03300_q104/lphi5/S7Pr1e2_surfmn/fgprofs.bin'


for j in range(npy[files[0]]):
    for i in range(qg[fgfile].size):    
        qcon1[i+1,j]=qg[fgfile][i]
    qcon1[i+2,j]=-5.
        

clines=301

pr_max = np.max(J0T[files[0]])
pr_min = np.min(J0T[files[0]])

if abs(pr_max)>abs(pr_min):
    pr_min=-pr_max
else:
    pr_max=-pr_min

lev,delta=np.linspace(pr_min,pr_max,clines,retstep=True)

qlev=np.array([-4,-3,-2])        

if plotfields==1:
    fig,ax=plt.subplots(1,2,figsize=(10,6),sharey=True)                
    CS=ax[0].contourf(R[files[0]],Z[files[0]],J0T[files[0]],levels=lev,cmap=cm.seismic)

    QS=ax[0].contour(R[files[0]],Z[files[0]],qcon1,levels=qlev,colors='k',linestyles='solid')

    fmt = {}
    strs=[r'$q=\!-4$',r'$q=\!-3$',r'$q=\!-2$']
    for l, s in zip(QS.levels, strs):
        fmt[l] = s

    manual_locations=[(1.36,0.85),(1.42,0.73),(1.48,0.59)]

    ax[0].clabel(QS,qlev,inline=1,fmt=fmt,manual=manual_locations,inline_spacing=15)

    ax[0].set_aspect('equal')

    plt.setp(ax[0].get_xticklabels(), fontsize=18)
    plt.setp(ax[0].get_yticklabels(), fontsize=18)
    ax[0].set_title(r'Real $(J_{\phi})$',fontsize=22)
    ax[0].set_xlabel(r'${\rm R\,(m)}$',fontsize=20)
    ax[0].set_ylabel(r'${\rm Z\,(m)}$',fontsize=20)
    cbar=plt.colorbar(CS,pad=0.01)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(16)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()

    CS=ax[1].contourf(R[files[0]],Z[files[0]],J0T[files[0]],levels=lev,cmap=cm.seismic)

    QS=ax[1].contour(R[files[0]],Z[files[0]],qcon1,levels=qlev,colors='k',linestyles='solid')

    ax[1].clabel(QS,qlev,inline=1,fmt=fmt,manual=manual_locations,inline_spacing=15)

    ax[1].set_aspect('equal')

    plt.setp(ax[1].get_xticklabels(), fontsize=18)
    plt.setp(ax[1].get_yticklabels(), fontsize=18)
    ax[1].set_title(r'Real $(J_{\phi})$',fontsize=22)
    ax[1].set_xlabel(r'${\rm R\,(m)}$',fontsize=20)

    plt.savefig('xy_J0T.png',bbox_inches='tight')                                

    fig,ax=plt.subplots(2,1,figsize=(5,10),sharex=True) 

    plt.setp(ax[0].get_xticklabels(), fontsize=16)
    plt.setp(ax[0].get_yticklabels(), fontsize=16)

    jj=0
    for j in slicesy:
      ax[0].plot(R[files[0]][:,j],-J0T[files[0]][:,j],c='red',lw=3,label=r'With SOL')
      ax[0].plot(R[files[0]][:,j],-J0T[files[0]][:,j],c='blue',lw=3,label=r'Without SOL')
        
      jj=jj+1

    ax[0].set_xlim(2.23,2.25)
    ax[0].set_ylim(-10000,80000)

    ax[0].yaxis.major.formatter._useMathText = True
    ax[0].ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax[0].yaxis.offsetText.set_fontsize(16)

    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax[0].set_ylabel(r'$J_{0,\phi}\,({\rm A/m^2})$',fontsize=18)
    #ax[0].set_xlabel(r'$R\,({\rm m})$',fontsize=18)

    plt.legend(loc=1,fontsize=14)

    plt.setp(ax[1].get_xticklabels(), fontsize=16)
    plt.setp(ax[1].get_yticklabels(), fontsize=16)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax[1].plot(R[files[0]][:,j],P0[files[0]][:,j],c='red',lw=3,label=r'With SOL')
      ax[1].plot(R[files[0]][:,j],P0[files[0]][:,j],c='blue',lw=3,label=r'Without SOL')
        
      jj=jj+1

    ax[1].set_xlim(2.23,2.25)
    ax[1].set_ylim(0,750)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax[1].set_ylabel(r'$P_0\,({\rm Pa})$',fontsize=18)
    ax[1].set_xlabel(r'$R\,({\rm m})$',fontsize=18)

    plt.legend(loc=1,fontsize=14)

    plt.savefig('xy_P0compare.png',bbox_inches='tight')  

    fig=plt.figure(figsize=(12,12))

    ax=fig.add_subplot(431)
    plt.title(r'$B_R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BRr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(432)
    plt.title(r'$B_Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BZr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(433)
    plt.title(r'$B_\phi$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BTr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_\phi$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(434)
    plt.title(r'$B_R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])

    for j in slicesy:
      ax.plot(R[files[0]][:,j],BRi[files[0]][:,j],label=r'$B_R$')

    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(435)
    plt.title(r'$B_Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BZi[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(436)
    plt.title(r'$B_\phi$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BTi[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_\phi$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    plt.tight_layout()
    filename=files[0].split('.')
    plt.savefig(filename[0]+"_B1_findex.png",bbox_inches='tight')                                                                                                
                                                                                                                                                                                                                                                                                                    
    fig=plt.figure(figsize=(12,12))

    ax=fig.add_subplot(431)
    plt.title(r'$B0R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],B0R[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B0R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(432)
    plt.title(r'$B0Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],B0Z[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B0Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(433)
    plt.title(r'$B0T$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],B0T[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B0T$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(434)
    plt.title(r'$V0T$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],V0T[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$V0T$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(435)
    plt.title(r'$P0$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],P0[files[0]][:,j],c=color[jj],label=r'$B_R$')
      
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$P0$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(436)
    plt.title(r'$n0$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],n0[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$n0$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(437)
    plt.title(r'$J0R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],J0R[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J0R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(438)
    plt.title(r'$J0Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],J0Z[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J0Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(439)
    plt.title(r'$J0T$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],J0T[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J0T$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    plt.tight_layout()
    filename=files[0].split('.')
    plt.savefig(filename[0]+"_eq_findex.png",bbox_inches='tight')


    fig=plt.figure(figsize=(12,12))
    #plt.subplots_adjust(left=0.10, bottom=0.12, right=0.95, top=0.92, wspace=0.175)
    ax=fig.add_subplot(431)
    plt.title(r'$B_R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BRr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1

    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(432)
    plt.title(r'$B_Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BZr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(433)
    plt.title(r'$B_\phi$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],BTr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$B_\phi$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)



    ax=fig.add_subplot(434)
    plt.title(r'$V_R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],VRr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$V_R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(435)
    plt.title(r'$V_Z$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],VZr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$V_Z$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(436)
    plt.title(r'$V_\phi$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],VTr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')

    sm = cm.ScalarMappable(cmap=cm.rainbow, norm=plt.Normalize(vmin=min(slicesy), vmax=max(slicesy)))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cx=plt.colorbar(sm)
    cx.set_ticks(slicesy)

    ax.set_ylabel(r'$V_\phi$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(437)
    plt.title(r'$J_R$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],JRr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J_R$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(438)
    plt.title(r'$J_\theta$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],JZr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J_\theta$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(439)
    plt.title(r'$J_\phi$ v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],JTr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'$J_\phi$',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(4,3,10)
    plt.title(r'Real($P$) v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],Pr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'Real($P$)',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(4,3,11)
    plt.title(r'Real($T_i$) v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],TIr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'Real($T_i$)',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)

    ax=fig.add_subplot(4,3,12)
    plt.title(r'Real($n$) v. $R$',fontsize=12)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    #sca=np.amax([abs(br1[:,13]),abs(bpi1[:,13]),abs(bbi1[:,13])])
    jj=0
    for j in slicesy:
      ax.plot(R[files[0]][:,j],Nr[files[0]][:,j],c=color[jj],label=r'$B_R$')
      jj=jj+1
    #ax.set_xlim(0.6,0.8)
    #ax.set_ylim(-1,0.25)
    #ax.legend(loc=2,prop={'size':8})
    #plt.axvline(x=rs, color='k')
    ax.set_ylabel(r'Real($n$)',fontsize=12)
    ax.set_xlabel(r'$r$',fontsize=12)
    plt.tight_layout()
    #plt.show()
    filename=files[0].split('.')
    if len(filename)==3:
      plt.savefig(filename[0]+'_findex.'+filename[1]+".png",bbox_inches='tight')
    else:
      plt.savefig(filename[0]+"_findex.png",bbox_inches='tight')

    plt.show()


