import time
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from fieldIndex import files,R,Z,BRr,BRi,B0R,B0Z,B0T,qcon1,J0T,pd
from fgProfs import rho,qg,irho2,irho3,irho4


fgfile='/home/research/ehowell/SCRATCH/166439/footpoint_03300_q104/lphi5/S7Pr1e2_surfmn/fgprofs.bin'

mrange=10

Z0=Z[files[0]][0,0]
R0=R[files[0]][0,0]

r_minor=np.sqrt((R[files[0]]-R0)**2+(Z[files[0]]-Z0)**2)

#calculate length along poloidal extent of axisymmetric poloidal flux contour
s=np.zeros((len(R[files[0]]),len(R[files[0]][0])))
ds=np.zeros((len(R[files[0]]),len(R[files[0]][0])))
for i in range(len(s[:,0])):
    for j in range(len(s[0,:])-1):
        ds[i,j]=math.sqrt((R[files[0]][i][j+1]-R[files[0]][i][j])**2+(Z[files[0]][i][j+1]-Z[files[0]][i][j])**2)
        s[i,j+1]=s[i,j]+ds[i,j]

#calculate equilibrium poloidal field for [r,pol] locations
B0P=np.sqrt((B0R[files[0]])**2+(B0Z[files[0]])**2)
if J0T[files[0]][0,0]>0: B0P=-B0P



dqds=(-B0T[files[0]])/(2*math.pi*B0P*R[files[0]]**2)
q1=np.trapz(dqds,s,axis=1)
jac=q1[:,None]*R[files[0]]**3*B0P/(-B0T[files[0]]) 

#calculate straight-field line theta (PEST coordinate, Jim's derivation)
theta_str=np.zeros((len(R[files[0]]),len(R[files[0]][0])))
dtheta_str=np.zeros((len(R[files[0]]),len(R[files[0]][0])))
for i in range(len(theta_str[:,0])):
    for j in range(len(theta_str[0,:])-1):
        theta_str[i,j+1]=theta_str[i,j]+1./(q1[i]+1.0e-11)*(ds[i,j]*(-B0T[files[0]][i,j])/(B0P[i,j]*R[files[0]][i,j]**2))
        dtheta_str[i,j]=1./(q1[i]+1.0e-11)*(ds[i,j]*(-B0T[files[0]][i,j])/(B0P[i,j]*R[files[0]][i,j]**2))
    for j in range(len(theta_str[0,:])):
        theta_str[i,j]=theta_str[i,j]-theta_str[i,pd[files[0]]]


dFSAAdth=2*math.pi*jac
FSArea=np.trapz(dFSAAdth,theta_str,axis=1)

drhodth=2*math.pi*jac*r_minor/(FSArea[:,None]+1.0e-11)
rho1=np.trapz(drhodth,theta_str,axis=1)

rholcfs=rho1[int(len(rho1)*.75)]    
                        
for i in range(len(rho1)):
    rho1[i]=rho1[i]/rholcfs

for i in range(len(q1)):
    mid2=(q1[i]+2.)*(q1[i+1]+2.)
    if mid2<0:
        irho12=i
    mid3=(q1[i]+3.)*(q1[i+1]+3.)
    if mid3<0:
        irho13=i
    mid4=(q1[i]+4.)*(q1[i+1]+4.)
    if mid4<0:
        irho14=i
    mid5=(q1[i]+5.)*(q1[i+1]+5.)
    if mid5<0:
        irho15=i
    mid6=(q1[i]+6.)*(q1[i+1]+6.)
    if mid6<0:
        irho16=i
    mid10=(q1[i]+10.)*(q1[i+1]+10.)
    if mid10<0:
        irho110=i
        break


mmax=mrange
mmin=-mrange
m=np.linspace(mmin,mmax,mmax-mmin+1)

bcnm=np.zeros((len(files),(mmax-mmin+1),len(R[files[0]])))
bsnm=np.zeros((len(files),(mmax-mmin+1),len(R[files[0]])))

for l in range(len(bcnm)):
    if files[l]=='xy_slice00.bin':
        multfact=1/5.e-5
    else:
        multfact=1
    for k in range(len(bcnm[0])):
        dbcnmdth=2*np.pi/(FSArea[:,None]+1e-11)*jac*multfact*(BRr[files[l]]*np.cos((mmin+k)*theta_str)-BRi[files[l]]*np.sin((mmin+k)*theta_str))
        dbsnmdth=2*np.pi/(FSArea[:,None]+1e-11)*jac*multfact*(-BRr[files[l]]*np.sin((mmin+k)*theta_str)-BRi[files[l]]*np.cos((mmin+k)*theta_str))
        
        bcnm[l,k]=np.trapz(dbcnmdth,theta_str,axis=1)
        bsnm[l,k]=np.trapz(dbsnmdth,theta_str,axis=1)

bnm=np.sqrt(bcnm**2+bsnm**2)
#

#plot choices

plotqcon=0
plotqlin=0
plotjac=0
plottheta=0
plotqcompare=0
plotFSArea=0
plotsurf=0
plotsurf2wall=0
plotmag=0
plotmag2wall=0

xyfile=5

#plotting routines

clines=301 #levels of filled contours
ibmax=len(rho[fgfile])+1

if plotqlin==1:
  
    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(R[files[0]][:,4],q1,color='r',label='surfmn',lw=3)
    ax.plot(R[files[0]][:,4],qcon1[:,4],color='b',label='Fluxgrid',lw=3)    

    ax.legend(loc=3,fontsize=20)

    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$R\,({\rm m})$',fontsize=24)
    ax.set_ylabel(r'$q$',fontsize=24)
    
    plt.savefig('q.png',bbox_inches='tight')

if plotmag2wall==1:

    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho1,bnm[xyfile,9,:].real,color='m',label='m=-1',lw=3)
    ax.plot(rho1,bnm[xyfile,8,:].real,color='r',label='m=-2',lw=3)
    ax.plot(rho1,bnm[xyfile,7,:].real,color='b',label='m=-3',lw=3)
    ax.plot(rho1,bnm[xyfile,6,:].real,color='g',label='m=-4',lw=3)
    ax.plot(rho1,bnm[xyfile,5,:].real,color='y',label='m=-5',lw=3)
    ax.plot(rho1,bnm[xyfile,4,:].real,color='lime',label='m=-6',lw=3)
    ax.plot(rho1,bnm[xyfile,11,:].real,color='cyan',label='m=1',lw=3)
    ax.axvline(x=rho1[irho12],lw=3,ls=(0,(3,2)),c='r',label=r'$q=-2$')
    ax.axvline(x=rho1[irho13],lw=3,ls=(0,(3,2)),c='b',label=r'$q=-3$')
    ax.axvline(x=rho1[irho14],lw=3,ls=(0,(3,2)),c='g',label=r'$q=-4$')
    ax.axvline(x=rho1[irho15],lw=3,ls=(0,(3,2)),c='y',label=r'$q=-5$')
    ax.axvline(x=rho1[irho16],lw=3,ls=(0,(3,2)),c='lime',label=r'$q=-5$')
    
    ax.legend(loc=1,ncol=2,fontsize=14)
    
    ax.yaxis.major.formatter._useMathText = True
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.locator_params(axis='x',nbins=5)

#    ax.set_xlim([.1,1]) 

#    ax.set_ylim([0,4e-3])    
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$<r_m>$',fontsize=24)
    ax.set_ylabel(r'$B_m$',fontsize=24)
    
    plt.savefig('Bm2wall16.png',bbox_inches='tight')


if plotmag==1:

    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho[fgfile],bnm[xyfile,9,1:ibmax].real,color='m',label='m=-1',lw=3)
    ax.plot(rho[fgfile],bnm[xyfile,8,1:ibmax].real,color='r',label='m=-2',lw=3)
    ax.plot(rho[fgfile],bnm[xyfile,7,1:ibmax].real,color='b',label='m=-3',lw=3)
    ax.plot(rho[fgfile],bnm[xyfile,6,1:ibmax].real,color='g',label='m=-4',lw=3)
    ax.plot(rho[fgfile],bnm[xyfile,11,1:ibmax].real,color='cyan',label='m=1',lw=3)
    ax.axvline(x=rho[fgfile][irho2],lw=3,ls=(0,(3,2)),c='r',label=r'$q=2$')
    ax.axvline(x=rho[fgfile][irho3],lw=3,ls=(0,(3,2)),c='b',label=r'$q=3$')
    ax.axvline(x=rho[fgfile][irho4],lw=3,ls=(0,(3,2)),c='g',label=r'$q=4$')
    
    ax.legend(loc=1,ncol=2,fontsize=14)
    
    ax.yaxis.major.formatter._useMathText = True
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.locator_params(axis='x',nbins=5)

#    ax.set_xlim([.1,1]) 

#    ax.set_ylim([0,4e-3])    
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$\rho$',fontsize=24)
    ax.set_ylabel(r'$B_m$',fontsize=24)
    
    plt.savefig('Bm16.png',bbox_inches='tight')

if plotsurf2wall==1:
    bmmax=np.amax(bnm[xyfile,:,1:])
    bmmin=0

#    bmmax=np.amax(bm2) 
#    bmmin=np.amin(bm2)    
            
#    if abs(bmmax)>abs(bmmin):
#        bmmin=-bmmax
#    else:
#        bmmax=-bmmin
    
#    bmmax=0.002792068   
    
    nlev=100
    levels=np.arange(bmmin,bmmax,(bmmax-bmmin)/nlev)
    
    fig,ax=plt.subplots(figsize=(10,6))
    
    CS = ax.contourf(m,rho1[1:],np.rot90(bnm[xyfile,:,1:],k=-1),levels,cmap=cm.nipy_spectral)
#    CS = ax.contourf(m,rho[fgfile],bm2,levels,cmap=cm.seismic)
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    cbar=fig.colorbar(CS)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(20)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()
    cbar.set_label(r'$B_m$',fontsize=24)
    ax.locator_params(axis='y',nbins=5)
    
    ax.plot(q1[1:irho110],rho1[1:irho110],c='white',lw=5,ls=(0,(10,2)),label=r'$m=qn$')
    
    ax.axhline(y=rho1[irho12],lw=3,ls=(0,(3,2)),c='r',label=r'$q=-2$')
    ax.axhline(y=rho1[irho13],lw=3,ls=(0,(3,2)),c='b',label=r'$q=-3$')
    ax.axhline(y=rho1[irho14],lw=3,ls=(0,(3,2)),c='g',label=r'$q=-4$')
    ax.axhline(y=rho1[irho15],lw=3,ls=(0,(3,2)),c='y',label=r'$q=-5$')
    ax.axhline(y=rho1[irho16],lw=3,ls=(0,(3,2)),c='lime',label=r'$q=-6$')
        
    ax.legend(loc=4,fontsize=18,framealpha=.75)
    
    ax.set_ylabel(r'$<r_m>$',fontsize=24)
    ax.set_xlabel(r'$m$',fontsize=24)
    
    ax.set_ylim([.1,1.45])    
        
    #plt.savefig('surfmn_comp_15000.png',bbox_inches='tight')
    plt.savefig('surfmn2wall00.png',bbox_inches='tight')

if plotsurf==1:
    bmmax=np.amax(bnm[xyfile,:,1:])
    bmmin=0

#    bmmax=np.amax(bm2) 
#    bmmin=np.amin(bm2)    
            
#    if abs(bmmax)>abs(bmmin):
#        bmmin=-bmmax
#    else:
#        bmmax=-bmmin
    
    bmmax=0.00278744622   
    
    nlev=100
    levels=np.arange(bmmin,bmmax,(bmmax-bmmin)/nlev)
    
    fig,ax=plt.subplots(figsize=(10,6))
    
    CS = ax.contourf(m,rho[fgfile],bnm[xyfile,:,1:ibmax],levels,cmap=cm.nipy_spectral)
#    CS = ax.contourf(m,rho[fgfile],bm2,levels,cmap=cm.seismic)
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    cbar=fig.colorbar(CS)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(20)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()
    cbar.set_label(r'$B_m$',fontsize=24)
    ax.locator_params(axis='y',nbins=5)
    
    ax.plot(q1[:479],rho[fgfile],c='white',lw=5,ls=(0,(10,2)),label=r'$m=qn$')
    
    ax.axhline(y=rho[fgfile][irho2],lw=3,ls=(0,(3,2)),c='r',label=r'$q=-2$')
    ax.axhline(y=rho[fgfile][irho3],lw=3,ls=(0,(3,2)),c='b',label=r'$q=-3$')
    ax.axhline(y=rho[fgfile][irho4],lw=3,ls=(0,(3,2)),c='g',label=r'$q=-4$')
    
    ax.legend(loc=4,fontsize=18,framealpha=.75)
    
    ax.set_ylabel(r'$\rho$',fontsize=24)
    ax.set_xlabel(r'$m$',fontsize=24)
    
    ax.set_ylim([.1,.99])    
        
    #plt.savefig('surfmn_comp_15000.png',bbox_inches='tight')
    plt.savefig('surfmn16.png',bbox_inches='tight')

if plotFSArea==1:
                        
    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho[fgfile],FSArea[1:256],color='r')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Flux Surface Area ${\rm m^2}$')

if plotqcompare==1:
                        
    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho[fgfile],q1[1:256],color='r',label='Mine')
    ax.plot(rho[fgfile],qg[fgfile],color='b',label='Fluxgrid')
    ax.legend()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$q$')

if plottheta==1:

    thmax=np.max(theta_str[1:,:])
    thmin=np.min(theta_str[1:,:])
    
#    if abs(thmax)>abs(thmin):
#        thmin=-thmax
#    else:
#        thmax=-thmin
    
    levels=np.arange(thmin,thmax,(thmax-thmin)/clines)
    
    fig,ax=plt.subplots()
    
    CS = ax.contourf(R[files[0]],Z[files[0]],theta_str,clines,cmap=cm.nipy_spectral)
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_aspect('equal')
    
    cbar=fig.colorbar(CS)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(20)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()
    cbar.set_label(r'$\Theta$',fontsize=24)

    ax.set_xlabel(r'$R\,{\rm (m)}$',fontsize=24)
    ax.set_ylabel(r'$Z\,{\rm (m)}$',fontsize=24)

    plt.savefig('Theta_straight.png',bbox_inches='tight')


if plotjac==1:
    pr_max = np.max(jac)
    pr_min = np.min(jac)
    
    #if abs(pr_max)>abs(pr_min):
    #    pr_min=-pr_max
    #else:
    #    pr_max=-pr_min
    
    lev,delta=np.linspace(pr_min,pr_max,clines,retstep=True)      
                                                    
    fig,ax=plt.subplots(figsize=(5,6))                
    CS=ax.contourf(R[files[0]],Z[files[0]],jac,levels=lev,cmap=cm.gnuplot2)
    
    ax.set_aspect('equal')
    
    plt.setp(ax.get_xticklabels(), fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=18)
    #plt.setp(ax.get_yticklabels(), fontsize=12)
    #ax.set_xlabel('$\sqrt{\psi} $',fontsize=30)
    #ax.set_ylabel('$q $',fontsize=30,rotation='horizontal')
    ax.set_title(r'Jacobian',fontsize=22)
    ax.set_xlabel(r'$R\,{\rm(m)}$',fontsize=20)
    ax.set_ylabel(r'$ Z\,{\rm(m)}$',fontsize=20)
    #cbar=plt.colorbar(CS,pad=0.01,format="%0.000f")
    cbar=plt.colorbar(CS,pad=0.01)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(16)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()
    
    plt.savefig('xy_jac_00000.png',bbox_inches='tight')      

if plotqcon==1:

    qcon=np.zeros((len(R[files[0]]),len(R[files[0]][0])))

    for i in range(len(qcon)):
        for j in range(len(qcon[0])):
            qcon[i,j]=q1[i]

    pr_max = np.max(qcon)
    pr_min = np.min(qcon)
    
    #if abs(pr_max)>abs(pr_min):
    #    pr_min=-pr_max
    #else:
    #    pr_max=-pr_min
    
    lev,delta=np.linspace(pr_min,pr_max,clines,retstep=True)
    
    qlev=np.array([-4,-3,-2])        
                                                    
    fig,ax=plt.subplots(figsize=(5,6))                
    CS=ax.contourf(R[files[0]],Z[files[0]],qcon,levels=lev,cmap=cm.gnuplot2)
    
    QS=ax.contour(R[files[0]],Z[files[0]],qcon,levels=qlev,colors='k',linestyles='solid')
    
    fmt = {}
    strs=[r'$q=\!-4$',r'$q=\!-3$',r'$q=\!-2$']
    for l, s in zip(QS.levels, strs):
        fmt[l] = s
    
    manual_locations=[(1.36,0.85),(1.42,0.73),(1.48,0.59)]
    
    ax.clabel(QS,qlev,inline=1,fmt=fmt,manual=manual_locations,inline_spacing=15)
    
    ax.set_aspect('equal')
    
    plt.setp(ax.get_xticklabels(), fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=18)
    #plt.setp(ax.get_yticklabels(), fontsize=12)
    #ax.set_xlabel('$\sqrt{\psi} $',fontsize=30)
    #ax.set_ylabel('$q $',fontsize=30,rotation='horizontal')
    ax.set_title(r'Safety Factor',fontsize=22)
    ax.set_xlabel(r'${\rm R\,(m)}$',fontsize=20)
    ax.set_ylabel(r'${\rm Z\,(m)}$',fontsize=20)
    #cbar=plt.colorbar(CS,pad=0.01,format="%0.000f")
    cbar=plt.colorbar(CS,pad=0.01)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.yaxis.set_offset_position('left')
    cbar.ax.yaxis.offsetText.set_fontsize(16)
    cbar.formatter.set_powerlimits((0,0))
    cbar.update_ticks()
    
    plt.savefig('xy_q_00000.png',bbox_inches='tight')   

plt.show()