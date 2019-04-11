#!/usr/bin/python

import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from surfmn import files,rho1,bnm,irho12,irho13,irho14,irho15,irho16,rho,fgfile,irho2,irho3,irho4,q1,irho110,m,mrange

bnmflip=bnm[:,::-1]

plotsurf=1
plotsurf2wall=0
plotmag=1
plotmag2wall=0
bylim=9e-4

savepath = "/home/research/ehowell/Documents/Talks/Sherwood2019/Sherwood2019_poster/Figures/"

basename="vac1"
nfour=1


xyfile=0

xyfilestr=files[xyfile][8:10]



vacn2=bnmflip[0,8,irho2]
vacn3=bnmflip[0,7,irho3]
vacn4=bnmflip[0,6,irho4]

vacn12=bnmflip[0,8,irho12]
vacn13=bnmflip[0,7,irho13]
vacn14=bnmflip[0,6,irho14]
vacn15=bnmflip[0,5,irho15]
vacn16=bnmflip[0,4,irho16]

#plotting routines

clines=301 #levels of filled contours
ibmax=len(rho[fgfile])+1

if plotmag2wall==1:

    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho1,bnmflip[xyfile,mrange+1].real,color='m',label='m=-1',lw=3)
    ax.plot(rho1,bnmflip[xyfile,mrange+2].real,color='r',label='m=-2',lw=3)
    ax.plot(rho1,bnmflip[xyfile,mrange+3].real,color='b',label='m=-3',lw=3)
    ax.plot(rho1,bnmflip[xyfile,mrange+4].real,color='g',label='m=-4',lw=3)
    ax.plot(rho1,bnmflip[xyfile,mrange+5].real,color='y',label='m=-5',lw=3)
    ax.plot(rho1,bnmflip[xyfile,mrange+6].real,color='lime',label='m=-6',lw=3)
    ax.axvline(x=rho1[irho12],lw=3,ls='dotted',c='r',label=r'$q=-2$')
    ax.axvline(x=rho1[irho13],lw=3,ls='dotted',c='b',label=r'$q=-3$')
    ax.axvline(x=rho1[irho14],lw=3,ls='dotted',c='g',label=r'$q=-4$')
    ax.axvline(x=rho1[irho15],lw=3,ls='dotted',c='y',label=r'$q=-5$')
    ax.axvline(x=rho1[irho16],lw=3,ls='dotted',c='lime',label=r'$q=-5$')

    ax.plot(rho1[irho12],vacn12,marker=7,color='r',markersize=16)    
    ax.plot(rho1[irho13],vacn13,marker=7,color='b',markersize=16)  
    ax.plot(rho1[irho14],vacn14,marker=7,color='g',markersize=16)  
    ax.plot(rho1[irho15],vacn15,marker=7,color='y',markersize=16)  
    ax.plot(rho1[irho16],vacn16,marker=7,color='lime',markersize=16)
                        
    ax.legend(loc=1,ncol=2,fontsize=14)
    
    ax.yaxis.major.formatter._useMathText = True
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.locator_params(axis='x',nbins=5)

#    ax.set_xlim([.1,1]) 

#    ax.set_ylim([0,4e-3])    
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$<r>_N}$',fontsize=24)
    ax.set_ylabel(r'$B_{r(m,n)}\,{\rm (T)}$',fontsize=24)
    
    plt.savefig('Bm2wall'+xyfilestr+'.png',bbox_inches='tight')


if plotmag==1:

    fig,ax=plt.subplots(figsize=(6,6))
    ax.plot(rho[fgfile],bnmflip[xyfile,mrange+1,1:ibmax].real,color='m',label='m=-1',lw=3)
    ax.plot(rho[fgfile],bnmflip[xyfile,mrange+2,1:ibmax].real,color='r',label='m=-2',lw=3)
    ax.plot(rho[fgfile],bnmflip[xyfile,mrange+3,1:ibmax].real,color='b',label='m=-3',lw=3)
#    ax.plot(rho[fgfile],bnmflip[xyfile,mrange+4,1:ibmax].real,color='g',label='m=-4',lw=3)
    ax.axvline(x=rho[fgfile][irho2],lw=3,ls='dotted',c='r',label=r'$q=2$')
    ax.axvline(x=rho[fgfile][irho3],lw=3,ls='dotted',c='b',label=r'$q=3$')
  #  ax.axvline(x=rho[fgfile][irho4],lw=3,ls='dotted',c='g',label=r'$q=4$')

#    ax.plot(rho[fgfile][irho2],vacn2,marker=7,color='r',markersize=16)    
#    ax.plot(rho[fgfile][irho3],vacn3,marker=7,color='b',markersize=16)  
#    ax.plot(rho[fgfile][irho4],vacn4,marker=7,color='g',markersize=16)     
            
    ax.legend(loc=2,ncol=2,fontsize=14)
    
    ax.yaxis.major.formatter._useMathText = True
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,-2))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.locator_params(axis='x',nbins=5)

#    ax.set_xlim([.1,1]) 

    ax.set_ylim([0,bylim])    
    
    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    
    ax.set_xlabel(r'$\rho_N$',fontsize=24)
    ax.set_ylabel(r'$B_{r(m,n)}\,{\rm (T)}$',fontsize=24)
    
#    plt.savefig('Bm'+xyfilestr+'.png',bbox_inches='tight')
    
    plt.savefig(savepath + basename + "Bm.png",bbox_inches='tight')
if plotsurf2wall==1:
    bmmax=np.amax(bnmflip[xyfile,:,1:])
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
    
    CS = ax.contourf(m,rho1[1:],np.rot90(bnmflip[xyfile,:,1:],k=-1),levels,cmap=cm.nipy_spectral)
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
    cbar.set_label(r'$B_{r(m,n)}\,{\rm (T)}$',fontsize=24)
    ax.locator_params(axis='y',nbins=5)
    
    ax.plot(q1[1:irho110],rho1[1:irho110],c='white',lw=5,ls='dashed',label=r'$m=qn$')
    
    ax.axhline(y=rho1[irho12],lw=3,ls='dotted',c='r',label=r'$q=-2$')
    ax.axhline(y=rho1[irho13],lw=3,ls='dotted',c='b',label=r'$q=-3$')
#    ax.axhline(y=rho1[irho14],lw=3,ls='dotted',c='g',label=r'$q=-4$')
#    ax.axhline(y=rho1[irho15],lw=3,ls='dotted',c='y',label=r'$q=-5$')
#    ax.axhline(y=rho1[irho16],lw=3,ls='dotted',c='lime',label=r'$q=-6$')
        
    ax.legend(loc=4,fontsize=18,framealpha=.75)
    
    ax.set_ylabel(r'$<r>_N$',fontsize=24)
    ax.set_xlabel(r'$m$',fontsize=24)
    ax.set_xlim([-mrange,mrange])
  
    ax.set_ylim([.1,1.45])    
        
    #plt.savefig('surfmn_comp_15000.png',bbox_inches='tight')
    plt.savefig('surfmn2wall'+xyfilestr+'.png',bbox_inches='tight')

if plotsurf==1:
    bmmax=np.amax(bnmflip[xyfile,:,1:])
    bmmax=bylim
    bmmin=0

#    bmmax=np.amax(bm2) 
#    bmmin=np.amin(bm2)    
            
#    if abs(bmmax)>abs(bmmin):
#        bmmin=-bmmax
#    else:
#        bmmax=-bmmin
    
#    bmmax=0.00278744622   
    
    nlev=100
    levels=np.arange(bmmin,bmmax,(bmmax-bmmin)/nlev)
    
    fig,ax=plt.subplots(figsize=(10,6))
    
    CS = ax.contourf(m,rho[fgfile],np.rot90(bnmflip[xyfile,:,1:ibmax],k=-1),levels,cmap=cm.nipy_spectral)
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
    cbar.set_label(r'$B_{r(m,n)}\,{\rm (T)}$',fontsize=24)
    ax.locator_params(axis='y',nbins=5)
    
    ax.plot(nfour*q1[:rho[fgfile].size],rho[fgfile],c='white',lw=5,ls='dashed',label=r'$m=qn$')
    
    ax.axhline(y=rho[fgfile][irho2],lw=3,ls='dotted',c='r',label=r'$q=-2$')
    ax.axhline(y=rho[fgfile][irho3],lw=3,ls='dotted',c='b',label=r'$q=-3$')
    #ax.axhline(y=rho[fgfile][irho4],lw=3,ls='dotted',c='g',label=r'$q=-4$')
    
    ax.legend(loc=4,fontsize=18,framealpha=.75)
    
    ax.set_ylabel(r'$\rho_N$',fontsize=24)
    ax.set_xlabel(r'$m$',fontsize=24)
    ax.set_xlim([-mrange,mrange])
    ax.set_ylim([.1,.99])    
        
    #plt.savefig('surfmn_comp_15000.png',bbox_inches='tight')
    #plt.savefig('surfmn'+xyfilestr+'.png',bbox_inches='tight')
    plt.savefig(savepath + basename + "surfmn.png",bbox_inches='tight')

plt.show()
