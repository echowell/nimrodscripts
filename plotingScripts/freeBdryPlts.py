import numpy as np
import matplotlib.pyplot as plt
import os

homeDir = os.environ['HOME']

fname=homeDir + "SCRATCH/NTM_TESTING/18102201/fbdry.txt"
data = np.loadtxt(fname)

ncoils = int(max(data[:,1]))
its = data.shape[0]/ncoils
print its


step = np.zeros(its)
ic0 = np.zeros([its,ncoils])
ic = np.zeros([its,ncoils])
chisq = np.zeros(its)

for ii in range(0, data.shape[0], ncoils):
    it = int(ii/ncoils)
    step[it] = data[ii,0]
    chisq[it] = data[ii,6]
    for jj in range(ncoils):
        ic0[it,jj] = data[ii+jj,2]/1e6
        ic[it,jj] = data[ii+jj,5]/1e6
        
ict=ic0+ic


maxplt = 100                
fig, ax1 = plt.subplots()
ax1.plot(step[0:maxplt], chisq[0:maxplt]/100, 'b-')
ax1.set_xlabel(r'GS Iteration',fontsize=16)
#ax1.axes.set_xlim(left=0,right=1)

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$\chi^2$',rotation=90,fontsize=16)
#ax1.axes.set_ylim(0,1)
plt.yscale("log")
ax1.tick_params(axis='both', which='major', labelsize=14)
#plt.locator_params(axis='y', nbins=6)
fig.tight_layout()

plt.show()      

#inital current

fig, ax1 = plt.subplots()
ax1.plot(step[0:maxplt], ic0[0:maxplt,0], '-')
#ax1.plot(step, ic0[:,1], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,2], '-')
#ax1.plot(step, ic0[0:maxplt,3], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,4], '-')
#ax1.plot(step, ic0[:,5], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,6], '-')
#ax1.plot(step, ic0[:,7], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,8], '-')
#ax1.plot(step, ic0[:,9], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,10], '-')
#ax1.plot(step, ic0[:,11], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,12], '-')
#ax1.plot(step, ic0[:,13], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,14], '-')
#ax1.plot(step, ic0[:,15], '-')
ax1.plot(step[0:maxplt], ic0[0:maxplt,16], '-')
#ax1.plot(step, ic0[:,17], '-')

ax1.set_xlabel(r'GS Iteration',fontsize=16)
#ax1.axes.set_xlim(left=0,right=1)

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'Inital $I_c$ [MA]',rotation=90,fontsize=16)
#ax1.axes.set_ylim(0,1)
ax1.tick_params(axis='both', which='major', labelsize=14)
plt.locator_params(axis='y', nbins=6)
fig.tight_layout()

plt.show()      

#inital current

fig, ax1 = plt.subplots()
ax1.plot(step[0:maxplt], ic[0:maxplt,0], '-')
#ax1.plot(step, ic[:,1], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,2], '-')
#ax1.plot(step, ic[:,3], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,4], '-')
#ax1.plot(step, ic[:,5], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,6], '-')
#ax1.plot(step, ic[:,7], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,8], '-')
#ax1.plot(step, ic[:,9], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,10], '-')
#ax1.plot(step, ic[:,11], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,12], '-')
#ax1.plot(step, ic[:,13], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,14], '-')
#ax1.plot(step, ic[:,15], '-')
ax1.plot(step[0:maxplt], ic[0:maxplt,16], '-')
#ax1.plot(step, ic[:,17], '-')

ax1.set_xlabel(r'GS Iteration',fontsize=16)
#ax1.axes.set_xlim(left=0,right=1)

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$\delta I_c$ [MA]',rotation=90,fontsize=16)
#ax1.axes.set_ylim(0,1)
ax1.tick_params(axis='both', which='major', labelsize=14)
plt.locator_params(axis='y', nbins=6)
fig.tight_layout()

plt.show()      

#totla current

fig, ax1 = plt.subplots()
ax1.plot(step[0:maxplt], ict[0:maxplt,0], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,2], '-')
#ax1.plot(step, ic[:,3], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,4], '-')
#ax1.plot(step, ic[:,5], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,6], '-')
#ax1.plot(step, ic[:,7], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,8], '-')
#ax1.plot(step, ic[:,9], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,10], '-')
#ax1.plot(step, ic[:,11], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,12], '-')
#ax1.plot(step, ic[:,13], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,14], '-')
#ax1.plot(step, ic[:,15], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,16], '-')
#ax1.plot(step, ic[:,17], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,1], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,3], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,5], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,7], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,9], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,11], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,13], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,15], '-')
ax1.plot(step[0:maxplt], ict[0:maxplt,17], '-')


ax1.set_xlabel(r'GS Iteration',fontsize=16)
#ax1.axes.set_xlim(left=0,right=1)

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r'$I_c$ [MA]',rotation=90,fontsize=16)
#ax1.axes.set_ylim(0,1)
ax1.tick_params(axis='both', which='major', labelsize=14)
plt.locator_params(axis='y', nbins=6)
fig.tight_layout()

plt.show()      