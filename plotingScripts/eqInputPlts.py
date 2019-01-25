import numpy as np
import matplotlib.pyplot as plt
import os


homeDir = os.environ['HOME']

fname=home + "/SCRATCH/NTM_TESTING/18102204_refine/eq_input_python.out"


data = np.loadtxt(fname)


rho = np.sqrt(data[:,1])
print rho
q = data[:,4]
q[0] = q[1]
f = data[:,2]
p = data[:,3]
print q



fig, ax1 = plt.subplots()
ax1.plot(rho, q, 'b-')
ax1.set_xlabel(r'$\rho$',fontsize=16)
ax1.axes.set_xlim(left=0,right=1)

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Safety Factor',rotation=90,fontsize=16)
ax1.axes.set_ylim(0,11)

ax2 = ax1.twinx()
ax2.plot(rho, p, 'r-')
ax2.set_ylabel('Pressure (Tm)',rotation=90,fontsize=16)
ax2.axes.set_ylim(0,0.14)

ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.locator_params(axis='y', nbins=6)
fig.tight_layout()

plt.show()