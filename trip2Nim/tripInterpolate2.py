''' Intepolate B(R,Z,PHI) from one set of RZ points onto a second set.

'''


import sys
sys.path.insert(0, "./")
import os
import numpy as np
import matplotlib.pyplot as plt

home_dir = os.environ['HOME']
scratch_dir = home_dir + "/SCRATCH"
old_file_dir = scratch_dir + "/nimruns/echowell_runs/heatwidthscaling/166439/03300/EF_GRID_18121801"
probe_gb_file = old_file_dir + "/166439.03300.probe_gb.out"
new_rz_dir = scratch_dir + "/166439/03300_2_equilbria/19061201"
new_rz_file = new_rz_dir + "/nimrod_bdry_rz.txt"

b_data = np.loadtxt(probe_gb_file,comments='%')
new_rz = np.loadtxt(new_rz_file,skiprows=1,delimiter=',')

phi_RZ_slices = []
phi_B_slices = []
phi_values = []
phi_last = -100
first_phi = True
for ii in range(b_data.shape[0]):
  if first_phi:
    first_phi = False
    start_index = ii
    phi_last = b_data[ii,0]
    count=1
    phi_values.append(phi_last)
  elif b_data[ii,0] == phi_last:
    count+=1
  else:
    thisRZ = np.zeros([count-1,2])
    thisB = np.zeros([count-1,3])
    thisRZ[:,0] = b_data[start_index:ii-1,1]
    thisRZ[:,1] = b_data[start_index:ii-1,2]
    thisB[:,0] = b_data[start_index:ii-1,4] #BR
    thisB[:,1] = b_data[start_index:ii-1,5] #BZ
    thisB[:,2] = b_data[start_index:ii-1,3] #Bphi
    phi_RZ_slices.append(thisRZ)
    phi_B_slices.append(thisB)
    phi_values.append(phi_last)
    start_index = ii
    phi_last = b_data[ii,0]
    count=1

ii=b_data.shape[0]
thisRZ = np.zeros([count-1,2])
thisB = np.zeros([count-1,3])
thisRZ[:,0] = b_data[start_index:ii-1,1]
thisRZ[:,1] = b_data[start_index:ii-1,2]
thisB[:,0] = b_data[start_index:ii-1,4] #BR
thisB[:,1] = b_data[start_index:ii-1,5] #BZ
thisB[:,2] = b_data[start_index:ii-1,3] #Bphi
phi_RZ_slices.append(thisRZ)
phi_B_slices.append(thisB)
phi_values.append(phi_last)
    
print(count)
print(phi_values)
print(b_data.shape)


for rz in phi_RZ_slices:
  plt.plot(rz[:,0],rz[:,1])
plt.show()