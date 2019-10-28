#!/usr/bin/env python3

''' Intepolate B(R,Z,PHI) from one set of RZ points onto a second set.

    Both set's of RZ points are assumed to be on the same surface.
    First: Read NIMROD wall file and the RZ points are Fouier Transformed
    Second: Read TRIP3D files, and Fourier transform B
    Third: Read new RZ locations, and find theta of this coordinates
    Finally: Interpolate B onto new coordinates
'''

import sys
sys.path.insert(0, "./")
import os
import nim_wall as wall
import numpy as np
import matplotlib.pyplot as plt

home_dir = os.environ['HOME']
scratch_dir = home_dir + "/SCRATCH"
file_dir = scratch_dir + "/166439/03300_2_equilbria/19061201"
wall_file =  "wall_file_2.txt"
full_wall_file_path = file_dir + "/" + wall_file
new_rz_file = "nimrod_bdry_rz.txt"
full_new_rz_file_path = file_dir + "/" + new_rz_file
full_wall_file_path = "/home/research/ehowell/SCRATCH/166439/03300_q104_flowtesting/flow_60/nimrod_bdry_rz.txt"
wall_modes = -1

this_wall = wall.NimrodWall(full_wall_file_path, wall_modes)

new_rz = np.loadtxt(full_new_rz_file_path, delimiter=',', skiprows=1)
#old_rz = np.loadtxt(full_old_rz_file_path, delimiter=',', skiprows=1)
new_theta = np.zeros(new_rz.shape[0])
err_theta = np.zeros(new_rz.shape[0])
for ix in range(new_rz.shape[0]):
  new_theta[ix], err_theta[ix] = this_wall.get_theta_rz(new_rz[ix,0],new_rz[ix,1])
print (new_theta)
print (err_theta)


plt.plot(this_wall.file_rz[:,0],this_wall.file_rz[:,1], 'bo', label="Orginal Wall")
plt.plot(this_wall.rz_wall[:,0],this_wall.rz_wall[:,1], 'g-', label="Smooth Wall")
plt.plot(new_rz[:,0],new_rz[:,1],'k', label="New Rz")
#plt.plot(old_rz[:,0],old_rz[:,1],'g', label="Old Rz")
plt.legend(loc=1)
plt.show()