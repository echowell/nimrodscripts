#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import scipy.interpolate as interp


def get_m_index(m, m_max):
  return m+m_max

plot_type = "nonlinear2"

#vacuum,linear, nonlinear


run_dir = "19092601_l5"
#run_dir = "1908163"
dump_num = "40000"
homeDir = os.environ['HOME']
scratchDir = homeDir + '/SCRATCH'
fileName = scratchDir+'/174446_novac_fl/eq26/1908013/08000/surfmn.08000.h5'
fileName = scratchDir+'/174446_novac_debug/vac_eq28_rmp/surfmn.00100.h5'
fileName = "/home/research/ehowell/SCRATCH/174446_novac_new_eq/surfmn_optimize/vac/surfmn.00100.h5"
fileName = "/home/research/ehowell/SCRATCH/174446_novac_new_eq/surfmn_optimize/vac/surfmn.00100.h5"
fileName = scratchDir+'/174446_novac_fl/eq26/'+run_dir+'/'+ dump_num +'/surfmn.'+ dump_num + '.h5'

#run_dir = "case_15_vac"
#dump_num = "00100"
#fileName = scratchDir+'/174446_novac_new_eq/surfmn_scan/'+run_dir+'/surfmn.'+ dump_num + '.h5'

#fileName = scratchDir+'/166439/03300_vac_eq/complexconj_rmp_vac5_fpsep2/surfmn.00005.h5'
#fileName = scratchDir+'/166439/03300_vac_eq/complexconj_rmp_vac5_fpsep2/surfmn.00005_rr.h5'
#fileName = scratchDir+'/166439//03300_vac_eq/complexconj_rmp_vac/surfmn.00005.h5'

if (plot_type=="vacuum"):
  fileName = scratchDir + "/166439/03300_2_equilbria/19091201_vac_lphi5_fp/surfmn.00100.h5"
  n1_title = r' Vacuum n=1 $|B_\psi|$ [A.U.]'
  n3_title = r' Vacuum n=3 $|B_\psi|$ [A.U.]'
elif (plot_type=="linear"):
  fileName = scratchDir + "/166439/03300_2_fl/19091702/lphi5_rmp_cfl_b/200000/surfmn.200000.h5"
  n1_title = r' Linear n=1 $|B_\psi|$ [A.U.]'
  n3_title = r' Linear n=3 $|B_\psi|$ [A.U.]'
elif (plot_type=="nonlinear"):
  fileName = scratchDir + "/166439/03300_2_fl//19091702/lphi5_nolinear_fresh/22000/surfmn.22000.h5"
  n1_title = r' Nonlinear n=1 $|B_\psi|$ [A.U.]'
  n3_title = r' Nonlinear n=3 $|B_\psi|$ [A.U.]'
elif (plot_type=="nonlinear2"):
  fileName = scratchDir + "/166439/03300_2_fl//19091702/lphi5_nolinear_fresh/10000/surfmn.10000.h5"
  n1_title = r' Nonlinear n=1 $|B_\psi|$ [A.U.]'
  n3_title = r' Nonlinear n=3 $|B_\psi|$ [A.U.]'

#n1_scale
plt.ylabel(r'$|B_\psi|$ [A.U.]')
sfac=1e4
vmax_n1 = 4.0
vmax_n3 = 4.0
level1=np.linspace(0,vmax_n1,301)
level3=np.linspace(0,vmax_n3,301)

n1_cbar_ticks = np.linspace(0,vmax_n1,11)
n3_cbar_ticks = np.linspace(0,vmax_n3,11)
print (n1_cbar_ticks)

m_plot_list = [-6,-5,-4,-3,-2,-1, 0]
q1_plot_list = [-2,-3,-4]
q3_plot_list = [-1.67,-2.00,-2.33]

profNames = ['Vprime','q']
with h5py.File(fileName,'r') as fc:
#  for aname, avalue in fc.attrs.items():
#    print(aname,avalue)
  mrGrid = fc['surfmnGrid'][:]
  bmn = fc['Bmn001'][:]
  bmn3 = fc['Bmn003'][:]
  rho = fc['rho'][:]
  profs = fc['prof'][:]
#  print(fc.keys())
  
psi_of_q = interp.interp1d(profs[1,:],rho)
m1_range = np.linspace(-1.043,-5.9)
q_of_m1 = psi_of_q(m1_range)
m3_range = np.linspace(-3.129,-15)
q_of_m3 = psi_of_q(m3_range/3.00)

m_max = int((bmn.shape[1]-1)/2)
#print(m_max)

#print(mrGrid.shape[2])
#print(bmn.shape)
fig = plt.figure()
ax=fig.add_subplot(111)
plt.set_cmap('nipy_spectral')
conf=plt.contourf(mrGrid[0,:,:],mrGrid[1,:,:],bmn*sfac,levels=level1,vmax=vmax_n1)
plt.plot(m1_range,q_of_m1,c='w')
plt.title(n1_title,fontsize=16)
plt.ylabel(r'<r>',fontsize=16)
plt.xlabel('m',fontsize=16)
cbar=fig.colorbar(conf,ticks=n1_cbar_ticks)
plt.show()
#plt.colorbar(ticks=n1_cbar_ticks)

fig = plt.figure()
ax=fig.add_subplot(111)
plt.set_cmap('nipy_spectral')
conf=plt.contourf(mrGrid[0,:,:],mrGrid[1,:,:],bmn3*sfac,levels=level3,vmax=vmax_n3)
plt.plot(m3_range,q_of_m3,c='w')
plt.title(n3_title,fontsize=16)
plt.ylabel(r'<r>',fontsize=16)
plt.xlabel('m',fontsize=16)
cbar=fig.colorbar(conf,ticks=n3_cbar_ticks)
plt.show()

fig = plt.figure(figsize=(6,6))
ax= fig.add_subplot(111)
for this_m in m_plot_list:
  this_i = get_m_index(this_m,m_max)
  plt_lbl = "m = " + str(this_m)
  ax.plot(mrGrid[1,:,1],bmn[:,this_i], label=plt_lbl)
for this_q in q1_plot_list:
  this_rho = psi_of_q(this_q) 
  this_lbl = "q = " + str(this_q)
  ax.axvline(this_rho,ls=':', label=this_lbl)

#ax.axvline(0.607025,ls=':',c='b', label="q=2")
#ax.axvline(0.75138,ls=':',c='g', label="q=3")
#ax.axvline(0.849892,ls=':',c='c', label="q=4")

ax.axhline(0,ls='-',c='k')
ax.legend(loc=0)
plt.title(r'Vacuum n=1 response')
plt.xlabel(r'<r>')
plt.ylabel(r'$|B_\psi|$ [A.U.]')
plt.tight_layout()
plt.show()


fig = plt.figure(figsize=(6,6))
ax= fig.add_subplot(111)
for this_m in m_plot_list:
  this_i = get_m_index(this_m,m_max)
  plt_lbl = "m = " + str(this_m)
  ax.plot(mrGrid[1,:,1],bmn3[:,this_i], label=plt_lbl)
for this_q in q3_plot_list:
  this_rho = psi_of_q(this_q) 
  this_lbl = "q = " + str(this_q)
  ax.axvline(this_rho,ls=':', label=this_lbl)

#ax.axvline(0.607025,ls=':',c='b', label="q=2")
#ax.axvline(0.75138,ls=':',c='g', label="q=3")
#ax.axvline(0.849892,ls=':',c='c', label="q=4")

ax.axhline(0,ls='-',c='k')
ax.legend(loc=0)
plt.title(r'Vacuum n=3 response')
plt.xlabel(r'<r>')
plt.ylabel(r'$|B_\psi|$ [A.U.]')
plt.tight_layout()
plt.show()