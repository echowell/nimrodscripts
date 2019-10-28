#!/usr/bin/env python3
# Input files:
# Ouput file:


import os
import shutil
import biotSavartFunction as bsf
import scipy.optimize as opt
import numpy as np

directory = os.getcwd()
dir_list=["gs","vac"]

home_dir = os.environ['HOME']

fgnimeq_cmd = home_dir + "/SCRATCH/nimall/builds/nimdevel/ser/bin/fgnimeq-ser"
nimrod_cmd = home_dir + "/SCRATCH/nimall/builds/nimdevel/ser/bin/nimrod-ser"
nimplot_cmd = home_dir + "/SCRATCH/nimall/builds/nimdevel/ser/bin/nimplot < nimplot_inputs"

gs_files=["a174446.03390","fluxgrid.in","g174446.03390","nimeq.in","nimrod.in",
"nimrod_bdry_rz.txt","oculus.in","p174446.3390.0_new_rot_fits","p174446.3390.0_new_rot_fits.smooth"]
vac_files=["nimrod.in_vac","nimplot_inputs"]
src_dir="nimrod_files"
files_dict={"gs":gs_files, "vac":vac_files}


def run_nimrod():
  os.chdir("gs")
  os.system(fgnimeq_cmd)
  os.chdir("..")
  src_path="gs/dumpgll.00000.h5"
  cp_path ="vac/dumpgll.00000.h5"
  shutil.copy(src_path, cp_path)
  os.chdir("vac")
  os.system(nimrod_cmd)
  os.system(nimplot_cmd)
  os.chdir("..")

def wrapper_fun(values):
  ncoil=int(len(values)/2)
  coil_theta=[]
  coil_tilt=[]
  distance_coil=1.5
  coil_r=0.6
  for ii in range(ncoil):    
    coil_theta.append(values[ii])
    coil_tilt.append(values[ii+ncoil])
  bsf.coil_calc("gs/",ncoil,coil_theta,coil_tilt,distance_coil,coil_r)
  run_nimrod()
  fileName="vac/surfmn.00100.h5"
  cost=bsf.surfmn_eval(fileName)
  with open("runlog.txt",'a+') as this_file:
    this_file.write("Status at step \n")
    this_file.write("Coil theta \n")
    this_file.write(str(coil_theta))
    this_file.write("\n Coil tilt \n")
    this_file.write(str(coil_tilt))
    this_file.write("\n Cost function \n")
    this_file.write(str(cost))
    this_file.write(" \n")
  return cost

for this_dir in dir_list:
  try:
    os.mkdir(this_dir)
  except FileExistsError:
    print("Folder "+ this_dir + " exists, removing contents")
    for this_file in os.listdir(this_dir):
      file_path = os.path.join(this_dir, this_file)
      os.remove(file_path)
  this_files=files_dict[this_dir]
  for this_file in this_files:
    src_path = os.path.join(src_dir,this_file)
    if (this_dir=="vac" and this_file=="nimrod.in_vac"):
      cp_path = os.path.join(this_dir,"nimrod.in")
    else:
      cp_path = os.path.join(this_dir,this_file)
    shutil.copy(src_path, cp_path)

x0 = [0, .15, .5, .85, 0.25, .65, .75, 0.85]

bnds = ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
tol=1.0e-6
res = opt.minimize(wrapper_fun,x0,method="TNC",bounds=bnds,tol=tol,options={'disp': True, "eps":1.0e-3, "maxiter":20})
wrapper_fun(res.x)
with open("result.txt",'w') as this_file:
  this_file.write("Optimization complete\n")
  this_file.write("Optimize success?:  ")
  this_file.write(str(res.success))
  this_file.write("\n The optimal parameters are:\n")
  this_file.write(str(res.x))
  this_file.write("\n Optimization iterations :  ")
  this_file.write(str(res.nit))
  this_file.write("\n The optimal cost :  ")
  this_file.write(str(res.fun))

print(res.x)

