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

this_coils=0

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

def wrapper_fun(currents):
  this_coils.coil_calc(currents)
  run_nimrod()
  fileName="vac/surfmn.00100.h5"
  cost=bsf.surfmn_eval(fileName)
  with open("runlog.txt",'a+') as this_file:
    this_file.write("Status at step \n")
    this_file.write("Coil currents \n")
    this_file.write(str(currents))
    this_file.write("\n Cost function \n")
    this_file.write(str(cost))
    this_file.write(" \n")
  return cost

### start main program

### set up coils
currents=[0.47528571, -0.57160186,  0.6809352,  -0.13174283]
this_coils=bsf.coil_opt("./")
this_coils.coil_calc(currents)

