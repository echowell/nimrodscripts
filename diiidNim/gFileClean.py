#! /usr/bin/env python3
''' This file reads an gFile, allows the user to modify the F and p profiles,
    and then writes a new gFile. This is useful for cleaning up the profiles 
    near the axis.
'''

import os
import datetime 
import numpy as np
import numpy as np
from scipy.interpolate import splev,splrep
import matplotlib.pyplot as plt


def readWriteEfitHeader(read_file, write_file):
  '''This function reads are writes the EFIT File Header 
     and returns the mesh mx and my''' 
  this_line = f.readline()
  this_words = this_line.split()
  efit_id = this_words[0]
  today = datetime.date.today().strftime("%m/%d/%Y") 
  shot_num = this_words[2]
  recon_time = this_words[3]
  no_clue = int(this_words[4])
  mx = int(this_words[5])
  my = int(this_words[6])
  header = "  {0:<9s}{1:<14s}{2:<9s}{3:<17s}{4:<2d}{5:<4d}{6:<4d}\n".format(efit_id,today,shot_num,recon_time,no_clue,mx,my)
  write_file.write(header)
  for ii in range(4):
    this_line = f.readline()
    write_file.write(this_line)
  return mx, my

def readWriteRestFile(read_file,write_file):
  while True:
    this_line=f.readline()
    if (len(this_line)==0): break
    write_file.write(this_line)

def readProfile(read_file, ma):
  profile = np.zeros(ma)
  floats_per_line = 5
  lines = ma // floats_per_line
  float_length=16
  for ii in range(0,lines):
    this_line = f.readline()
    for jj in range(0,floats_per_line):
      profile[ii*floats_per_line + jj] = float(this_line[jj*16:(jj+1)*16])
  this_line = f.readline()
  ii = lines
  for jj in range(0, ma % floats_per_line):
    profile[ii*floats_per_line + jj] = float(this_line[jj*16:(jj+1)*16])
  return profile

def writeProfile(write_file, profile):
  floats_per_line = 5  
  float_length=16
  for ii in range(profile.size//floats_per_line):
    this_line=''
    for jj in range(floats_per_line):
      this_line += '{: 16.9E}'.format(profile[ii*floats_per_line+jj])
    write_file.write(this_line+'\n')
  this_line=''
  ii = profile.size//floats_per_line
  for jj in range(0, profile.size % floats_per_line):
      this_line += '{: 16.9E}'.format(profile[ii*floats_per_line+jj])
  write_file.write(this_line+'\n')


class modField:
  def __init__ (self, mode, psi_1, fit_args):
    self.mode = mode
    self.psi_1 = psi_1
    self.fit_args = fit_args
  def quadFit(self, value_1, dvalue_1):
    ''' Fit a quadratic function from psi=0 to psi_i to the profile

        The first constraint is psi_1 is continuious. This first argument
        determines which other constraints are used. The options are
          0: profile at psi_0 and derivative at psi = 0
          1: profile at psi_0 and derivative at psi_1
          2: derivative at psi=0 and psi_1
        Second argument is the first constraint 
          profile at 0 for fit_args[0]= 0,1
          derivative at 0 for fit_args[0] = 2
        Third argument is the second constraint in needed
    '''
    if self.fit_args[0]==0:
      self.quadC = self.fit_args[1]
      self.quadB = self.fit_args[2]
      self.quadA = (value_1-self.quadB * self.psi_1 - self.quadC)/(self.psi_1**2)
    elif self.fit_args[0]==1:
      self.quadC = self.fit_args[1]
      self.quadA = (dvalue_1 * self.psi_1 +self.quadC - value_1)/self.psi_1**2
      self.quadB = dvalue_1 - 2.0 * self.quadA * self.psi_1
    elif self.fit_args[0]==2:
    # do something
      print("Todo")
    else:
      raise ValueError("Mode not defined for quad fit")
  def quadEval(self, psi):
    return self.quadA * psi**2 + self.quadB * psi + self.quadC
  def smooth (self, psi, field):
    ''' This function smooths a field field'''
    splineField = splrep(psi,field,k=3)
    field_psi_1 = splev(self.psi_1,splineField)
    dField_psi_1 = splev(self.psi_1,splineField,der=1)
    if self.mode=="quad":
      self.quadFit(field_psi_1,dField_psi_1)

    tempField = np.zeros(field.size)
    for ix, ipsi in enumerate(psi):
      if (ipsi < self.psi_1):
        if self.mode=="quad":
          tempField[ix] = self.quadEval(ipsi)
      else:
        tempField[ix] = field[ix]

    newSplineField = splrep(psi,tempField[:],k=3)

# plot fields
    x2 = np.linspace(0, 1, 200)
    y2 = splev(x2, splineField)
    y3 = splev(x2, newSplineField)
    plt.plot(x2, y2, x2, y3)
    plt.show()
    
    return tempField



home_dir = os.environ['HOME']
scratch_dir = "/SCRATCH"
local_dir = "/174446_novac_debug/eq21"
gfile = "g174446.03390"
new_gfile = gfile + ".mod3"
full_read_file = home_dir + scratch_dir + local_dir +"/"+ gfile
full_write_file = home_dir + scratch_dir + local_dir +"/"+ new_gfile

smooth_f = True
new_f_profile = modField("quad",  0.1, (1, -3.479))
smooth_p = True
new_p_profile = modField("quad",  0.1, (0, 86000., 0.0))

with open(full_read_file) as f:
  with open(full_write_file,'w') as out_file:
    mx, my = readWriteEfitHeader(f, out_file)
    print(mx,my)
    
    psi = np.zeros(mx)
    for ii in range(mx):
      psi[ii]= ii/(mx-1)
    f_profile = readProfile(f,mx)
    p_profile = readProfile(f,mx)
    q_profile = readProfile(f,mx)
    m_profile = readProfile(f,mx)
    if smooth_f:
      f_profile = new_f_profile.smooth(psi,f_profile)
    if smooth_p:
      p_profile = new_p_profile.smooth(psi,p_profile)


    writeProfile(out_file,f_profile)
    writeProfile(out_file,p_profile)
    writeProfile(out_file,q_profile)
    writeProfile(out_file,m_profile)
    readWriteRestFile(f,out_file)

#  with(full_write_file,'w') as out_file:
# read the file header


