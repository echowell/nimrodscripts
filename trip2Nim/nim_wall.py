#!/usr/bin/env python3

import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import scipy.optimize as opt
class NimrodWall:
  ''' Base class for dealing with NIMROD wall files '''
  wall_file = ''
  wall_modes = -1 # Number of Fouier modes in wall
  plt_wall = False
  check_ft = False # check my Fourier transform
  wall_points = 500 # points to use in high res wall
  file_rz = np.zeros([1,2])
  rm_wall = np.zeros(1)
  zm_wall = np.zeros(1)
  rz_wall = np.zeros([1,2])
  rz2_wall = np.zeros([1,2])
  def __init__(self, file, modes, plot_wall=False, check_ft=False,*args, **kwargs):
    ''' Initalize NimrodWall class '''
    self.wall_file=file
    self.wall_modes=modes
    self.plt_wall=plot_wall
    self.check_ft=check_ft
    self.read_file()
    self.get_fourier()
    self.get_wall_rz()
    if self.check_ft:
      self.rz2_wall =np.zeros([self.wall_points,2])
      for itheta, theta in enumerate(np.linspace(0,2*np.pi,num=self.wall_points)):
        self.rz2_wall[itheta,0],self.rz2_wall[itheta,1] = self.get_rz_theta(theta)
    if self.plt_wall:
      self.plot_rz()
  def read_file(self):
    ''' Read a nimrod wall file '''
    self.file_rz = np.loadtxt(self.wall_file, delimiter=',', skiprows=1)
  def get_fourier(self):
    ''' Calculate the Fourier coffiecents of R and Z 

        The fourier modes are stored as followes:
          index 0 : m=0
          index 1 : real m=1
          index 2 : imag m=1
          index 2l-1 : real m = l
          index 2l   : imag m = l
    '''
    rm = fft.rfft(self.file_rz[:,0])/self.file_rz.shape[0]
    zm = fft.rfft(self.file_rz[:,1])/self.file_rz.shape[0]
    if self.wall_modes < 0: # use all the Fourier modes
      self.rm_wall = rm
      self.zm_wall = zm
    else: # Truncate Fourier spectrum 
      max_index = min(2*self.wall_modes+1,rm.size)
      self.rm_wall = rm[:max_index]
      self.zm_wall = zm[:max_index]
  def get_wall_rz(self):
    ''' Calculate a high resolution wall from the fourier RZ '''
    r_wall = fft.irfft(self.rm_wall,self.wall_points)*self.wall_points
    z_wall = fft.irfft(self.zm_wall,self.wall_points)*self.wall_points
    self.rz_wall = np.stack((r_wall,z_wall),axis=-1) 
  def get_rz_theta(self, theta):
    ''' Return a given RZ location on the wall as a function of theta 

        Here theta is an angle in radians
    '''
    rt = self.rm_wall[0]
    zt = self.zm_wall[0]
    m = self.rm_wall.size
    if self.rm_wall.size % 2: #odd
      m_max = int((self.rm_wall.size-1)/2)
      m_modes = range(1,m_max+1)
    else: # even
      m_max = int(self.rm_wall.size)/2+1
      m_modes = range(1,m_max)
      rt += self.rm_wall[-1]*np.cos(theta * m_max)
      zt += self.zm_wall[-1]*np.cos(theta * m_max)
    for im, m in enumerate(m_modes,1):
      rt += 2.0*(self.rm_wall[2*im-1]*np.cos(theta*m)-self.rm_wall[2*im]*np.sin(theta*m))
      zt += 2.0*(self.zm_wall[2*im-1]*np.cos(theta*m)-self.zm_wall[2*im]*np.sin(theta*m))
    return rt, zt
  def get_drz_dtheta(self, theta):
    ''' Return a the theta derivative of RZ along the wall at theta 

        Here theta is an angle in radians
    '''
    rt = 0
    zt = 0
    m = self.rm_wall.size
    if self.rm_wall.size % 2: #odd
      m_max = int((self.rm_wall.size-1)/2)
      m_modes = range(1,m_max+1)
    else: # even
      m_max = int(self.rm_wall.size)/2+1
      m_modes = range(1,m_max)
      rt += self.rm_wall[-1]*np.sin(theta * m_max)*m_max
      zt += self.zm_wall[-1]*np.sin(theta * m_max)*m_max
    for im, m in enumerate(m_modes,1):
      rt += 2.0*m*(self.rm_wall[2*im-1] * np.sin(theta*m) + self.rm_wall[2*im]*np.cos(theta*m))
      zt += 2.0*m*(self.zm_wall[2*im-1] * np.sin(theta*m) + self.zm_wall[2*im]*np.cos(theta*m))
    return rt, zt
  def get_theta_rz(self,r,z):
    ''' Return the theta value at r,z on the surface '''
    #opts = {"maxiter" : 10000}
    opts = {}
    dist = (self.rz_wall[:,0]-r)**2 + (self.rz_wall[:,1]-z)**2
    guess = np.argmin(dist)/self.rz_wall.shape[0] * 2*np.pi
    opt_theta = opt.minimize(self.dist,guess,(r,z),method='CG',tol=1e-5,options=opts)
    if not opt_theta.success:
      print("get_theta_rz optimze failed with message: " + opt_theta.message)
    return opt_theta.x[0],self.dist(opt_theta.x[0],r,z)

  def dist(self,theta, *args):
    ''' Return the distance between a point on the surface and r,z '''
    ro = args[0]
    zo = args[1]
    r,z = self.get_rz_theta(theta)
    dist = np.sqrt((ro-r)**2 + (zo-z)**2)
    return dist
  def plot_rz(self):
    plt.plot(self.file_rz[:,0],self.file_rz[:,1], 'bo', label="Orginal Wall")
    plt.plot(self.rz_wall[:,0],self.rz_wall[:,1], 'g-', label="Smooth Wall")
    if self.check_ft:
      plt.plot(self.rz2_wall[:,0],self.rz_wall[:,1], 'k-', label="2nd Wall")
    plt.legend(loc=1)
    plt.show()