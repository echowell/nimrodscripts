#!/usr/bin/env python3

#
# Input files:
# Ouput file:


import numpy as np
import sys
import matplotlib.pyplot as plt


def bump_function(amp,time,time_scale,time_offset):
  t_norm = (time-time_offset)/time_scale
  t_norm_sq = t_norm * t_norm
  if t_norm_sq >= 1.0:
    return 0.0
  else:
    return amp*np.exp(- 1.0/(1.0-t_norm_sq))


amp = 1.0
t_scale = 0.50
t_offset = 0.505
npts = 1000

time = np.linspace(0,10*t_offset,npts)
bump = np.zeros([npts])
for it, tt in enumerate(time):
  bump[it] = bump_function(amp,tt,t_scale,t_offset)

fig =plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.set_title("Magnetic Perturbation Pulse Shape",fontsize=18)
ax.plot(time,bump)
ax.set_xlabel("Time [ms]",fontsize=18)
ax.set_ylabel("Amplitude [A.U.]",fontsize=18)
plt.tight_layout()
plt.show()
