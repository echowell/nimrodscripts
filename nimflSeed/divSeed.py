#!/usr/bin/env python3
import numpy as np
import itertools as iter
import random as ran

SEGMENTS = [
    ((1.016, 0.000), (1.016,-1.223)),
    ((1.016,-1.223), (1.153,-1.363)),
    ((1.153,-1.363), (1.372,-1.363)),
    ((1.372,-1.250), (1.682,-1.250))
]

def write_start(file_name, xy, phi):
    npts = xy.shape[0] * phi.shape[0]
    with open(file_name,'w') as thisFile:
        thisFile.write(str(npts)+"\n")
        for [r,z],p in iter.product(xy,phi):
            thisLine = '{: 16.16e}'.format(r) + ", "
            thisLine+= '{: 16.16e}'.format(z) + ", "
            thisLine+= '{: 16.16e}'.format(p) + "\n"
            thisFile.write(thisLine)

def div_seed():
    nphi = 360
    nl = 100
    dtan = 0.01
    BASENAME = "start_positions.dat"
    counter = 0
    for segment in SEGMENTS:
        start = np.array(segment[0])
        end = np.array(segment[1])
        tan = np.array([-end[1] + start[1], end[0] - start[0]])
        length = np.linalg.norm(tan)
        tan = dtan/length * tan
        l_arr = np.linspace(start,end,nl)
        l_arr = l_arr + tan
        phi = np.linspace(0, 2*np.pi, nphi)
        file_name = BASENAME + f"_{counter}"
        write_start(file_name, l_arr, phi)
        counter += 1
if __name__ == "__main__":
    div_seed()
