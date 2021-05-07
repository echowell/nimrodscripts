#!/usr/bin/env python3
import numpy as np
import itertools as iter
import argparse
import pathlib
import subprocess
from shutil import copy2
import os

SEGMENTS = [
    ((1.016, 0.000), (1.016,-1.223)),
    ((1.016,-1.223), (1.153,-1.363)),
    ((1.153,-1.363), (1.372,-1.363)),
    ((1.372,-1.250), (1.682,-1.250))
]

def write_start(file_name, rz, phi):
    npts = rz.shape[0]
    with open(file_name,'w') as thisFile:
        thisFile.write(str(npts)+"\n")
        for idx, in np.ndindex(rz.shape[0]):
            thisLine = '{: 16.16e}'.format(rz[idx,0]) + ", "
            thisLine+= '{: 16.16e}'.format(rz[idx,1]) + ", "
            thisLine+= '{: 16.16e}'.format(phi) + "\n"
            thisFile.write(thisLine)

def div_seed_runner(nimfl, dump, **kwargs):
    nphi = 360
    nl = 100
    dtan = 0.01
    phi_array = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    BASENAME = "start_positions.dat"
    DIRNAME = "temp_dir"
    counter = 0
    for segment in SEGMENTS:
        start = np.array(segment[0])
        end = np.array(segment[1])
        tan = np.array([-end[1] + start[1], end[0] - start[0]])
        length = np.linalg.norm(tan)
        tan = dtan/length * tan
        l_arr = np.linspace(start,end,nl)
        l_arr = l_arr + tan
        for phi in np.nditer(phi_array):
            dir = DIRNAME + f"_{counter}"
            print(f"Makeing dir {dir}")
            try:
                os.mkdir(dir)
            except:
                print(f"{dir} exists")
            copy2(dump,dir)
            copy2('nimrod.in',dir)
            copy2('nimfl.in',dir)
            copy2('dmap_rz.dat',dir)
            #enter dir
            os.chdir(dir)
            file_name = BASENAME
            write_start(BASENAME, l_arr, phi)
            subprocess.run(nimfl)
            os.remove('nimrod.in')
            os.remove('nimfl.in')
            os.remove('dmap_rz.dat')
            os.remove(dump.name)
            os.chdir('../')
            counter += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NIMFL divertor plate runner')
    parser.add_argument('nimfl_exe',
                        action='store',
                        type=pathlib.Path,
                        help="nimfl exectable path" )
    parser.add_argument('dumpfile',
                        action='store',
                        type=pathlib.Path,
                        help="dumpfile path" )

    args = vars(parser.parse_args())

    div_seed_runner(args['nimfl_exe'], args['dumpfile'], **args)
