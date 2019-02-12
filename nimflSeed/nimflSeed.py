#!/usr/local/bin/python3
''' This script randomly generates a collection of seed locations of NIMFL
field line integration'''

################################################################################
#  Set up envirment and import modules
################################################################################
import sys
sys.path.insert(0, "./")
import os
import startPosClass as sp
pwd = os.getcwd()
homeDir = os.environ['HOME']
################################################################################
#  User defined input
################################################################################
fileName = "start_positions.dat"
writeDirectory = pwd
writeDirectory = homeDir + "/SCRATCH/test_vac_fl"
nPoints = 100
geom = 'd3dlower'
phiZero = 0.0
randomPhi = True # Ture
################################################################################
#  Set up auxiliary variables 
################################################################################
writeFile = writeDirectory + "/" + fileName
################################################################################
#  Run code
################################################################################
#initalize start position object
startPos = sp.startPosClass(nPoints,geom,randomPhi,phiZero)
# generate starting points
startPos.calculateRZPhi()
# write output
startPos.writeStartPos(writeFile)