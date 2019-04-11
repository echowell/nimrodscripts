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
writeDirectory = homeDir + "/SCRATCH//166439/footpoint_03300_q104/lphi5/S7Pr1e2_nimfl2"
nPoints = 3000
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