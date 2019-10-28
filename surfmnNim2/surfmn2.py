#!/usr/local/bin/python3

import sys
import xySliceClass as xy
import surfmnClass as mn

def main(argv):
  thisSurfmn=mn.SurfmnClass(argv)
  fields=xy.xyClass(thisSurfmn.xyFile,thisSurfmn.mx,thisSurfmn.my,thisSurfmn.pd,False)
  thisSurfmn.calcBmn(fields)
  thisSurfmn.plotBmn()
  print(thisSurfmn.xyFile)
  print(fields.ix.shape)

main(sys.argv[1:])