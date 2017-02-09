"""
Test program for Lennard Jones force and potential calculations

"""

from Particle3D import Particle3D as P
import numpy as np
import math
import sys
from MIC import vecsep
import LennardJones as LJ


#Read name of input file
if len(sys.argv)!=5:
    print "Wrong number of arguments."
    quit()
else:
    infile1 = sys.argv[1] #particle1
    infile2 = sys.argv[2] #particle2
    L = float(sys.argv[3]) #box length
"""
    r_c = float(sys.argv[4]) #cutoff radius
"""
#Create particles
p1_file = open(infile1, "r")
p1 = P.createparticle(p1_file)

p2_file = open(infile2, "r")
p2 = P.createparticle(p2_file)

print LJ.ljforce(p1,p2,L,r_c)
print LJ.ljpotential(p1,p2,L,r_c)


 
