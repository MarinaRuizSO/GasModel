""" 
Lennard Jones: 
returns the force and potential energy between two particles interacting through the Lennard-Jones pair potential 

Marina Ruiz Sanchez-Oro

07/02/2017
"""
import sys
import math
import matplotlib.pyplot as pyplot
import numpy as np
from Particle3D import Particle3D as P
from MIC import vecsep 


# Calculates the force between two particles based on the separation, the cutoff radius and the box dimensions 

def ljforce (p1,p2,L,r_c):
    # Computes the vector separation of the particles
    rvector = mic.vecsep(p1,p2,L)
    # Computes the magnitue of the particle separation 
    rscalar = math.sqrt(np.inner(rvector,rvector))
    if rscalar > r_c:
        return 0
    else: 
        return (48*((1/rscalar**14)-(1/(2*rscalar**8))))*rvector

# Calculates the Lennard Jones potential for the particle pair. 			
def ljpotential(p1,p2,L,r_c):
    # Computes the vector separation of the particles
    rvector = mic.vecsep(p1,p2,L)
    # Computes the magnitude of the particle separation 
    rscalar = math.sqrt(np.inner(rvector,rvector)
    if rscalar > r_c:
        return 0
    else: 
	return (4*((1/rscalar**12)-(1/(rscalar**6))))*rvector



	

 
