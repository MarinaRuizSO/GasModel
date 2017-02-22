""" 
Lennard Jones: 
returns the force and potential energy acting on the ith particle through the Lennard-Jones pair potential 

Marina Ruiz Sanchez-Oro

Cara Lynch

21/02/2017
"""

import math
import numpy as np
from ParticleList import ParticleSyst as P


def ljforce(system,boxdim,Rc) :
    """
    Computes the Lennard-Jones force acting on each particle in the system as an (N,3) Numpy array.
	
    :param system: ParticleSyst object representing the system of N particles
    :param boxdim: Box dimensions as an (1,3) Numpy array
    :param Rc: cutoff radius as a float
    :return force: an (N,3) Numpy array where the ith row is the force experienced by the ith particle
    """
    force = np.empty(shape={P.N,3})
	
    # Compute for i th particle
    for i in range(0, P.N):
        # With relation to j th particle ---would this count the particles twice?---
        for j in range(0, P.N):
            if j == i:
                # not computing the effect the particle has on itself
                pass
            else:
                vecsep = np.empty(shape={1,3})
                for k in range(0,3):
                    
                    # Calculate the vector separation of the particles adhering to the minimum image convention
                    vecsep[k] = P.position[i,k] - P.position[j,k]
                    
                    # Minimum image convention
                    while abs(vecsep[k]) > abs(boxdim[k]/2) :
                        
                        # If the separation is negative, then the image will be -L/2 in the direction considered
                        if vecsep[k] < 0:
                            vecsep[k] = vecsep[k] - abs(boxdim[k])/2
                            
                        # If the separation is positive, then the image will be +L/2 in the direction considered
                        else:
                            vecsep[k] = vecsep[k] + abs(boxdim[k])/2
            # Magnitude of vector separation
            r = math.sqrt(np.inner(vecsep, vecsep))
			
		# Cutoff radius condition
            if r > Rc :
                # If above cutoff radius set the force to 0
                for k in range(0,3):
                    force[i,k] = force[i,k] + 0
		# If not, add this contribution to the force
            else:
               for k in range (0,3):
                   force[i,k] = force[i,k] + (48*((1/r**14) - (1/ (2 * r**8)) ) ) * vecsep[k]
				      
    return force
	
def ljpotential(system,boxdim,Rc) :
    """
    Computes the Lennard-Jones potential of each particle in the system as an (N,1) Numpy array.
    
    :param system: ParticleSyst object representing the system of N particles
    :param boxdim: Box dimensions as a (1,3) Numpy array
    :param Rc: cutoff radius
    :return force: an (N,1) Numpy array where the ith row is the potential of the ith particle
    """
    
    potential = np.empty(shape={P.N,1})
    
    # Compute for i th particle
    for i in range(0, P.N):
        
        # With relation to j th particle
        for j in range(0, P.N):
            
            if j == i:
                # not computing the effect the particle has on itself
                pass
            
            else:
                vecsep = np.empty(shape={1,3})
                for k in range(0,3):
                    
                    # Calculate the vector separation of the particles
                    vecsep[k] = P.position[i,k] - P.position[j,k]
                    
                    # Minimum image convention
                    while abs(vecsep[k]) > abs(boxdim[k]/2) :
                        # If the separation is negative, then the image will be -L/2 in the direction considered
                        if vecsep[k] < 0:
                            vecsep[k] = vecsep[k] - abs(boxdim[k])/2
                        # If the separation is positive, then the image will be +L/2 in the direction considered
                        else:
                            vecsep[k] = vecsep[k] + abs(boxdim[k])/2
                # Magnitude of vector separation
                r = math.sqrt(np.inner(vecsep, vecsep))
                
                # Cutoff radius
                if r > Rc :
                    potential[i] = potential[i] + 0
			# Add this contribution to the overall potential array
                else:
                    potential[i] = potential[i] + (4*((1/r**12) - (1/ (r**6)) ) )
				      
    return potential

def totPE(system, boxdim, Rc):
    """
    Computes the total potential energy of the system according to the Lennard-Jones potential.
    
    :param system: ParticleSyst object representing the system of N particles
    :param boxdim: Box dimensions
    :param Rc: cutoff radius
    :return totalPE: a float representing the potential of the system
    """
    # Compute the potential of each individual particle
    potential = ljpotential(system,boxdim,Rc)
    
    # Set initial potential as a float
    doubletotalPE = 0.0
    
    # Add all potentials
    for i in range(0,P.N):
        doubletotalPE += potential[i]
    # Every interaction is counted twice so must divide potential calculated by 2
    return doubletotalPE/2

 
