# -*- coding: utf-8 -*-
"""
Radial Distribution Function algorithm

calculates the probability to find a particle at a given distance from the reference particle using a histogram

Created on Sun Feb 26 11:43:43 2017

@author: Cara Lynch
"""

from ParticleList import ParticleSyst as P
import matplotlib.pyplot as plt
import numpy as np



def particledistances(syst, fileName, k):
    """
    Collects distances between pairs of particles at a specific point in time with respect to the k th particle and exports into a file.
    
    Has to be integrated into VV code
    :param syst: N body system represented as a ParticleSyst instance
    :param fileName: name of file which the trajectories will be exported to as a string
    :param k: index of reference particle
    :param dt: time interval
    """
    
    # Create empty array to hold radial distances from reference particle at each timestep (should this be 0,0?)what 
    radialdistance = np.empty([P.N(syst) - 1,1])
    """
    # Open output file for writing
    outfile = open(fileName, "w")
    outfile.write("{0:f} \n".format(radialdistance[0]))
    """
    for i in range(0, P.N(syst)):
        
        # If reference particle, do not calculate distance
        if i == k:
            pass
        else:
            
            #Calculate magnitude of vector separation, assuming PBCs have already been applied
            squaredDistance = 0.0
            for u in range(0,3):
                squaredDistance += (syst.position[i,u] - syst.position[k,u])**2
            radialdistance[i] = squaredDistance
        
        outfile.write("{0:f} \n".format(radialdistance[i]))
        
def histogram(fileName):
    
    # Open the file of radial distances for reading
    fileIn = open(fileName, "r")
    distances = fileIn.readlines()
    
    # Count number of elements in distances list to create array
    num = len(distances)
    rdfArray = np.empty([1,num])
    
    # Assign list values as floats into array
    for i in range(0,num +1):
        rdfArray[i] = float(distances[i])
    
    # Create and plot normalised histogram
    plt.histogram(rdfArray, bins = 'auto', density = True)
    plt.title("Radial distribution function")
    plt.show()
