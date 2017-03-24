"""
 CMod Ex3: Particle3D, a class to describe 3D particles
"""

import numpy as np
import math


class Particle3D(object):

   

    # Initialise a Particle3D instance
    def __init__(self, pos, vel, mass):
        self.position = pos
	self.velocity = vel
	self.mass = mass
    
    
    # Formatted output as String
    def __str__(self):
        return "x = " + str(self.position[0]) + ", y = " + str(self.position[1]) + "z = " + str(self.position[2]) 
    

    # Kinetic energy
    def kineticEnergy(self):
        sq= sum(self.velocity*self.velocity)
        return 0.5*self.mass*sq

    # Time integration methods
    # First-order velocity update
    def leapVelocity(self, dt, force):
        dvel= (dt/self.mass)*force
        self.velocity = np.add(self.velocity, dvel)

    # First-order position update
    def leapPos1st(self, dt):
        dpos = dt*self.velocity
        self.position = np.add(self.position, dpos)

    # Second-order position update
    def leapPos2nd(self, dt, force):
	dpos1= dt*self.velocity
        dvel1= (0.5*dt**2/self.mass)*force
        change = np.add(dpos1, dvel1)
        self.position = np.add(self.position, change)

    def from_file(file_handle):
        line = inFile.readline()
        tokens = line.split()
        posit = np.array[float(tokens[0]), float(tokens[1]), float(tokens[2])]
        line = inFile.readline()
        tokens = line.split()
        veloc = np.array[float(tokens[0]), float(tokens[1]), float(tokens[2])]
        line = infile.readline()
        tokens = line.split()
        mass = float(tokens[0])
        return Particle3D(posit, veloc, mass)

    def separation(x, y):
        return np.subtract(x.position, y.position)
       
    def LJpotential(r):
        magr = math.sqrt(sum(r*r))
        return 4*((1/(magr**12))-(1/(magr**6)))      

    def LJforce(r):
        magr = math.sqrt(sum(r*r))
        return 48*((1/(magr**14))-(1/(2*(magr**8))))*(r)

