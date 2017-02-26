"""
A class representing a list of N particles moving in 3D space, complete with various methods.

Authors: Cara Lynch

Date: 07/02/2017

"""

import numpy as np
import math
from LennardJones import LennardJones as lj


class ParticleSyst(object) :

    def __init__(self, name, label, pos, vel, mass, N):
        """
        Initialise a ParticleSyst instance
        :param name: name of system as a string
        :param label: labels of particles as an (N,1) Numpy array
        :param N: number of particles in system
        :param pos: positions of the particles as an (N,3) Numpy array
        :param vel: velocities of the particle as an (N,3) Numpy array
        :param mass: mass of the particles as an (N,1) Numpy array
        """
        self.name = name
        self.label = label
        self.position = pos
        self.velocity = vel
        self.N = N
        self.mass = mass


    def __str__(self):
        """
        Print the label of the system and relevant information

        :param ParticleSyst: ParticleSyst instance
        :return: string with label, number of particles
        """     
        return  str(self.name) + " represents a system of " + str(self.N) + " Particle3D instances"
 

    def velmag(self, i):
        """
        Calculates the magnitude of the velocity of the particle of index i
	
        :param ParticleSyst: ParticleSyst instance
        :return: speed of particle as float
        """
        speedsquared = 0.0
        for k in range(0,3):
            speedsquared = speedsquared + (self.velocity[i,k])**2
            speed = math.sqrt(speedsquared)
        return speed
	
    # Kinetic energy
    def kineticEnergy(self):
        """
        Return the kinetic energy of the system
	
        :param ParticleSyst: ParticleSyst instance
        :return: kinetic energy of system as float
        """
        energy = 0.0
        for i in range(0, self.N):
            energy = energy + 0.5*self.mass[i]*self.velmag(i)**2
        return energy

# Method to compute total energy of the system
    def energy(self, boxdim, r_c):
        return self.kineticEnergy() + lj.totPE(self, boxdim, r_c)

    # Time integration methods
    # First-order velocity update
    def leapVelocity(self, dt, force):
        """
        Update the velocity to the first order

        :param ParticleSyst: ParticleSyst instance
        :param dt: timestep as float
        :param force: force between every particle represented by an (N,3) Numpy array
        :return: updated velocity of particle represented by an (N,3) Numpy array
        """
        self.velocity = self.velocity + dt*np.divide(force, self.mass)

    # First-order position update
    def leapPos1st(self, dt):
        """
        Update the position to the first order
	
        :param ParticleSyst: ParticleSyst instance
        :param dt: timestep as float
        :return: first order update of position of particle represented by a (N,3) Numpy array
        """
        self.position = self.position + dt*self.velocity

    # Second-order position update
    def leapPos2nd(self,dt,force):
        """
        Update the position to the second order
	
        :param ParticleSyst: ParticleSyst instance
        :param dt: timestep as float
        :param force: force as a vector represented by an (N,3) Numpy array
        :return: second order update of position of particle represented by an (N,3) Numpy array
        """
        self.position = self.position + dt*self.velocity + 0.5*dt**2*np.divide(force, self.mass)
	
    # Prints in format necessary for VMD
    def printVMD(self, k):
        """
        Prints strings in the right format for a VMD trajectory file

        :param ParticleSyst: ParticleSyst instance
        :param k: timestep number as float
        :return: N strings
        """
        total = str()
        total = total + str(self.N) + "/n" + "Point = " + str(k) + "/n"
        for i in range(0,self.N):
            total = total + str(self.label[i]) + " " + str(self.position[i,0]) + " " + str(self.position[i,1]) + " " + str(self.position[i,2]) + "/n"
        return total
		 

    # Create a particle from a file entry
    @staticmethod
    def createsystem(fileIn):
        """
        Create a system of N particles from a file entry.
        The first line of the file contains the system name and the number of particles, N, as a float.
        Each subsequent line of the file must have the label, mass, position and velocity coordinates of the relevant particle in that order, separated by spaces.
	
        :param fileIn: file opened for reading
        :return: ParticleSyst as an instance
        """
        alllines = fileIn.readlines()
        line0 = alllines[0].split()
        name = str(line0[0])
        N = float(line0[1])
	
        # Create Numpy arrays to hold system information
        label = np.empty([N,1], dtype = str)
        mass = np.empty([N,1], dtype = float)
        pos = np.empty([N,3], dtype = float)
        vel = np.empty([N,3], dtype = float)
	
        for i in range(1,N+1):
		
		# Define list of elements of line i corresponding to particle of index i
              line = alllines[i].split()
		
		# Add label of particle i to label array
              label[0] = str(line[0])
		
		# Add mass of particle i to mass array
              mass[0] = float(line[1])
		
		# Add position of particle i to row i of position array
              pos[i,0] = float(line[2])
              pos[i,1] = float(line[3])
              pos[i,2] = float(line[4])
		
		# Add velocity of particle i to 
              vel[i,0] = float(line[5])
              vel[i,1] = float(line[6])
              vel[i,2] = float(line[7])

        return ParticleSyst(name, label, pos, vel, mass, N)
