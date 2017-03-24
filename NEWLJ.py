import sys
import matplotlib.pyplot as pyplot
from Particle3Dmodd import Particle3D
import numpy as np
import math
import MDUtilities as MD 


def vecsepmin(particle1, particle2, boxSize):
    xsep= particle1.position[0]-particle2.position[0] #separation between two particles on the x y and z axis.
    ysep= particle1.position[1]-particle2.position[1]
    zsep= particle1.position[2]-particle2.position[2]
    if xsep>(boxSize/2):
        xsep=xsep-boxSize
    if xsep<((-1)*(boxSize/2)):
        xsep=xsep+boxSize
    if ysep>(boxSize/2):
        ysep=ysep-boxSize
    if ysep<((-1)*(boxSize/2)):
        ysep=ysep+boxSize
    if zsep>(boxSize/2):
        zsep=zsep-boxSize
    if zsep<((-1)*(boxSize/2)):
        zsep=zsep+boxSize
    r= np.array([xsep, ysep, zsep])#separation between particle and an image particle 
    return r

def forceTwoParticles(part1, part2, r, rcut):
	magr = math.sqrt(sum(r*r))
	magrcutoff = math.sqrt(3*(rcutoff**2))
	# Set up Lennard Jones force
	if magr < magrcutoff: 
		force = 48*((1/(magr**14))-(1/(2*(magr**8))))*(r)
	else:
		force = 0*r
	return force

def computeForces(particlelist):
	force_list = []
	for l in range(len(particlelist)):  #the number of particles in the list
		force =0
        	for m in range(len(particlelist)): #for each particle acting on particle l
			if l != m:
				force += forceTwoParticles(particlelist[l], particlelist[m], rcut)
			force_list.append(force)
	return force_list

def writetofile(particlelist, outfile, time):
	# Write initial condition to output file
	n = len(particlelist)
	outfile.write(str(n) + "\n" + "Point = " + str(time) + "\n")
	for i in range(len(particlelist)):
    		part1=i+1
    		outfile.write(str(part1) + " " + str(particlelist[i].position[0]) + " " + str(particlelist[i].position[1]) + " " + str(particlelist[i].position[2]) + "\n")


def updateVelocities(dt, particle, force_before, force_new):
	force=0.5*(force_before+force_new)
	particle.leapVelocity(dt, force)
	



def UpdatePosition(particle, dt, Force, boxSize):
	particle.leapPos2nd(dt, Force)
	if particlelist[i].position[0]>boxSize:
    		xpos = particlelist[i].position[0]%boxSize
       		particlelist[i].position[0]=xpos
	if particlelist[i].position[0]<=0.0:
        	xps=particlelist[i].position[0]%boxSize
        	xpos=boxSize-xps
        	particlelist[i].position[0] = xpos #Applys periodic boundary condition, assuming however far out the box, it will be replaced with a particle inside the box and appaends that new postion to particle i for x,y,z.
	if particlelist[i].position[1]>boxSize:
        	ypos=particlelist[i].position[1]%boxSize
        	particlelist[i].position[1]=ypos
	if particlelist[i].position[1]<=0.0:
        	yps=particlelist[i].position[1]%boxSize
        	ypos=boxSize-yps
        	particlelist[i].position[1]=ypos
	if particlelist[i].position[2]>=boxSize:
        	zpos=particlelist[i].position[2]%boxSize
        	particlelist[i].position[2]=zpos
	if particlelist[i].position[2]<=0.0:
        	zps=particlelist[i].position[2]%boxSize
        	zpos= boxSize-zps
        	particlelist[i].position[2]=zpos

#Set up energies
def PotentialCalc (r, rcutoff) :
	magr= math.sqrt(sum(r*r))
	magrcutoff= math.sqrt(3*(rcutoff**2))
	if magr < magrcutoff :
		pot = 4*((1/(magr**12))-(1/(magr**6)))
		return pot
	else:
		return 0.0
	

#Calculate the Total system kinetic energy
def TotalKE(particlelist):
	kinetic = 0
	for i in range(len(particlelist)):
		newKE= particlelist[i].kineticEnergy()
		kinetic = kinetic+newKE
	return kinetic

#Update MSD values for all particles
def MSDUpdate(particlenumber, initpos, position, step):
	N= step+1
	prevN = step
	if prevN == 0:
		prevN=1
	disp= abs(position-initpos)
	msdlist[particlenumber]=((prevN/float(N))*msdlist[i])+((1./N)*(disp)**2)

#Calculate the average Mean Squared Distribution over all particles
def AvgMSD(particlelist, msdlist):
	avmsd = 0.0
	for d in range(len(particlelist)):
	    magmsd=math.sqrt(sum(msdlist[d]**2))
            avmsd= avmsd+((1.0/number_particles)*(magmsd))
	avgmsdlist.append(avmsd)

#put RDF here
def RDF(separation):
    countlist= []
    radiilist= []
    maxr= math.sqrt(3*(boxSize**2)) #ditance from center to corner to give max distance
    bins= int(len(separation)**0.5) #Determines number of bins for histogram
    flbins= float(bins) #converts bins into float number for width calculation
    width= maxr/bins
    dr= 0.0
    p=0
    while dr<maxr:
        countlist.append(0)
        leftlimit=dr
        rightlimit=leftlimit+width
        for i in range(len(separation)):
            sep=separation[i]
            magsep= math.sqrt(sum(sep*sep))
            if magsep>leftlimit:
                if magsep<rightlimit:
                    countlist[p]+= 1
        rval= 0.5*(leftlimit+rightlimit)
        normfactor= 4.0*math.pi*(rval**2)*rho*width
        countlist[p]=countlist[p]/normfactor
        radiilist.append(rval)
        dr= dr+width
        p += 1
    pyplot.figure(1)
    pyplot.plot(radiilist, countlist)
    pyplot.xlabel('Radius')
    pyplot.ylabel('No. of particles')
    pyplot.title('Radial distribution function')
    pyplot.show

	
# Read name of output file from command line
if len(sys.argv)!=3:
    print "Wrong number of arguments."
    #print "Usage: " + sys.argv[0] + " <input file> <output file>"
    quit()
else:
    outfileName = sys.argv[2]
    inputfileName = sys.argv[1]

# Open output file for writing
inputfile = open(inputfileName, 'r')
outfile = open(outfileName, "w")


#Set up parameters (read input file and set values accordingly)
line = inputfile.readline()
tokens = line.split()
number_particles = int(tokens[0])
T= float(tokens[1])
rho= float(tokens[2])
# Set up simulation parameters
numstep = int(tokens[3])
time = 0.0
dt = float(tokens[4])
resolution = float(tokens[5])
rcutoff= float(tokens[6])

# Initialize the particle list
particlelist =[]
for i in range(number_particles):
	particlelist.append(Particle3D(np.array([0.0,0.0,0.0]),np.array([0.0,0.0,0.0]),1.0))

#uses MDUtilities to initialise each particle in the list    
boxsize = MD.setInitialPositions(rho, particlelist)
MD.setInitialVelocities(T, particlelist)

boxSize = boxsize[0]


#Create lists to store energy, time and average MSD
kineticlist= []
potentiallist= []
totallist= []
avgmsdlist= []

#Create two time lists, one to store time values for energy graph, one for time values when Average MSD is being calculated.
timelist=[]
msdtimelist= []

#Create a list of empty values to store forces on particles and MSD.
forcelist=[]
msdlist= []
for i in range(len(particlelist)):
	forcelist.append(0)
	msdlist.append(0)


# write initial conditions
writetofile(particlelist, outfile, 0)

#Store initial positions in a list
initposlist=[]
for i in range(len(particlelist)):
	initposlist.append(particlelist[i].position)

# Start the time integration loop

for k in range(numstep):
	#Update time value and append it to timelist
	time += dt
	timelist.append(time)
	#Work out systems kinetic energy and store to list
	kinetic= TotalKE(particlelist)
        kineticlist.append(kinetic)
	#create bin for system's potential energy
	potential=0.0    
	# Update particle positions
	for i in range(len(particlelist)):
		forcelist[i]= 0
		for j in range(len(particlelist)):
			if i != j:
				r = vecsepmin(particlelist[i], particlelist[j], boxSize)
	#calculate potential on particle and add to system bin
				pot= PotentialCalc(r, rcutoff)
				potential += pot
	#Calculate force and update positions
				force = forceTwoParticles(particlelist[i], particlelist[j], r, rcutoff)
				forcelist[i] += force
	#Update position and apply PBC
		UpdatePosition(particlelist[i], dt, forcelist[i], boxSize)
	#Append system potential to potential list
	potentiallist.append(potential)
	#Work out system total energy and append to list
	Total= kineticlist[k]+potentiallist[k]
	totallist.append(Total)	
	#After all particle positions have been updated, update velocity using force values from last and current positions
	for i in range(len(particlelist)):
		force= 0		
		for j in range(len(particlelist)):
			if i != j:
				r = vecsepmin(particlelist[i], particlelist[j], boxSize)
				force += forceTwoParticles(particlelist[i], particlelist[j], r, rcutoff)
		updateVelocities(dt, particlelist[i], forcelist[i], force)				
	if (k+1)%resolution == 0:
		writetofile(particlelist, outfile, (k+1))

	#Update MSD
	for i in range(len(particlelist)):	
		MSDUpdate(i, initposlist[i], particlelist[i].position, k)
	if k%resolution == 0:
		msdtimelist.append(time)
		AvgMSD(particlelist, msdlist)	


#Calculate RDF and produce graph
	if k%(numstep-1) == 0:
		separationslist=[]
		for i in range(len(particlelist)):
			for j in range(len(particlelist)):
				if i<j:
					r = vecsepmin(particlelist[i], particlelist[j], boxSize)
					separationslist.append(r)
		RDF(separationslist)
	

#Plot graph of Average MSD against time


pyplot.figure(2)
pyplot.plot(msdtimelist, avgmsdlist)
pyplot.ylabel('Average MSD')
pyplot.xlabel('Time')
pyplot.title('Average Mean Squared Displacement against Time')
pyplot.show

#Plot graph of energies against time


pyplot.figure(3)
pyplot.plot(timelist, totallist, 'r')
pyplot.plot(timelist, kineticlist, 'b')
pyplot.plot(timelist, potentiallist, 'g')
pyplot.xlabel('Time')
pyplot.ylabel('Energy')
pyplot.title('System energy against time')
pyplot.show()


# Close output file
outfile.close()


# Plot graph of x and y position amplitude vs time
#pyplot.plot(xaxis, yaxis)
pyplot.show()
