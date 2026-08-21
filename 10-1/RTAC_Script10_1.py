#
# This python script solves for the position, mass, radius and mass loss rate of an externally irradiated globule. 
# T. J. Haworth. 
#

import numpy as np
import os
import matplotlib.pyplot as plt


#constants
au=1.5e13       #1au in cm
pc=3.08e18      #parsec in cm
msol=1.989e33   #solar mass in g
pi=3.14159
Myr=3.14e13
kyr=3.14e10
yr=3.14e7       #year in seconds
cii=12.8e5      #ionised gas sound speed in cm/s
mu=0.5
mH=1.67e-24     #hydrogen mass in grams


###initial conditions###

#initial globule radius
iniRglob=0.6*pc

#initial globule mass
iniMglob=10.3*msol

#timescale of the calculation
tstart=0.0
tend = 0.5*Myr

#initial distance of the globule from the UV source
distance = 2.0*pc

#UV source ionising photons per second
StellarNly=1.e49

#case B recombination coefficient
alphaB=2.7e-13

#initialise some parameters
Mglob=iniMglob
Rglob=iniRglob
Rdotdot=0.0
Rdot=0.0
Mdot=0.0
globVel=0.0
Rdotlast=0.0
Rlast=Rglob
Mlast=Mglob

t=tstart
#we will do 100 thousand time steps (you can try changing this) 
dt = (tend-tstart)/1.e5
step=0

Marray=[]
Rarray=[]
distArray=[]
velArray=[]
tArray=[]
MdotArray=[]
rhoGlob= Mglob*3./4./pi/Rglob**3
firstTime=False

twoPhase = False

#evolve the cloud over time
while t < tend:
    #see Reiter et al. 2020MNRAS.497.3351R
    #unattenuated ionising flux
    Fo = StellarNly/4.0/pi/distance**2
    #mass loss rate (from Mellema et al. 2008)
    Fprime = Fo / (1.0+ (alphaB*Fo*Rglob/3.0/cii**2))**0.5
    Mdot=-Fprime*mu*mH*pi*(pi*Rglob/2.0)**2

    #update mass first
    Mlast=Mglob
    Mglob=Mlast+(Mdot*dt)

    Rdotlast=Rdot
#    Rdot = (Rdot - dM*(dRdot-cii))/Mglob
    Rdot = Rdot - ((dt *  Mdot * (Rdot -- cii))/Mglob)
    #move the cloud
    distance = distance + dt*(Rdot)

    Rlast=Rglob
    #can have two phase collapse if desired (see Reiter + )
    if(twoPhase): 
        if(t < 0.3*Myr):
            Rglob = iniRglob*(Mglob/iniMglob)**(3.)
            r_RDI = Rglob
            rhoGlob = Mglob*3./4./pi/Rglob**3
        else:
            Rglob = (Mglob *3./4./pi/rhoGlob)**(1./3.)
    else:
        Rglob = (Mglob *3./4./pi/rhoGlob)**(1./3.)
    step=step+1
    #append the results to an array that will be plotted every 100 time steps
    if(step%100):
        tArray.append(t/Myr)
        MdotArray.append(-Mdot*yr/msol)
        velArray.append(Rdot/1.e5)
        distArray.append(distance/pc)
#        Rarray.append(Rglob/au/1.e3)
        Rarray.append(Rglob/pc)        
        Marray.append(Mglob/msol)
    t = t + dt


#Plot the evolution of various globule properties 
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Time, Myr", fontsize=16)
plt.ylabel("Globule mass loss rate, log$_{10}M_\odot\,yr^{-1}$", fontsize=16)
plt.plot(tArray[:], np.log10(MdotArray[:]), linewidth = 3, label="model")
plt.legend()
plt.savefig("Mdot_2part.pdf")
plt.show()



plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Time, Myr", fontsize=16)
plt.ylabel("Globule velocity, km/s", fontsize=16)
plt.plot(tArray[:], velArray[:], linewidth = 3, label="model")
plt.legend()
plt.savefig("velocity_2part.pdf")
plt.show()


plt.ylim(0,10)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Time, Myr", fontsize=16)
plt.ylabel("Globule distance from UV source, pc", fontsize=16)
plt.plot(tArray[:], distArray[:], linewidth = 3, label="model")
plt.legend()
plt.savefig("position_2part.pdf")
plt.show()


plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Time, Myr", fontsize=16)
plt.ylabel("Globule Mass, $M_\odot$", fontsize=16)
plt.plot(tArray[:], Marray[:], linewidth = 3, label="model")
plt.legend()
plt.savefig("Mass_2part.pdf")
plt.show()


plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel("Time, Myr", fontsize=16)
plt.ylabel("Globule Radius, pc", fontsize=16)
plt.plot(tArray[:], Rarray[:], linewidth = 3, label="model")
plt.legend()
plt.savefig("Radius_2part.pdf")
plt.show()


