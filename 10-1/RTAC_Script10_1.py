#
# Script 10-1 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# This python script solves for the position, mass, radius and mass loss rate of an
# externally irradiated globule over time
# 

#import libraries
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
mu=0.5          #mean particle mass in ionised hydrogen gas
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



#evolve the cloud over time
while t < tend:
    #unattenuated ionising flux
    Fo = StellarNly/4.0/pi/distance**2

    #mass loss rate (from Mellema et al. 2008)
    Fprime = Fo / (1.0+ (alphaB*Fo*Rglob/3.0/cii**2))**0.5
    Mdot=-Fprime*mu*mH*pi*(pi*Rglob/2.0)**2

    #update mass first
    Mlast=Mglob
    Mglob=Mlast+(Mdot*dt)

    #check if the globule has been almost completely dispersed
    if Mglob <= 0.001*iniMglob:
        print("Globule fully photoevaporated at t = %.3f Myr" % (t/Myr))
        break

    #update globule bulk velocity
    Rdotlast=Rdot
    Rdot = Rdot - ((dt *  Mdot *  cii)/Mglob)

    #move the cloud
    distance = distance + dt*(Rdot)

    #update the globule radius
    Rlast=Rglob
    Rglob = (Mglob *3./4./pi/rhoGlob)**(1./3.)
    step=step+1
    #append the results to an array that will be plotted every 100 time steps
    if(step%100):
        tArray.append(t/Myr)
        MdotArray.append(-Mdot*yr/msol)
        velArray.append(Rdot/1.e5)
        distArray.append(distance/pc)
        Rarray.append(Rglob/pc)        
        Marray.append(Mglob/msol)
    t = t + dt

#Plot the evolution of various globule properties in one figure
fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
panels = [
    (axs[0,0], np.log10(MdotArray), r"Mass loss rate, log$_{10}(M_\odot\,{\rm yr}^{-1})$"),
    (axs[0,1], velArray,  "Velocity, km/s"),
    (axs[0,2], distArray, "Distance from UV source, pc"),
    (axs[1,0], Marray,    r"Mass, $M_\odot$"),
    (axs[1,1], Rarray,    "Radius, pc"),
]
for ax, y, label in panels:
    ax.plot(tArray, y, linewidth=2)
    ax.set_ylabel(label, fontsize=13)
    ax.tick_params(labelsize=11)
for ax in (axs[1,0], axs[1,1], axs[0,2]):
    ax.set_xlabel("Time, Myr", fontsize=13)
axs[0,2].xaxis.set_tick_params(labelbottom=True)
axs[1,2].axis("off")   #hide the unused sixth panel
fig.tight_layout()
fig.savefig("GlobuleEvolution.pdf")
plt.show()
