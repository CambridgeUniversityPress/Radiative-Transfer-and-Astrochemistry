#
# Script 5-1 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# This script illustrates the extent of the habitable zone, defined as the region in which water is in liquid form
# This is based on radiative equilibrium and ignores factors such as the greenhouse effect. 
#

#load libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

#1au in m 
autom = 1.5e11

#solar radius in m
rsol=696.24e6


#set up the range of stellar temperatures for the slider
starTMin = 2000.0
starTMax = 8000.0
delta_starT = 10.0
starT0=5780.0

#set up the range of planetary albedos for the slider
AlphaMin = 0.0
AlphaMax = 0.99
delta_Alpha = 0.01
Alpha0=0.3

#set up the plot
fig, ax1 = plt.subplots(figsize=(8,7))
plt.subplots_adjust(left=0.25, bottom=0.35)

axcolor = 'lightgoldenrodyellow'
axstarT = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axAlpha = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sstarT=Slider(axstarT, 'Stellar Temperature, K', starTMin, starTMax, valinit=starT0, valstep=delta_starT)
sAlpha=Slider(axAlpha, 'Albedo, \u03B1', AlphaMin, AlphaMax, valinit=Alpha0, valstep=delta_Alpha)

#set up the initial extent of the habitable zone 
Rstar = 1.0*rsol
planetD0=1.0
Tplanet = starT0*(Rstar/2.0/planetD0/autom)**(0.5)*(1.0-Alpha0)**(0.25)

#Outer edge of habitable zone, where T=100 degrees centigrade (373.15K)
DistOut=0.5*Rstar/(((373.15/starT0)**4/(1.0-Alpha0))**0.5)/autom

#Inner edge of habitable zone, where T=0 degrees centrigrate (273.15K)
DistIn=0.5*Rstar/(((273.15/starT0)**4/(1.0-Alpha0))**0.5)/autom

n, radii = 100, [DistIn, DistOut]
theta = np.linspace(0, 2*np.pi, n, endpoint=True)
xs = np.outer(radii, np.cos(theta))
ys = np.outer(radii, np.sin(theta))
xs[1,:] = xs[1,::-1]
ys[1,:] = ys[1,::-1]

ax1.set_xlim(-5.5, 5.5)
ax1.set_ylim(-5.5, 5.5)
ax1.set_xlabel("x, au", fontsize=12)
ax1.set_ylabel("y, au", fontsize=12)
ax1.text(-3.8,3.0,"Radiative Equilibrium Habitable Zone, No Atmospheres")

#place markers for the star and the positions of some solar system planets
ax1.scatter(0,0, marker='*', s=300, color="orange")
ax1.scatter(0.3,0.0, marker='o', s=5, color="gray", label="Mercury")
ax1.scatter(1.0,0.0, marker='o', s=10, color="blue", label="Earth")
ax1.scatter(5.2,0.0, marker='o', s=50, color="red", label="Jupiter")

ax1.legend(loc=3)

#add the habitable zone ring in green
ring=ax1.fill_between(np.ravel(xs), np.ravel(ys), color="green", alpha=0.5, linewidth=0.0)

#update the plot as the slider(s) are updated
def update(val):
    Rstar=1.0*rsol*(sstarT.val/5780.)**(0.5)
    DistOut=0.5*Rstar/(((373.15/sstarT.val)**4/(1.0-sAlpha.val))**0.5)/autom
    DistIn=0.5*Rstar/(((273.15/sstarT.val)**4/(1.0-sAlpha.val))**0.5)/autom    

    ax1.cla()
    ax1.set_xlim(-5.5, 5.5)
    ax1.set_ylim(-5.5, 5.5)
    ax1.set_xlabel("x, au", fontsize=12)
    ax1.set_ylabel("y, au", fontsize=12)
    ax1.scatter(0,0, marker='*', s=300, color="orange")
    ax1.scatter(0.3,0.0, marker='o', s=5, color="gray", label="Mercury")
    ax1.scatter(1.0,0.0, marker='o', s=10, color="blue", label="Earth")
    ax1.scatter(5.2,0.0, marker='o', s=50, color="red", label="Jupiter")
    ax1.text(-3.8,3.0,"Radiative Equilibrium Habitable Zone, No Atmospheres")
    n, radii = 100, [DistIn, DistOut]
    theta =	np.linspace(0, 2*np.pi, n, endpoint=True)
    xs = np.outer(radii, np.cos(theta))
    ys = np.outer(radii, np.sin(theta))
    xs[1,:] = xs[1,::-1]
    ys[1,:] = ys[1,::-1]
    ring=ax1.fill(np.ravel(xs), np.ravel(ys), color="green", alpha=0.5, linewidth=0.0)
    ax1.legend(loc=3)
    
sstarT.on_changed(update)
sAlpha.on_changed(update)

plt.show()


