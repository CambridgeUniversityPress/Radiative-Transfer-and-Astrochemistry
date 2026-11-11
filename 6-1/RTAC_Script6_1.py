#
# Script 6-1 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# Plots synthetic CO datacubes of a protoplanetary disc at various inclinations
# It also plots the spatially unresolved line profile
# Includes sliders that change the inclination and velocity channel
#
# Note that the synthetic datacubes take time to calculate so were pre-computed
# with the TORUS radiative transfer code (Harries et al. 2019, A&C, 27, 63)

#import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from astropy.io import fits
import matplotlib.colors as colors


#this is the initial dataset that will be plotted
cube70="intensity_dusttest.fits"

#open the fits cube
cube70_hdu_list = fits.open(cube70)
cube70_image_data = cube70_hdu_list[0].data

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,7))
#plt.subplots_adjust(left=0.25, bottom=0.35)
plt.subplots_adjust(left=0.08, right=0.97, bottom=0.35, wspace=0.35)
axcolor = 'lightgoldenrodyellow'

axvel = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axinc = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

#set the velocity slider range and interval
velMin=-5.0
velMax=5.0
nV=80.
delta_vel=(velMax-velMin)/nV
vel0 = 0.0

#set the inclination slider range and interval
incMin=0.0
incMax=90.0
delta_inc=10.0
inc0=70.0

#set up the sliders
svel=Slider(axvel, 'Velocity, km/s', velMin, velMax, valinit=vel0, valstep=delta_vel)
sinc=Slider(axinc, 'Inclination, degrees', incMin, incMax, valinit=inc0, valstep=delta_inc)

ax1.imshow(cube70_image_data[40,:,:], norm=colors.LogNorm())
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("CO channel map, v = %.1f km/s" % vel0)
ax2.set_title("Spatially unresolved line profile")


lineProf=np.zeros(80)
vel=np.zeros(80)
for i in range(0,80):
    #This is doing a crude continuum subtraction by subtracting the emission in a ~line free channel
    lineProf[i]=np.sum(cube70_image_data[i,:,:])-np.sum(cube70_image_data[2,:,:])
    vel[i]=-5.0+delta_vel*i
lineProf[:]=lineProf[:]/max(lineProf[:])
ax2.set_ylim(0.25, 1.1)    
ax2.set_xlim(-4.9, 4.9)    
ax2.set_xlabel("Velocity, km/s")
ax2.set_ylabel("Normalised intensity")
ax2.plot(vel, lineProf)

#update the plot as the slider is moved
def update(val):
    ax1.cla()
    ax2.cla()
    incInt=int(sinc.val)
    fileString="intensity_dusttest_"+str(incInt)+".fits"
    cube_hdu_list = fits.open(fileString)
    cube_image_data = cube_hdu_list[0].data
    index = min(int( (svel.val-velMin)/delta_vel), 79)
    ax1.imshow(cube_image_data[index,:,:], norm=colors.LogNorm())
    ax1.set_ylabel("")
    ax1.set_xlabel("")
    lineProf=np.zeros(80)
    vel=np.zeros(80)    
    for i in range(0,80):
        lineProf[i]=np.sum(cube_image_data[i,:,:])-np.sum(cube_image_data[2,:,:])
        vel[i]=-5.0+delta_vel*i
    lineProf[:]=lineProf[:]/max(lineProf[:])
    ax2.set_ylim(0.25, 1.1)
    ax2.set_xlim(-4.9, 4.9)    
    ax2.set_xlabel("Velocity, km/s")
    ax2.set_ylabel("Normalised intensity")
    ax1.set_title("CO channel map, v = %.1f km/s" % svel.val)
    ax2.set_title("Spatially unresolved line profile")
    ax2.plot(vel, lineProf)

svel.on_changed(update)
sinc.on_changed(update)

plt.show()

