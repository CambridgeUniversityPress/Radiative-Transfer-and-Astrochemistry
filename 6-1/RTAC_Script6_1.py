
#
# EMRA L4/L5 CO synthetic observations of a protoplanetary disc. 
#
# Started by T. J. Haworth October 2020
#
#
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from astropy.io import fits
import matplotlib.colors as colors

#cube40="intensity_test_40.fits"
cube40="intensity_dusttest.fits"

cube40_hdu_list = fits.open(cube40)
cube40_image_data = cube40_hdu_list[0].data

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,7))
plt.subplots_adjust(left=0.25, bottom=0.35)

axcolor = 'lightgoldenrodyellow'

axvel = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axinc = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

velMin=-5.0
velMax=5.0
nV=80.
delta_vel=(velMax-velMin)/nV
vel0 = 0.0

incMin=0.0
incMax=90.0
delta_inc=10.0
inc0=40.0

svel=Slider(axvel, 'Velocity, km/s', velMin, velMax, valinit=vel0, valstep=delta_vel)
sinc=Slider(axinc, 'Inclination, degrees', incMin, incMax, valinit=inc0, valstep=delta_inc)

ax1.imshow(cube40_image_data[40,:,:], norm=colors.LogNorm())
ax1.set_ylabel("")
ax1.set_xlabel("")
#plt.colorbar()


lineProf=np.zeros(80)
vel=np.zeros(80)
for i in range(0,80):
    lineProf[i]=np.sum(cube40_image_data[i,:,:])-np.sum(cube40_image_data[2,:,:])
    vel[i]=-5.0+delta_vel*i
lineProf[:]=lineProf[:]/max(lineProf[:])
ax2.set_ylim(0.25, 1.1)    
ax2.set_xlim(-4.9, 4.9)    
ax2.set_xlabel("velocity")
ax2.set_ylabel("I")
ax2.plot(vel, lineProf)

def update(val):
    ax1.cla()
    ax2.cla()
    incInt=int(sinc.val)
    fileString="intensity_dusttest_"+str(incInt)+".fits"
    cube_hdu_list = fits.open(fileString)
    cube_image_data = cube_hdu_list[0].data
    index = int( (velMin + svel.val)/delta_vel)
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
    ax2.set_xlabel("velocity")
    ax2.set_ylabel("I")
    ax2.plot(vel, lineProf)

svel.on_changed(update)
sinc.on_changed(update)

plt.show()

