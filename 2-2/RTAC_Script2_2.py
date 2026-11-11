#
# Script 2-2 for the textbook "Foundations of Radiative Transfer & Astrochemistry"
# This script is a simple ray tracing code for the example of an optically thin emitting sphere. 
# The script includes a slider allowing the user to change the grid resolution and the density 
# structure of the sphere. 
#


#load libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

#constants/parameters
#ncells=80              #number of cells in our grid
hConst=6.626205e-34     #Planck constant
kB=1.380626e-23         #Boltzmann constant
c=2.99792458e8         #Speed of light


#setup a 3D density structure

def setupGrid(index, sphereRadius, ncells):
#have a sphere of radius "sphereRadius", in units of parsecs
#size of a cell for ncells is dL
    dL = (2.4*sphereRadius)/float(ncells)

#set up some coordinate arrays (x,y,z) and some empty density/temperature arrays    
    x=np.arange(-1.2*sphereRadius, 1.2*sphereRadius, dL)
    y=np.arange(-1.2*sphereRadius, 1.2*sphereRadius, dL)
    z=np.arange(-1.2*sphereRadius, 1.2*sphereRadius, dL)
    rho = np.zeros(shape=(int(ncells),int(ncells),int(ncells)))
    Temperature = np.zeros(shape=(int(ncells),int(ncells),int(ncells)))
    
    #loop through all cells and populate with a density and temperature                         
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                thisR=(x[i]**2+y[j]**2+z[k]**2)**0.5
                if thisR <= sphereRadius:
                    #inside the sphere
                    rho[i,j,k]=10.0*(0.001+thisR/sphereRadius)**(index) #density profile scales with radius to the power of some index
                    Temperature[i,j,k]=1.e4 #Temperature in K, this is ionised gas
                else:
                    #outside the sphere
                    rho[i,j,k]=1.e-10
                    Temperature[i,j,k]=10.0
    return x, y, z, rho, Temperature, dL


#return free-free Bremsstrahlung emissivity for a given density and temperature
def freefreeEmissivity(thisRho, thisT, freq):
    thisAlpkk = 1.0
    eta= thisRho**2*thisAlpkk*np.exp(-hConst*freq/kB/thisT)*(2.0*hConst*freq**3/c**2)
    return eta

#Perhaps just make this a number that the user can specify
#Or maybe get them to specify the source function?
def getKappa(density):
    Kappa = 0.1
    return Kappa
    
#routine to do the ray tracing
def TraceRaysToObserver(pixIDx, pixIDy, pixIDz, density, temperature, freq, ncells, dL):
    #set up the intensity array with zero everywhere. 
    intensity = np.zeros(shape=(int(ncells),int(ncells)))
    #loop over x/y and trace a ray through z to work out the total intensity in each cell

    for i in range(len(pixIDx)):
        for j in range(len(pixIDy)):
            for k in range(len(pixIDz)):
                intensity[i,j] = intensity[i,j] + dL * freefreeEmissivity(density[i,j,k], temperature[i,j,k], freq)
    maxval=np.amax(intensity)
    intensity[:,:] = intensity[:,:]/maxval
    return intensity



#main execution

ncells=64

sphereRadius=1.0
freq=1.5e9 #Hz (20cm)
pi = 3.14159

EmissionTimesVolume = (4.0/3.0)*pi*sphereRadius**3*freefreeEmissivity(10.0, 1.e4, freq)

print("Emission times volume ", EmissionTimesVolume)

#setup the grid
xcoord, ycoord, zcoord, density, T, dL=setupGrid(0.0, sphereRadius, ncells)

#do the ray tracing
iArray=TraceRaysToObserver(xcoord, ycoord, zcoord, density, T, freq, ncells, dL)


#lets also compare a cut along the mid-plane with the analytic solution
analyticEmissivity=freefreeEmissivity(10., 1.e4, freq)
analytic=np.zeros(len(xcoord))
for i in range(0, len(xcoord)):
    if((xcoord[i]**2)**0.5< sphereRadius):        
        analytic[i]=2.0*analyticEmissivity*sphereRadius*(1.0-(xcoord[i]**2)/(sphereRadius**2))**0.5

maxval=np.amax(analytic)
analytic=analytic/maxval        

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
plt.subplots_adjust(left=0.25, bottom=0.35)


#
#SLIDERS
#
axcolor = 'lightgoldenrodyellow'

pixMin=10
pixMax=128
delta_pix=2
pix0=64

axpix = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
spix = Slider(axpix, 'Pixels per side', pixMin, pixMax, valinit=pix0, valstep=delta_pix)

powMin=-0.5
powMax=3.0
delta_pow=0.1
pow0=0.0 #uniform

axpow = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
spow = Slider(axpow, r'$\gamma$ in $\rho = \rho_0(R/R_s)^\gamma$', powMin, powMax, valinit=pow0, valstep=delta_pow)
ax1.set_xlabel("")
ax1.set_ylabel("")

toplot=ax1.pcolor(xcoord, ycoord, iArray)
ax1.set_xlim(-1.15,1.15)
ax1.set_ylim(-1.15,1.15)
#cbar_ax = fig.add_axes()
cb=plt.colorbar(toplot, ax=[ax1], location='left')

analyticplot, =ax2.plot(xcoord, analytic, color='red', label="Analytic Uniform Sphere")
profile, =ax2.plot(xcoord, iArray[:,int(ncells/2)], label="Numerical Model")
ax2.legend()
ax2.set_xlabel("Offset from sphere centre", fontsize=12)
ax2.set_ylabel("Normalized intensity", fontsize=12)


def update(val):
    #setup the grid
    xcoord, ycoord, zcoord, density, T, dL=setupGrid(spow.val, sphereRadius, spix.val)
    #do the ray tracing
    iArray=TraceRaysToObserver(xcoord, ycoord, zcoord, density, T, freq, spix.val, dL)
    profile.set_xdata(xcoord)
    profile.set_ydata(iArray[:,int(spix.val/2)])
    toplot=ax1.pcolor(xcoord, ycoord, iArray)
    

spix.on_changed(update)
spow.on_changed(update)
plt.show()


