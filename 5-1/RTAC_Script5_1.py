#
# Script 5-1 for the textbook "Foundations of Radiative Transfer & Astrochemistry" 
#
# A simple script that plots the Boltzmann distribution with and without statistical weights.
# Includes an interactive slider to modify the temperature
#
# data from LAMDA database (https://home.strw.leidenuniv.nl/~moldata/)
#

#load libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

#Boltzmann constant
kB = 1.38e-23

levelsCO=np.arange(0,41,1)

#
#energies are in cm^{-1} units in LAMDA database
#in the script these are converted to Joules using the factor 1.9863e-23 which is the Planck constant * speed of light
#
energiesCO=[0.0, 3.845, 11.5350, 23.0695, 38.4481, 57.6704, 80.7354, 107.6424, 138.3904, 172.9780, 211.4041, 253.6672, 299.7656, 349.6975, 403.4612, 461.0545, 522.4751, 587.7209, 656.7892, 729.6774, 806.3828, 886.9024, 971.2332, 1059.3718, 1151.3150, 1247.0592, 1346.6008, 1449.9360, 1557.0607, 1667.9709, 1782.6624, 1901.1308, 2023.3715, 2149.3799, 2279.1511, 2412.6801, 2549.9620, 2690.9913, 2835.7627, 2984.2707, 3136.5095]

#statistical weights of the CO levels
gCO=[1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0,15.0,17.0,19.0,21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0, 57.0, 59.0, 61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0, 79.0, 81.0]


#Works out the partition function, including the statistical weights (degeneracy factors)
def calcPartitionFunction(energies, degeneracies, temperature):
    PartitionFunction = 0.0
    for i in range(0,len(energies)):
        PartitionFunction=PartitionFunction + (degeneracies[i]*np.exp(-energies[i]*1.98630e-23/kB/temperature))
    return PartitionFunction

#calculate the level populations including the statistical weights (degeneracy factors)
def calcPops(energies, degeneracies, temperature, partitionFunction):
    Pop=np.zeros(len(energies))
    for i in range(0,len(energies)):    
        Pop[i]=(1.0/partitionFunction)*degeneracies[i]*np.exp(-energies[i]*1.98630e-23/kB/temperature)  
    return Pop

#works out the partition function ignoring the statistical weights/degeneracy factors
def calcPartitionFunctionND(energies, degeneracies, temperature):
    PartitionFunction = 0.0
    for i in range(0,len(energies)):
        PartitionFunction=PartitionFunction + (1.0*np.exp(-energies[i]*1.98630e-23/kB/temperature))
    return PartitionFunction

#calculate the level populations ignoring the statistical weights/degeneracy factors
def calcPopsND(energies, degeneracies, temperature, partitionFunction):
    Pop=np.zeros(len(energies))
    for i in range(0,len(energies)):    
        Pop[i]=(1.0/partitionFunction)*1.0*np.exp(-energies[i]*1.98630e-23/kB/temperature)  
    return Pop

#initial temperature in the script
temperature = 10.0

#Calculate the CO level populations including statistical weights
COpart = calcPartitionFunction(energiesCO, gCO, temperature)
COpops = calcPops(energiesCO, gCO, temperature, COpart)

#Calculate the CO level populations ignoring statistical weights 
COpartND = calcPartitionFunctionND(energiesCO, gCO, temperature)
COpopsND = calcPopsND(energiesCO, gCO, temperature, COpartND)


#plotting with an interactive slider
fig, ax1 = plt.subplots(figsize=(8,6))
plt.subplots_adjust(left=0.25, bottom=0.35)
axcolor = 'lightgoldenrodyellow'
axtemp = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

#these values describe the range and interval of the temperature slider
delta_temp=5.0
tempMin=10.0
tempMax=1000.0
temp0=10.0
sTemp=Slider(axtemp, 'Temperature, K', tempMin, tempMax, valinit=temp0, valstep=delta_temp)

ax1.set_xlabel(r"CO rotational energy level, $J$", fontsize=16)
ax1.set_ylabel(r"$p_J$", fontsize=16)
COlevPlot1, =ax1.plot(levelsCO, COpops, label="With statistical weights")
COlevPlot2, =ax1.plot(levelsCO, COpopsND, label="Without statistical weights")
ax1.set_yscale("log")
ax1.set_ylim(0.001, 1.0)
ax1.set_xlim(0,30)

#update the plot as the slider is updated
def update(val):
    COpart = calcPartitionFunction(energiesCO, gCO, sTemp.val)
    COpops = calcPops(energiesCO, gCO, sTemp.val, COpart)
    COpartND = calcPartitionFunctionND(energiesCO, gCO, sTemp.val)
    COpopsND = calcPopsND(energiesCO, gCO, sTemp.val, COpartND)    
    COlevPlot1.set_ydata(COpops)
    COlevPlot2.set_ydata(COpopsND)
    
sTemp.on_changed(update)
ax1.legend()
plt.show()

