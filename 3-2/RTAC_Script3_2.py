#
# Script 3-2 for the textbook "Foundations of Radiative Transfer & Astrochemistry"  
# This script is a possible solution for Quick Exercise 3.3. # 
# It estimates pi using a Monte Carlo approach (see 3.4.1)
#

import random


#number of monte carlo samples (smaller value = quicker run time, but lower accuracy).
nTests = 10000000
#initialise the number of random points that land in or outside the circle (in the square)
nCirc = 0.0

for i in range(0, nTests):
    #draw a random point
    randomX = random.uniform(-1, 1)
    randomY = random.uniform(-1, 1)
    #determine whether it is in or outside of the circle.
    R = (randomX**2 + randomY**2)**0.5
    if(R <= 1.0):
        #the point is inside the circle
        nCirc = nCirc + 1.0

print("pi estimate is ", (4.0*nCirc/nTests))
