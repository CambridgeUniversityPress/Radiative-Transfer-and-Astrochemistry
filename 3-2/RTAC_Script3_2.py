import numpy as np
import matplotlib.pyplot as plt
import random


#number of monte carlo samples 
nTests = 100000
#initialise the number of random points that land in or outside the circle (in the square)
nCirc = 0.0
nSquare = 0.0

for i in range(0, nTests):
    #draw a random point
    randomX = random.uniform(-1, 1)
    randomY = random.uniform(-1, 1)
    #determine whether it is in or outside of the circle.
    R = (randomX**2 + randomY**2)**0.5
    if(R <= 1.0):
        nCirc = nCirc + 1.0
    else:
        nSquare = nSquare + 1.0

print("pi estimate is ", (4.0*nCirc/nSquare))
