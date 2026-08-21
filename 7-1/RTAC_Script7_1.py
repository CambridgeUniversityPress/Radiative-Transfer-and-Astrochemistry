
import numpy as np
import matplotlib.pyplot as plt

#initial abundances
F1old = 1.0
F2old = 0.0
F3old = 0.0

#time step
dt = 0.1
#dt = 1.e-3

#these arrays will store the time, and abundances, over time
tarray = []
F1array = []
F2array = []
F3array = []

#start time
tstart=0.0

#end time
tend= 10.0 

#set up an array of times to evaluate the abundances over
t = np.arange(tstart, tend, dt)

print(t)

#t = 0.0 

#now loop through all the times
for i in range(0, len(t)):
    #these are our discretised abundance differential equations
    F1new = F1old - dt*F1old
    F2new = F2old + dt*(-0.5*F2old + F1old) 
    F3new = F3old + (dt/2.0)*F2old
    #we have update all abundances, so those now become the "old" ones
    F1old = F1new
    F2old = F2new
    F3old = F3new
    #OPTIONAL: add our values to the arrays so we can plot them
    tarray.append(t[i]+dt)
    F1array.append(F1new)
    F2array.append(F2new)
    F3array.append(F3new)
    #print out the time and various abundances
    print(f"{tarray[i]:.3f}" +  " & " + f"{F1new:.5f}" +  " & " + f"{F2new:.5f}" +  " & " + f"{F3new:.5f}" + " \\\ ")

#This loads in the KROME result. Just remove this if using for some other application    
data=np.loadtxt("fort.66", skiprows=1)


#plot results
plt.plot(data[:,0], data[:,1],label="FK_1")
plt.plot(data[:,0], data[:,2],label="FK_2")
plt.plot(data[:,0], data[:,3],label="FK_3")

plt.ylabel("Abundance [arbitrary]", fontsize=16)
plt.xlabel("Time [arbitrary]", fontsize=16)
    
plt.scatter(tarray, F1array)
plt.scatter(tarray, F2array)
plt.scatter(tarray, F3array)
plt.legend()
plt.xlim(0,5)
plt.savefig("Comparison_0p5step.pdf")
plt.show()



