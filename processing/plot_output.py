#plot_output.py

'''
Samuel Ward
25/1/2018
----
plotting routines for LOCUST output data
---
notes:
---
'''

##################################################################
#Preamble
import matplotlib.pyplot as plt
import numpy as np



##################################################################
#Main Code


def plot_orbits(some_orbits,particles):
    """
    simple orbits plot in the R,Z plane
     
    notes:
        particles needs to be iterable list of particle numbers
    """
    
    for particle in particles:

        plt.plot(some_orbits['orbits'][:,particle,0],some_orbits['orbits'][:,particle,2])

    plt.show()


#################################

##################################################################

###################################################################################################