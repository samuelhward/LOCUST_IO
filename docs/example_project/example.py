#!/usr/bin/env python

#example.py

'''
Samuel Ward
24/01/2018
----
very simple example project script for manipulating the sample input/output data
---
notes:

---
'''



##################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import sys
    import numpy
    import context
    import matplotlib.pyplot as plt
    from classes.input_classes import number_density
    from classes.output_classes import orbits
    import matplotlib
    from matplotlib import cm
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main Code

print('''\n\n\n\
88                                                              88              
88                                                ,d            ""              
88                                                88                        
88  ,adPPYba,   ,adPPYba, 88       88 ,adPPYba, MM88MMM         88  ,adPPYba,   
88 a8"     "8a a8"     "" 88       88 I8[    ""   88            88 a8"     "8a  
88 8b       d8 8b         88       88  `"Y8ba,    88            88 8b       d8  
88 "8a,   ,a8" "8a,   ,aa "8a,   ,a88 aa    ]8I   88,           88 "8a,   ,a8"  
88  `"YbbdP"'   `"Ybbd8"'  `"YbbdP'Y8 `"YbbdP"'   "Y888 ________88  `"YbbdP"' \n\n\n''')



#read in the sample input data
my_Ne=number_density.Number_Density(ID='LOCUST_IO sample Ne profile',data_format='LOCUST',filename='sample_ne_profile.file',species='electrons')
#read in the sample output data
my_orbit=orbits.Orbits(ID='LOCUST_IO sample orbits',data_format='LOCUST',filename='sample_orbits.file')

#add some new data to our orbits object
my_orbit.set(number_of_orbits=1)

#what do the orbits look like?
my_orbit.plot(colmap=settings.cmap_g) #use default green colour map - see settings.py

#what if I want to plot my orbits and the number density in one figure?

#start by creating some axis objects
fig,(ax1,ax2)=plt.subplots(1,2)
#then let us add an x,y plot of the orbits to this figure on the first axis
mycmap=matplotlib.cm.get_cmap('jet') #make a new colourmap
my_orbit.plot(axes=['R','Z'],colmap=mycmap,ax=ax1,fig=fig)
#now add a plot of number density on the second axis
my_Ne.plot(ax=ax2,fig=fig)
#if I do not like the default axis labels, I can always overwrite them
ax1.set_xlabel('a nice trajectory')
#display the plot
plt.show()

#################################

##################################################################

###################################################################################################