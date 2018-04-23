#!/usr/bin/env python

#example.py

"""
Samuel Ward
24/01/2018
----
very simple example project script for manipulating the sample input/output data
---
notes:

---
"""



##################################################################
#Preamble

import sys
import numpy

try:
    import context
except:
    raise ImportError("ERROR: local context.py not found - please create one and add it to this directory\nreturning\n")
    sys.exit(1)

try: 
    from classes.input_classes import number_density
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/number_density.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    from classes.output_classes import orbits
except:
    raise ImportError("ERROR: LOCUST_IO/classes/output_classes/orbits.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    from processing import plot_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_input.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    from processing import plot_output
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_output.py could not be imported!\nreturning\n")
    sys.exit(1)



##################################################################
#Main Code

print('''\n\n\n\
88                                                          88              
88                                                ,d        ""              
88                                                88                    
88  ,adPPYba,   ,adPPYba, 88       88 ,adPPYba, MM88MMM     88  ,adPPYba,   
88 a8"     "8a a8"     "" 88       88 I8[    ""   88        88 a8"     "8a  
88 8b       d8 8b         88       88  `"Y8ba,    88        88 8b       d8  
88 "8a,   ,a8" "8a,   ,aa "8a,   ,a88 aa    ]8I   88,       88 "8a,   ,a8"  
88  `"YbbdP"'   `"Ybbd8"'  `"YbbdP'Y8 `"YbbdP"'   "Y888 ____88  `"YbbdP"' \n\n\n''')




#read in the sample input data
my_Ne=number_density.Number_Density(ID='LOCUST_IO sample Ne profile',data_format='LOCUST',filename='sample_ne_profile.file',properties='electrons')

#see what it looks like
plot_input.plot_number_density(my_Ne)

#read in the sample output data
my_orbit=orbits.Orbits(ID='LOCUST_IO sample orbits',data_format='LOCUST',filename='sample_orbits.file')

#plot particle 0's trajectory in the R,Z plane
plot_output.plot_orbits(my_orbit)


#################################

##################################################################

###################################################################################################