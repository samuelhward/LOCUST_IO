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
    import input_classes
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    import output_classes
except:
    raise ImportError("ERROR: LOCUST_IO/classes/output_classes.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    import plot_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_input.py could not be imported!\nreturning\n")
    sys.exit(1)

try: 
    import plot_output
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_output.py could not be imported!\nreturning\n")
    sys.exit(1)



##################################################################
#Main Code


#read in the sample input data
my_Ne=input_classes.Number_Density(ID='LOCUST_IO sample Ne profile',data_format='ASCII',input_filename='sample_ne_profile.file',properties='electrons')

#see what it looks like
plot_input.plot_number_density(my_Ne)



#read in the sample output data
my_orbits=output_classes.Orbits(ID='LOCUST_IO sample orbits',data_format='ASCII',input_filename='sample_orbits.file')

#plot particle 0's trajectory in the R,Z plane
plot_output.plot_orbits(my_orbits,[0])


#################################

##################################################################

###################################################################################################