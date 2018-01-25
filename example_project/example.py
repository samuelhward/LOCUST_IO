#!/usr/bin/env python

#example.py

"""
Samuel Ward
24/01/2018
----
very simple example project script for converting a set of IDSs into the LOCUST inputs
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





##################################################################
#Main Code



#read in all the required data

equilibrium=input_classes.Equilibrium(ID='that awesome Jet equilibrium',data_format='IDS',shot=1,run=7) 

particles=input_classes.Beam_Deposition(ID='that awesome ASDEX-U NBI deposition',data_format='IDS',shot=2,run=8)

temperature_i=input_classes.Temperature(ID='that awesome KSTAR ion temperature profile',data_format='IDS',shot=3,run=9,properties='ions') 

temperature_e=input_classes.Temperature(ID='that awesome MAST-U electron temperature profile',data_format='IDS',shot=4,run=10,properties='electrons')

density_i=input_classes.Number_Density(ID='that awesome DIII-D ion density profile',data_format='IDS',shot=5,run=11,properties='ions')

density_e=input_classes.Number_Density(ID='that awesome ITER electron density profile',data_format='IDS',shot=6,run=12,properties='electrons')


#dump all the data to the LOCUST formats 

equilibrium.dump_data(data_format='GEQDSK',output_filename='my_equilibrium.eqdsk') 

particles.dump_data(data_format='ASCII',output_filename='ptcles.dat')

temperature_i.dump_data(data_format='ASCII',output_filename='profile_Ti.dat')

temperature_e.dump_data(data_format='ASCII',output_filename='profile_Te.dat')

density_i.dump_data(data_format='ASCII',output_filename='profile_ni.dat')

density_e.dump_data(data_format='ASCII',output_filename='profile_ne.dat')



#################################

##################################################################

###################################################################################################