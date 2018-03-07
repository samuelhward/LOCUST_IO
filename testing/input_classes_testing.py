#!/usr/bin/env python

#input_classes_testing.py

"""
Samuel Ward
02/11/2017
----
Small examples to test the input_classes.py methods
---
notes:
---
"""



##################################################################
#Preamble

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
	import context
except:
	raise ImportError("ERROR: local context.py not found - please create one and add it to this directory\nreturning\n")
	sys.exit(1)

try: 
	import input_classes as ic
except:
	raise ImportError("ERROR: LOCUST_IO/classes/input_classes.py could not be imported!\nreturning\n")
	sys.exit(1)

try: 
	import process_input as prin
except:
	raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")
	sys.exit(1)

pi=np.pi

##################################################################
#Main Code


################################# test basic functionality




################################# test calculator functions with analytical cases
#QTP_calc()

	#define an analytical solution
	P=np.linspace(0,8,num=1000) 
	T=np.sin(P)
	Q=np.cos(P)

	#QTP results
	P_recover=QTP_calc(T=T,Q=Q)
	Q_recover=QTP_calc(T=T,P=P)
	T_recover=QTP_calc(Q=Q,P=P)

	plt.plot(P,'r-')
	plt.plot(Q,'g-')
	plt.plot(T,'b-')
	plt.show()

	plt.plot(P_recover,'r-')
	plt.plot(Q_recover,'g-')
	plt.plot(T_recover,'b-')
	plt.show()

	#residuals
	plt.plot(P_recover-P,'r-')
	plt.plot(Q_recover-Q,'g-')
	plt.plot(T_recover-T,'b-')
	plt.show()

#interp1D()

	x=np.linspace(0,10,9)
	y=x+5.0
	x_new=pi
	y_new=prin.interp1D(x_new,x,y)
	print('y_new={}'.format(y_new))

#interp2D()

#fpolrz_calc()

#B_calc()

#transform_marker_velocities()


#################################

##################################################################

###################################################################################################