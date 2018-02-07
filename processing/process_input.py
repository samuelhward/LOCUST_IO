#process_input.py

'''
Samuel Ward
25/1/2018
----
processing routines for LOCUST input data
---
notes:
    https://stackoverflow.com/questions/9455111/python-define-method-outside-of-class-definition
---
'''

##################################################################
#Preamble

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d #import 3D plotting axes
from matplotlib import cm #get colourmaps
import scipy.integrate



##################################################################
#Main Code
'''
def calc_gradient(x=None,y=None):
	"""
	function to calculate gradient for arbitrarily-spaced x,y data

	notes:
		first order at sides, second order at in the centre
	"""

	dx=np.diff(x)
	dy=np.diff(y)
	gradient=[]

	gradient.append(dy[0]/dx[0])

	for i+1 in range(len(x)):

		gradient.append((dy[i]-dy[i-1])/(2*(dx[i]-dx[i-1])))

	gradient.append(dy[-1]/dx[-1])

	return gradient
'''
def calc_Q_tor_pol(Q=None,T=None,P=None):
    """
    calculates the missing quantity out of Q, toroidal or poloidal flux
    notes:
        http://theory.ipp.ac.cn/~yj/research_notes/tokamak_equilibrium/node11.html

        Q - first order at sides, second order at in the centre
        T - second order composite trapezium rule
        P - second order composite trapezium rule
    """

    if Q is None: #need to calculate Q
        
        dP=np.diff(P)
        dP=np.append(dP,dP[0]) #use same difference as the start of array, since this usually set smaller to counteract np.gradient's higher error
        Q=np.gradient(T,dP)

        return Q

    elif T is None: #need to calculate T

        T=scipy.integrate.cumtrapz(y=Q,x=P,initial=0) #0 here should be toroidal flux at the magnetic axis

        return T

    elif P is None: #need to caclulate P

    	P=scipy.integrate.cumtrapz(y=1/Q,x=T,initial=0) #0 here should be poloidal flux at the magnetic axis

        return P


#################################

##################################################################

###################################################################################################