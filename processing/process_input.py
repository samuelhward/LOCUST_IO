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

def QTP_calc(Q=None,T=None,P=None):
    """
    generic script to solve Q=dT/dP when given 2/3 variables (no integration constants)

    notes:
        used to calculate the missing quantity out of Q, toroidal or poloidal flux
        http://theory.ipp.ac.cn/~yj/research_notes/tokamak_equilibrium/node11.html
        returns profiles normalised to zero at origin 

        Q - first order at sides, second order at in the centre
        T - second order composite trapezium rule
        P - second order composite trapezium rule
    """

    if Q is None: #need to calculate Q
        
        ''' old numpy version
        dP=np.diff(P)
        dP=np.append(dP,dP[0]) #use same difference as the start of array, since this usually set smaller to counteract np.gradient's higher error
        Q=np.gradient(T,dP)'''

        P[P==0.0]=0.0001 #replace zero values with small numbers to stop divide by zero
        Q=np.gradient(T,P)

        return Q

    elif T is None: #need to calculate T

        T=scipy.integrate.cumtrapz(y=Q,x=P,initial=0) #0 here should be toroidal flux at the magnetic axis

        return T

    elif P is None: #need to caclulate P

        Q[Q==0.0]=0.0001 #replace zero values with small numbers to stop divide by zero
        P=scipy.integrate.cumtrapz(y=1/Q,x=T,initial=0) #0 here should be poloidal flux at the magnetic axis

        return P






#################################

##################################################################

###################################################################################################