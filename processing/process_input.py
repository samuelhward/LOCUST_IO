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
import scipy



##################################################################
#Main Code


def calc_Q_tor_pol(Q=None,T=None,P=None):
    """
    calculates the missing quantity out of Q, toroidal or poloidal flux
    notes:
        feed this function a string 
        http://theory.ipp.ac.cn/~yj/research_notes/tokamak_equilibrium/node11.html
    """

    if not Q: #need to calculate Q
        
        Q=np.gradient(T,P)

        return Q

    elif not T: #need to calculate T

        T=scipy.integrate.cumtrapz(y=Q,x=P,initial=0) #XXX  or this could be just Q=T/P because our segments are straight?

        return T

    elif not P: #need to caclulate P

            #this is calculating the dx
        return P



#################################

##################################################################

###################################################################################################