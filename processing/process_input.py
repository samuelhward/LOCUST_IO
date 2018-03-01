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

###################################################################################################
#Preamble

import scipy.integrate
import random
import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d #import 3D plotting axes
from matplotlib import cm #get colourmaps

pi=np.pi


###################################################################################################
#Main Code


################################################################## Supporting functions

def interpolate_2D(x,y,X_grid,Y_grid,Z_grid):
    """
    interpolate a 2D grid to positions x,y

    notes:
        X_grid, Y_grid are 1D grid axes
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    """

    #XXX need this? X_grid,Y_grid=np.meshgrid(X_grid,Y_grid)
    interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function='cubic',smooth=0)
    interpolated_values=interpolator(x,y)

    return interpolated_values

def interpolate_1D(x,X_line,Y_line):
    """
    interpolate a 1D function to position x

    notes:
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    """

    interpolator=scipy.interpolate.Rbf(X_line,Y_line,function='cubic',smooth=0)
    interpolated_values=interpolator(x)

    return interpolated_values





################################################################## Specific functions

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
        Q[Q==0.0]=0.000001 #replace zero values with small numbers to stop divide by zero
        P=scipy.integrate.cumtrapz(y=1/Q,x=T,initial=0) #0 here should be poloidal flux at the magnetic axis

        return P

def fpolrz_calc(fpol_1d,psirz,simag,sibry,B_vacuum_toroidal_centre,R_centre):
    """
    interpolates 1D poloidal flux function 'fpol' onto an r,z grid given  

    notes:
        assumes fpol is defined from magnetic axis (simag) to plasma boundary (sibry)

        if the grid extends outside of the plasma radius, where poloidal flux function is not defined, then use a given vacuum toroidal field and work
        backwardsto get the value of the flux function (since f is constant in vacuum)
    """

    flux_pol_1d=np.linspace(simag,sibry,len(fpol_1d)) #build uniform flux grid

    nw=len(psirz[:,0]) #determine size of computational domain
    nh=len(psirz[1,:])

    fpolrz=np.zeros((nw,nh)) #initialise 2D grid

    for w in np.arange(nw): #loop over 2D grid
        for h in np.arange(nh):
            
            if np.abs(psirz[w,h])<np.abs(sibry) and np.abs(psirz[w,h])>np.abs(simag): #if within boundaries of 1D fpol
                fpolrz[w,h]=interpolate_1D(psirz[w,h],flux_pol_1d,fpol_1d) #calling this interpolater function constantly like this is very inefficient
            else:
                fpolrz[w,h]=B_vacuum_toroidal_centre*R_centre

    return fpolrz

def B_calc(psirz,fpolrz,R_1D,Z_1D): #XXX check all the ordering
    """
    calculates r, phi, z axisymmetric magnetic field at coordinates R_1D, Z_1D

    notes:
        returns B(r,z,i), a 3D array holding component i of B field at grid indices r, z
        i=0 - B_r
        i=1 - B_toroidal
        i=2 - B_z
    """
    
    gradient=np.gradient(psirz,R_1D,Z_1D) #calculate gradient along both axes (gradient[i] is 2D)

    one_R=1.0/R_1D

    B_z=gradient[0]*(-1.0)*one_R #XXX check here that we're dividing by R in the correct order (should break if not since length of one_R should only be equal to nw)...maybe one of these references to one_R should be transposed?
    B_r=gradient[1]*one_R
    B_tor=fpolrz*one_R

    B_field=np.array([[[B_r[w,h],B_tor[w,h],B_z[w,h]] for w in range(len(R_1D))] for h in range(len(Z_1D))],ndmin=3)

    return B_field

def transform_marker_velocities(r,phi,z,pitch,speed,R_1D,Z_1D,B_field,conversion='guiding_centre'):
    """
    converts marker velocities to r,phi,z coordinates using random sampling

    notes:
        currently uses 2D B_field (toroidally symmetric)
        pitch here is v|| / v
    """

    if conversion=='guiding_centre':

        R_max=np.max(R_1D) #calculate these in advance to speed up
        R_min=np.min(R_1D)
        Z_max=np.max(Z_1D)
        Z_min=np.min(Z_1D)
        escapees=[] #to hold indices of rogue particles 

        if np.max(r)>R_max or np.min(r)<R_min or np.max(z)>Z_max or np.min(z)<Z_min: #check particles are within domain
            print("WARNING: markers found outside of computational boundary - removing them from list")

            for index,(r_particle,z_particle) in enumerate(zip(r,z)): #if some are rogue, begin search and removal
                if r_particle>R_max or r_particle<R_min or z_particle>Z_max or z_particle<Z_min:
                
                    escapees.extend([index]) #mark index for deletion (do not delete yet otherwise dynamically changing array length in loop!)
                
            r_trim=np.delete(r,escapees) #new arrays with rogue particles trimmed out (still works if escapees=[])
            phi_trim=np.delete(phi,escapees)
            z_trim=np.delete(z,escapees)
            pitch_trim=np.delete(pitch,escapees)
            speed_trim=np.delete(speed,escapees)

        if escapees: 
            B_field_r=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,0]) #interpolate B field components to each particle 
            B_field_phi=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,1])
            B_field_z=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,2])
        else:
            B_field_r=interpolate_2D(r,z,R_1D,Z_1D,B_field[:,:,0]) #interpolate B field components to each particle 
            B_field_phi=interpolate_2D(r,z,R_1D,Z_1D,B_field[:,:,1])
            B_field_z=interpolate_2D(r,z,R_1D,Z_1D,B_field[:,:,2])
        
        #calculate velocity magnitudes (speed=magnitude of velocity)
        v_par_magnitude=pitch*speed
        v_perp_magnitude=np.sqrt(1-((v_par_magnitude**2)/(speed**2)))

        #resolve the parallel components (V||i = some_constant * Bi for all i)
        some_constant=

        #resolve the perpendicular components (requires random sampling)

        rand_1=[random.uniform(0,1) for particle in range(len(r)-len(escapees))] #need to generate two random numbers per particle
        rand_2=[random.uniform(0,1) for particle in range(len(r)-len(escapees))]

        return v_r,v_phi,v_z

#################################

##################################################################

###################################################################################################