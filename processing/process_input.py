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

import sys 

try:
    import scipy.integrate
    import numpy as np
    import matplotlib.pyplot as plt
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
    
try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from constants import *
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


###################################################################################################
#Main Code

################################################################## Specific functions

def QTP_calc(Q=None,T=None,P=None):
    """
    generic script to solve Q=dT/dP when given 2/3 variables (no integration constants)

    notes:
        used to calculate the missing quantity out of Q, toroidal or poloidal flux
        http://theory.ipp.ac.cn/~yj/research_notes/tokamak_equilibrium/node11.html
        returns profiles normalised to zero at origin 
    args:
        Q - first order at sides, second order at in the centre
        T - second order composite trapezium rule
        P - second order composite trapezium rule
    """

    if Q is None: #need to calculate Q
        
        ''' old numpy version
        dP=np.diff(P)
        dP=np.append(dP,dP[0]) #use same difference as the start of array, since this usually set smaller to counteract np.gradient's higher error
        Q=np.gradient(T,dP)'''
        print("QTP_calc - calculating Q profile")
        P[P==0.0]=0.0001 #replace zero values with small numbers to stop divide by zero
        Q=np.gradient(T,P)
        print("QTP_calc - finished calculating Q profile")
        return Q

    elif T is None: #need to calculate T
        print("QTP_calc - calculating toroidal flux")
        T=scipy.integrate.cumtrapz(y=Q,x=P,initial=0) #0 here should be toroidal flux at the magnetic axis
        print("QTP_calc - finished calculating toroidal flux")
        return T

    elif P is None: #need to caclulate P
        print("QTP_calc - calculating poloidal flux")
        Q[Q==0.0]=0.000001 #replace zero values with small numbers to stop divide by zero
        P=scipy.integrate.cumtrapz(y=1/Q,x=T,initial=0) #0 here should be poloidal flux at the magnetic axis
        print("QTP_calc - finished calculating poloidal flux")
        return P

def fpolrz_calc(some_equilibrium):
    """
    interpolates 1D poloidal flux function fpol onto an r,z grid given  

    notes:
        fpol is defined over the poloidal flux points flux_pol

        assumes flux_pol and fpol are defined from magnetic axis to plasma boundary
        
        if the grid extends outside of the plasma radius, where poloidal flux function is not defined, then use a given vacuum toroidal field and work
        backwards to get the value of the flux function (since f is constant in vacuum)
    """
    print("fpolrz_calc - calculating 2D flux function")

    fpolrz=np.zeros((some_equilibrium['nR_1D'],some_equilibrium['nZ_1D'])) #initialise 2D grid
    fpolrz_interpolator=processing.utils.interpolate_1D(some_equilibrium['flux_pol'],some_equilibrium['fpol']) #generate the interpolator

    for w in np.arange(some_equilibrium['nR_1D']): #loop over 2D grid
        for h in np.arange(some_equilibrium['nZ_1D']):
            
            if np.abs(some_equilibrium['psirz'][w,h])<=np.max(np.abs(some_equilibrium['flux_pol'])) and np.abs(some_equilibrium['psirz'][w,h])>=np.min(np.abs(some_equilibrium['flux_pol'])): #if the poloidal flux function is defined for this poloidal flux
                fpolrz[w,h]=fpolrz_interpolator(some_equilibrium['psirz'][w,h])
            else:
                fpolrz[w,h]=some_equilibrium['bcentr']*some_equilibrium['rcentr']

    print("fpolrz_calc - finished calculating 2D flux function")
    return fpolrz

def B_calc(some_equilibrium): #ordering OK - for Z[X,Y].shape=(2,3), gradient(Z,X,Y) --> gradient[0] along X, gradient[1] along Y
    """
    calculates r, phi, z axisymmetric magnetic field at coordinates R_1D, Z_1D

    notes:
        see Wesson, based on RH coordinate system [R,Phi,Z]
        returns B(r,z,i), a 3D array holding component i of B field at grid indices r, z
        i=0 - B_r
        i=1 - B_toroidal
        i=2 - B_z

        remember - A[i,j,k] is an array nested such as i*[j*[k*[some_number]]]
                 - [row,column] in numpy
    """
    
    if 'fpolrz' not in some_equilibrium.data: #could have a try except statement here to accommodate for both dictionaries (which may not contain .data) and LOCUST_IO objects 
        print("WARNING: B_calc - fpolrz missing in equilibrium object")
        some_equilibrium['fpolrz']=fpolrz_calc(some_equilibrium) #if fpolrz is missing calculate it

    print("B_calc - calculating 2D magnetic field")

    gradient=np.gradient(some_equilibrium['psirz'],some_equilibrium['R_1D'],some_equilibrium['Z_1D']) #calculate gradient along both axes (gradient[i] is 2D)

    one_R=1.0/some_equilibrium['R_1D']
    one_R=one_R[:,np.newaxis]

    B_field_Z=(gradient[0])*one_R
    B_field_R=(gradient[1])*one_R*(-1.0)
    B_field_tor=some_equilibrium['fpolrz']*one_R

    B_field=np.array([[[B_field_R[w,h],B_field_tor[w,h],B_field_Z[w,h]] for h in range(len(some_equilibrium['Z_1D']))] for w in range(len(some_equilibrium['R_1D']))],ndmin=3)
    
    print("B_calc - finished calculating 2D magnetic field")
    
    return B_field

def mag_axis_calc(some_equilibrium,index=False):
    """
    returns psirz index of magnetic axis 
    
    args:
        index - if True, returns index of axis otherwise return location in metres 
    notes:
        works by calculating a contour at 10% poloidal flux and averaging the position of the contour
        XXX does not work well
    """

    Z_2D,R_2D=np.meshgrid(some_equilibrium['Z_1D'],some_equilibrium['R_1D'])
    contour=plt.contour(R_2D,Z_2D,np.abs(some_equilibrium['psirz']),levels=[0.0,np.amin(np.abs(some_equilibrium['psirz']))+0.1*np.amax(np.abs(some_equilibrium['psirz']))])  
    contour=contour.collections[1].get_paths()[0]
    contour=contour.vertices
    r=contour[:,0]
    z=contour[:,1]
    r_av=np.mean(r)
    z_av=np.mean(z)

    if index is True: #convert to psirz indices
        r_av=(np.abs(some_equilibrium['R_1D']-r_av)).argmin()
        z_av=(np.abs(some_equilibrium['Z_1D']-z_av)).argmin()
    
    return r_av,z_av

def transform_marker_velocities(r=None,phi=None,z=None,pitch=None,speed=None,R_1D=None,Z_1D=None,B_field=None,conversion='guiding_centre'):
    """
    returns 6D marker positions in r,phi,z coordinates

    XXX untested

    notes:
        currently uses 2D B_field (toroidally symmetric)
        pitch means v|| / v
        if some particles lie outside the computational domain then these are removed from the calculation
    """

    if conversion=='guiding_centre':

        if all(arg is not None for arg in [r,phi,z,pitch,speed,R_1D,Z_1D,B_field]): #check we have enough data to do this type of conversion

            print("transform_marker_velocities - converting beam_deposition from guiding_centre coordinates")

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
            
            r_trim=np.delete(r,escapees) #delete the rogue particles (if escapees empty this does nothing)
            phi_trim=np.delete(phi,escapees)
            z_trim=np.delete(z,escapees)
            pitch_trim=np.delete(pitch,escapees)
            speed_trim=np.delete(speed,escapees)

            print('transform_marker_velocities - generating B field interpolators')
            B_field_r_interpolator=processing.utils.interpolate_2D(R_1D,Z_1D,B_field[:,:,0]) #generate the interpolator functions
            B_field_tor_interpolator=processing.utils.interpolate_2D(R_1D,Z_1D,B_field[:,:,1])
            B_field_z_interpolator=processing.utils.interpolate_2D(R_1D,Z_1D,B_field[:,:,2])
            print('transform_marker_velocities - interpolating B field to particle positions')
            B_field_r=B_field_r_interpolator(r_trim,z_trim) #interpolate B field components to each particle 
            B_field_tor=B_field_tor_interpolator(r_trim,z_trim)
            B_field_z=B_field_z_interpolator(r_trim,z_trim)
            print('transform_marker_velocities - done interpolating')
            B_field_magnitude=np.sqrt((B_field_r**2)+(B_field_tor**2)+(B_field_z**2)) #B magnitude at each particle

            #calculate velocity magnitudes (speed)
            print('transform_marker_velocities - calculating parallel velocities')
            v_par_magnitude=np.abs(pitch_trim*speed_trim)
            v_perp_magnitude=np.sqrt((speed_trim**2)-(v_par_magnitude**2))

            #resolve the parallel components (V||i = some_constant * Bi for all i=r,phi,z)
            some_constant=np.abs(v_par_magnitude/B_field_magnitude) #keep this absolute, since v sign introduced in multiplication with B components
            v_par_r=some_constant*B_field_r
            v_par_tor=some_constant*B_field_tor
            v_par_z=some_constant*B_field_z
            
            print('transform_marker_velocities - sampling perpendicular velocities')
            rand_1=np.array([np.random.uniform(0.0,v) for v in v_perp_magnitude]) #need to generate two random velocities per particle
            rand_2=np.array([np.random.uniform(0.0,np.sqrt((v**2)-(rand**2))) for v,rand in zip(v_perp_magnitude,rand_1)])
            rand_3=np.sqrt((v_perp_magnitude**2)-((rand_1**2)+(rand_2**2))) 

            rand_1*=np.random.choice([-1,1],size=len(rand_1)) #generate the magnitude and sign separately (.uniform(-v,v) would not produce 0.0 with correct probability)
            rand_2*=np.random.choice([-1,1],size=len(rand_2))
            rand_3*=np.random.choice([-1,1],size=len(rand_3))

            v_perp_r=rand_1 #NOTE could randomise this order of assignment of velocities to each rand_ number here
            v_perp_tor=rand_2
            v_perp_z=rand_3

            v_r=v_par_r+v_perp_r
            v_tor=v_par_tor+v_perp_tor
            v_z=v_par_z+v_perp_z

            print("transform_marker_velocities - finished transforming marker velocities")

            return r_trim,phi_trim,z_trim,v_r,v_tor,v_z

        else:
            print("ERROR: insufficient data supplied to transform_marker_velocities for {} conversion".format(conversion))
            return

#################################

##################################################################

###################################################################################################