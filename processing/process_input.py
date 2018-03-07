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

try:
    import random
    import scipy.interpolate
    import scipy.integrate
    import numpy as np
except:
    raise ImportError("ERROR: initial imported modules not found!\nreturning\n")
    sys.exit(1)
    
pi=np.pi


###################################################################################################
#Main Code


################################################################## Supporting functions

def interpolate_2D(x,y,X_axis,Y_axis,Z_grid):
    """
    interpolate a 2D grid to positions x,y

    notes:
        X_axis, Y_axis are 1D grid axes
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    """

    X_grid,Y_grid=np.meshgrid(X_axis,Y_axis)
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
        print("QTP calculating Q profile")
        P[P==0.0]=0.0001 #replace zero values with small numbers to stop divide by zero
        Q=np.gradient(T,P)
        print("QTP finished calculating Q profile")
        return Q

    elif T is None: #need to calculate T
        print("QTP calculating toroidal flux")
        T=scipy.integrate.cumtrapz(y=Q,x=P,initial=0) #0 here should be toroidal flux at the magnetic axis
        print("QTP finished calculating toroidal flux")
        return T

    elif P is None: #need to caclulate P
        print("QTP calculating poloidal flux")
        Q[Q==0.0]=0.000001 #replace zero values with small numbers to stop divide by zero
        P=scipy.integrate.cumtrapz(y=1/Q,x=T,initial=0) #0 here should be poloidal flux at the magnetic axis
        print("QTP finished calculating poloidal flux")
        return P

def fpolrz_calc(pol_flux_1d,pol_flux_func_1d,psirz,B_vacuum_toroidal_centre,R_centre):
    """
    interpolates 1D poloidal flux function 'pol_flux_func_1d' onto an r,z grid given  

    notes:
        pol_flux_func_1d is defined over the poloidal flux points 'pol_flux_1d'

        assumes pol_flux_1d and pol_flux_func_1d are defined from magnetic axis to plasma boundary
        
        if the grid extends outside of the plasma radius, where poloidal flux function is not defined, then use a given vacuum toroidal field and work
        backwards to get the value of the flux function (since f is constant in vacuum)
    """
    print("fpolrz_calc 2D flux function")

    nw=len(psirz[:,0]) #determine size of computational domain
    nh=len(psirz[1,:])

    fpolrz=np.zeros((nw,nh)) #initialise 2D grid

    for w in np.arange(nw): #loop over 2D grid
        for h in np.arange(nh):
            
            if np.abs(psirz[w,h])<=np.max(np.abs(pol_flux_1d)) and np.abs(psirz[w,h])>=np.min(np.abs(pol_flux_1d)): #if the poloidal flux function is defined for this poloidal flux
                fpolrz[w,h]=interpolate_1D(psirz[w,h],pol_flux_1d,pol_flux_func_1d) #NOTE calling this interpolater function constantly like this is very inefficient
            else:
                fpolrz[w,h]=B_vacuum_toroidal_centre*R_centre

    print("finished fpolrz_calc 2D flux function")

    return fpolrz

def B_calc(psirz,fpolrz,R_1D,Z_1D): #XXX check all the ordering- if i swap ordering in B_field_i calcs with that break things? I think the order is bad - but our data is 101x101 so need a GEQDSK which is rectangular in domain to test this
    """
    calculates r, phi, z axisymmetric magnetic field at coordinates R_1D, Z_1D

    notes:
        returns B(r,z,i), a 3D array holding component i of B field at grid indices r, z
        i=0 - B_r
        i=1 - B_toroidal
        i=2 - B_z
    """
    
    print("B_calc 2D magnetic field")

    gradient=np.gradient(psirz,R_1D,Z_1D) #calculate gradient along both axes (gradient[i] is 2D) XXX is R_1D and Z_1D in the right order?

    one_R=1.0/R_1D

    B_field_z=gradient[0]*(-1.0)*one_R #XXX check here that we're dividing by R in the correct order (should break if not since length of one_R should only be equal to nw)...maybe one of these references to one_R should be transposed?
    B_field_r=gradient[1]*one_R
    B_field_tor=fpolrz*one_R

    B_field=np.array([[[B_field_r[w,h],B_field_tor[w,h],B_field_z[w,h]] for w in range(len(R_1D))] for h in range(len(Z_1D))],ndmin=3)
    
    print("finished B_calc 2D magnetic field")
    
    return B_field

def transform_marker_velocities(r=None,phi=None,z=None,pitch=None,speed=None,R_1D=None,Z_1D=None,B_field=None,conversion='guiding_centre'):
    """
    returns 6D marker positions in r,phi,z coordinates

    notes:
        currently uses 2D B_field (toroidally symmetric)
        pitch means v|| / v
        if some particles lie outside the computational domain then these are removed from the calculation
    """
    print("transform_marker_velocities")

    if conversion=='guiding_centre':

        if all(arg is not None for arg in [r,phi,z,pitch,speed,R_1D,Z_1D,B_field]): #check we have enough data to do this type of conversion

            print("converting beam_deposition from guiding_centre coordinates")

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

            print('interpolating B field')

            B_field_r=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,0]) #interpolate B field components to each particle 
            B_field_tor=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,1])
            B_field_z=interpolate_2D(r_trim,z_trim,R_1D,Z_1D,B_field[:,:,2])
            B_field_magnitude=np.sqrt((B_field_r**2)+(B_field_tor**2)+(B_field_z**2)) #B magnitude at each particle

            print('done interpolating')
            print('calculating parallel velocities')

            #calculate velocity magnitudes (speed)
            v_par_magnitude=np.abs(pitch_trim*speed_trim)
            v_perp_magnitude=np.sqrt((speed_trim**2)-(v_par_magnitude**2))

            #resolve the parallel components (V||i = some_constant * Bi for all i=r,phi,z)
            some_constant=np.abs(v_par_magnitude/B_field_magnitude) #keep this absolute, since v sign introduced in multiplication with B components
            v_par_r=some_constant*B_field_r
            v_par_tor=some_constant*B_field_tor
            v_par_z=some_constant*B_field_z
            
            print('sampling perpendicular velocities')

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

            print("finished transform_marker_velocities")

            return r_trim,phi_trim,z_trim,v_r,v_tor,v_z

        else:
            print("ERROR: insufficient data supplied to transform_marker_velocities for {} conversion".format(conversion))


#################################

##################################################################

###################################################################################################