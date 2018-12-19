#processing.utils.py
 
"""
Samuel Ward
02/11/2017
----
supporting functions for LOCUST input/output classes
---
usage:
    see README.md for usage
 
notes:         
---
"""


###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway
try:
    import numpy as np
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import processing.process_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")

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

################################################################################################### Supporting functions
 
def none_check(ID,LOCUST_input_type,error_message,*args):
    """
    generic function for checking if a None value appears in *args
 
    notes:
        message should be something specific to the section of code which called none_check
    """
    if all(arg is not None for arg in args):
        return False
    else:
        print("WARNING: none_check returned True (LOCUST_input_type={LOCUST_input_type}, ID={ID}): {message}".format(LOCUST_input_type=LOCUST_input_type,ID=ID,message=error_message))
        return True
          
def get_next(obj):
    """
    generic object iterator
 
    notes:
        every time get_next() is called, it will advance the counter one
    """
    pyVer = sys.version_info[0]
    if pyVer == 2:
        return obj.next() #NOTE how exactly does .next() work? will it accept tok in toklist values from file_numbers() and allow file_numbers to then read in the next bunch (since yield will pause the While True loop?)
    else:
        return next(obj)
     
def dict_set(data,**kwargs):
    """
    generalised function (based upon safe_set) to set values in dictionary 'data' after checking source data exists 
    
    notes:
    usage:
        dict_set(data,some_key=some_source,some_other_key=[1,2,3,4]) to set multiple values simultaneously
        dict_set(data,**{'some_key':100,'some_other_key':200}) equally
    """
    keys=kwargs.keys()
    values=kwargs.values()
    
    for key,value in zip(keys,values): #loop through kwargs
        if key is not None and value is not None: #use 'is' since we may want to set data to 0
            data[key]=value

def safe_set(target,source):
    """
    generalised function to set a target value to source if it exists
    
    notes:
    """
    if source is not None: #use 'is' since we may want to set data to 0
        target=source

 
def fortran_string(number_out,length,decimals=None,exponential=True):
    """
    produces a string stream of specified length, padded with space at the start

    notes:
        supposed to be a quick fix to fortran format descriptors
    args:
        number_out - variable holding number to write out 
        length - total number of characters for output string stream including white space 
        decimals - required for non-ints, how many decimal places the number is written to 
        exponential - true for output to exponential format, otherwise formatted like float
    """

    if decimals is not None:
        if exponential is True:
            output_string='{:.{decimals}e}'.format(number_out,decimals=decimals)        
        else:
            output_string='{:.{decimals}f}'.format(number_out,decimals=decimals)        
    else: #assume integer if decimals not specified
        output_string=str(number_out)

    number_spaces=length-len(output_string) #caclulate amount of padding needed from length of the number_out
    if number_spaces>=0:
        string_stream='{}{}'.format(' '*number_spaces,output_string) #construct the string stream
        return string_stream
    else:
        print('ERROR: fortran_string() too many decimal places for requested string length!')


def sort_arrays(main_array,*args):
    """
    sort an arbitrary number of arrays in parallel (*args) according to main_array

    notes:
        assumes all arrays are numpy arrays
    usage:
        to sort ascending according to array x, 
        x_sorted,y_sorted,z_sorted,...=sort_arrays(x,y,z,...)
    """

    sorted_indices=main_array.argsort() #get the order of sorted indices
    
    main_array=main_array[sorted_indices] #convert all the arrays
    returned_arrays=[main_array] #add arrays to list
    for array in args:
        array=array[sorted_indices]
        returned_arrays.append(array)

    return returned_arrays

def angle_pol(R_major_width,R,Z):
    """
    returns poloidal angle

    notes:
    """

    R_minor=R-R_major_width
    angle=np.arctan2(Z,R_minor)

    return angle

def interpolate_2D(X_axis,Y_axis,Z_grid,function='multiquadric',type='RBS',smooth=0,rect_grid=True):
    """
    generate a 2D grid interpolator
    notes:
        keep as separate functions so can freely swap out interpolation method
        RBF - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
            - high memory overhead, most accurate
        RBS - https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html#scipy.interpolate.RectBivariateSpline
              https://scipython.com/book/chapter-8-scipy/examples/two-dimensional-interpolation-with-scipyinterpolaterectbivariatespline/
    args:
        X_axis - 1D x-axis
        Y_axis - 1D y-axis
        Z_grid - 2D z-axis
        function - spline function
        type - name of interpolation function to use
        smooth - level of smoothing
        rect_grid - toggle whether X_axis, Y_axis are regular rectangular grid edges or arbitrary cordinates 
    usage:
        my_interpolator=interpolate_2D(X_axis,Y_axis,data)
        interpolated_value=my_interpolator(x,y)
    """

    try:
        import scipy.interpolate 
    except:
        raise ImportError("ERROR: interpolate_2D could not import scipy.interpolate module!\nreturning\n")
        return

    if type=='RBF':
        if rect_grid:
            Y_grid,X_grid=np.meshgrid(Y_axis,X_axis) #swap since things are defined r,z 
        else:
            Y_grid,X_grid=Y_axis,X_axis
        interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function=function,smooth=smooth)
    
    elif type=='RBS':
        interpolator=scipy.interpolate.RectBivariateSpline(X_axis,Y_axis,Z_grid) #normally order is other way in RBS but I have swapped my axes

    return interpolator

def interpolate_1D(X_axis,Y_axis,function='cubic',type='RBF',smooth=0):
    """
    generate a 1D line interpolator

    notes:
        keep as separate functions so can freely swap out interpolation method
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    args:
        X_axis - 1D x-axis
        Y_axis - 1D y-axis
        function - spline function
        type - name of interpolation function to use
        smooth - level of smoothing
    usage:
        my_interpolator=interpolate_1D(X_axis,data)
        interpolated_value=my_interpolator(x)        
    """

    try:
        import scipy.interpolate 
    except:
        raise ImportError("ERROR: interpolate_1D could not import scipy.interpolate module!\nreturning\n")
        return

    if type=='RBF':
        interpolator=scipy.interpolate.Rbf(X_axis,Y_axis,function=function,smooth=smooth)
    elif type=='interp1d':
        interpolator=scipy.interpolate.interp1d(X_axis,Y_axis,kind=function)

    return interpolator

def RphiZ_to_XYZ(R,phi,RH=True):
    """
    converts R,phi positions to X,Y 

    args:
        R - R coordinates to convert
        phi - Toroidal angle coordinates to convert
        RH - toggle right-handed R phi Z coordinate system
    notes:
    """

    if RH is not True:
        phi=2.*pi-phi

    X=R*np.cos(phi)
    Y=R*np.sin(phi)

    return X,Y

def XYZ_to_RphiZ(X,Y):
    """
    converts X,Y positions to R,phi 
    """

    phi=np.arctan2(Y,X)
    R=X*np.cos(phi)+Y*np.sin(phi)

    return R,phi

def RZ_to_Psi(R_coordinate,Z_coordinate,equilibrium):
    """
    maps RZ point to corresponding poloidal flux

    args:
        R_coordinate - corresponding 1D R coordinates  
        Z_coordinate - corresponding 1D Z coordinates
        equilibrium - equilibrium object with psi over 2D grid
    notes:
        returned flux is in [Wb/rad]
    """

    interpolator_psi=interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['psirz'])
    psi_at_coordinate=[]

    for R,Z in zip(R_coordinate,Z_coordinate): #do element-wise to avoid implicitly interpolating 2x1D arrays onto a single 2D grid
        psi_point=interpolator_psi(R,Z)
        psi_at_coordinate.append(np.squeeze(psi_point))

    return np.asarray(psi_at_coordinate)

def dot_product(vec1,vec2):
    """
    performs dot products of two vectors sharing same orthogonal basis
    
    notes:
    vec1,vec2 - lists of orthogonal vector components  
    """

    dot_product=0.
    for vec_1_component,vec_2_component in zip(vec1,vec2):
        dot_product+=vec_1_component*vec_2_component

    return dot_product

def pitch_calc_2D(some_particle_list,some_equilibrium):
    """
    calculates the pitch angles of a given set of particle positions for axisymmetric B_field 

    notes:
        XXX untested
    args:
        some_equilibrium - equilibrium object
        some_particle_list - generic particle list with positions some_particle_list['R'],some_particle_list['Z']
    """

    eq=copy.deepcopy(some_equilibrium)

    if 'fpolrz' not in some_equilibrium.data.keys(): #first go about calculating the magnetic field if missing
        eq.set(fpolrz=processing.process_input.fpolrz_calc(eq)) 
    if 'B_field' not in some_equilibrium.data.keys():
        eq.set(B_field=processing.process_input.B_calc(eq))

    B_field_R_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,0]) #generate interpolators for B
    B_field_tor_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,1])
    B_field_Z_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,2])

    print("pitch_calc_2D - interpolating B_field to particles\n")

    #to stop memory errors, doing element-wise

    Br_at_particles=[]
    Btor_at_particles=[]
    Bz_at_particles=[]
    B_at_particles=[]

    for R,Z in zip(some_particle_list['R'],some_particle_list['Z']):

        Br=B_field_R_interpolator(R,Z)
        Btor=B_field_tor_interpolator(R,Z)
        Bz=B_field_Z_interpolator(R,Z)
        Br_at_particles.append(Br) #interpolate B to particles
        Btor_at_particles.append(Btor)
        Bz_at_particles.append(Bz)
        B_at_particles.append(np.sqrt(Br**2+Btor**2+Bz**2)) #calculate |B| at particles - quicker than interpolating

    print("pitch_calc_2D - finished interpolating B_field to particles\n")

    Br_at_particles=np.asarray(Br_at_particles)
    Btor_at_particles=np.asarray(Btor_at_particles)
    Bz_at_particles=np.asarray(Bz_at_particles)

    Br_at_particles/=B_at_particles #calculate B_hat
    Btor_at_particles/=B_at_particles
    Bz_at_particles/=B_at_particles

    V_parallel=dot_product([some_particle_list['V_R'],some_particle_list['V_tor'],some_particle_list['V_Z']],[Br_at_particles,Btor_at_particles,Bz_at_particles])

    V=np.sqrt(some_particle_list['V_R']**2+some_particle_list['V_tor']**2+some_particle_list['V_Z']**2)
    V_pitch=V_parallel/V

    return V_pitch

def get_dfn_point(dfn,type='LOCUST',**kwargs):
    """
    returns magnitude of dfn at a point closest to that supplied

    args:
        dfn - distribution_function object
        type - type of distribution function
        kwargs - should define desired points in all possible dimensions to sample at
    usage:
        my_values=get_dfn_point(my_dfn,E=[1,2,3],V_pitch=[-1,0,1],R=[1,2,3],Z=[0,0,0],P=[-pi,-pi,-pi])
        my_values=get_dfn_point(my_dfn,'TRANSP',E=[1,2,3],V_pitch=[-1,0,1],R=[1,2,3],Z=[0,0,0]) (for TRANSP)
    notes:
    """

    if type=='TRANSP': #TRANSP has non-standard way of defining dimensions

        dfn_values=[]

        for E,V_pitch,R,Z in zip(kwargs['E'],kwargs['V_pitch'],kwargs['R'],kwargs['Z']):

            diff_R=np.abs(dfn['R2D']-R)**2 #irregular grid so need to minimise distance to every point
            diff_Z=np.abs(dfn['Z2D']-Z)**2

            index_RZ=(diff_R+diff_Z).argmin()
            index_V_pitch=np.abs(dfn['V_pitch']-V_pitch).argmin()
            index_E=np.abs(dfn['E']-E).argmin()

            dfn_values.append(dfn['dfn'][index_RZ,index_V_pitch,index_E])

            print('R') #diagnostic printing out to show the point is in right ball park
            print(dfn['R2D'][index_RZ])
            print('Z') #diagnostic printing out to show the point is in right ball park
            print(dfn['Z2D'][index_RZ])
            print('V_pitch') #diagnostic printing out to show the point is in right ball park
            print(dfn['V_pitch'][index_V_pitch])  
            print('E') #diagnostic printing out to show the point is in right ball park
            print(dfn['E'][index_E])

    else:

        number_dimensions=len(kwargs.keys())
        for key in kwargs:
            number_points=len(kwargs[key])

        coordinate_array=np.ndarray(shape=(number_points,number_dimensions),dtype=int) #hold nD indices of all the N points we want to know in an NxnD array for an nD distribution function

        for key in kwargs: #go through supplied dimensions

            for dimension in dfn['dfn_index']: #find the dfn index of the dimension we're dealing with
                if key == dimension:
                    dimension_index=dfn['dfn_index'].tolist().index(key)

            for counter,point in enumerate(kwargs[key]):
                point_index=np.abs(dfn[key]-point).argmin() #find index of nearest grid point to the point we're dealing with along this dimension 
                coordinate_array[counter,dimension_index]=point_index 
                #print(key) #diagnostic printing out to show the point is in right ball park
                #print(dfn[key][point_index])

        dfn_values=[]
        for coordinate in coordinate_array:
            dfn_values.append(dfn['dfn'][tuple(coordinate)])


    dfn_values=np.asarray(dfn_values)
    return dfn_values