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
    import os
    import re
    import time
    import itertools
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from scipy.io import netcdf as ncdf
    import scipy.interpolate
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import h5py
except:
    print("WARNING: h5py could not be imported!\n") 
try:
    from classes import support
except:
    raise ImportError("ERROR: LOCUST_IO/classes/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from processing import plot_output
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_output.py could not be imported!\nreturning\n")
try:
    from processing import process_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")


np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
plot_style_LCFS='m-' #set plot style for LCFS
plot_style_limiters='w-' #set plot style for limiters

pi=np.pi
e_charge=1.60217662e-19 #define electron charge
mass_neutron=1.674929e-27 #define mass of neutron
amu=1.66053904e-27
mass_deuterium=2.0141017781*amu



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
         
def file_numbers(ingf):#instance of generators are objects hence can use .next() method on them
    """
    generator to read numbers in a file
 
    notes:
        originally written by Ben Dudson
 
        calling get_next() on a generator produced with this function will permanently advance the iterator in this generator
        when generators reach the end, that's permanently it!
    """
    toklist = []
    while True: #runs forever until break statement (until what is read is not a line) below will cause a return - which permanently stops a generator
        line = ingf.readline() #read in line from INGF
        if not line: break #break if readline() doesnt return a line i.e. end of file
        line = line.replace("NaN","-0.00000e0") #replaces NaNs with 0s
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?' #regular expression to find numbers
        toklist = re.findall(pattern,line) #toklist now holds all the numbers in a file line (is essentially a single file line)
        for tok in toklist: #yield every number in that one line individually
            yield tok #so in terms of iteration, using get_next() will cause another line to be read
 
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

def interpolate_2D(X_axis,Y_axis,Z_grid,type='RBS',rect_grid=True):
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
        rect_grid - toggle whether X_axis, Y_axis are regular rectangular grid edges or arbitrary cordinates 

    usage:
        my_interpolator=interpolate_2D(X_axis,Y_axis,data)
        interpolated_value=my_interpolator(x,y)
    """
    
    if type=='RBF':
        if rect_grid:
            Y_grid,X_grid=np.meshgrid(Y_axis,X_axis) #swap since things are defined r,z 
        else:
            Y_grid,X_grid=Y_axis,X_axis
        interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function='multiquadric',smooth=0)
    
    elif type=='RBS':
        interpolator=scipy.interpolate.RectBivariateSpline(X_axis,Y_axis,Z_grid) #normally order is other way in RBS but I have swapped my axes

    return interpolator

def interpolate_1D(X_axis,Y_axis):
    """
    generate a 1D line interpolator

    notes:
        keep as separate functions so can freely swap out interpolation method
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    """

    interpolator=scipy.interpolate.Rbf(X_axis,Y_axis,function='cubic',smooth=0)

    return interpolator

def RphiZ_to_XYZ(R,phi):
    """
    converts R,phi positions to X,Y 
    """

    X=R*np.sin(phi)
    Y=R*np.cos(phi)

    return X,Y

def XYZ_to_RphiZ(X,Y):
    """
    converts X,Y positions to R,phi 
    """

    phi=np.arctan2(Y,X)
    R=X*np.cos(phi)+Y*np.sin(phi)

    return R,phi

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
        eq.set(fpolrz=process_input.fpolrz_calc(eq)) 
    if 'B_field' not in some_equilibrium.data.keys():
        eq.set(B_field=process_input.B_calc(eq))

    B_field_R_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,0]) #generate interpolators for B
    B_field_tor_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,1])
    B_field_Z_interpolator=interpolate_2D(eq['R_1D'],eq['Z_1D'],eq['B_field'][:,:,2])

    Br_at_particles=B_field_R_interpolator(some_particle_list['R'],some_particle_list['Z']) #interpolate B to particles
    Btor_at_particles=B_field_tor_interpolator(some_particle_list['R'],some_particle_list['Z'])
    Bz_at_particles=B_field_Z_interpolator(some_particle_list['R'],some_particle_list['Z'])
    B_at_particles=np.sqrt(Br_at_particles**2+Btor_at_particles**2+Bz_at_particles**2) #calculate |B| at particles - quicker than interpolating

    Br_at_particles/=B_at_particles #calculate B_hat
    Btor_at_particles/=B_at_particles
    Bz_at_particles/=B_at_particles

    V_parallel=dot_product([some_particle_list['V_R'],some_particle_list['V_tor'],some_particle_list['V_Z']],[Br_at_particles,Btor_at_particles,Bz_at_particles])

    V=np.sqrt(some_particle_list['V_R']**2+some_particle_list['V_tor']**2+some_particle_list['V_Z']**2)
    V_pitch=V_parallel/V

    return V_pitch