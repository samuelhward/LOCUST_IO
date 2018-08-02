#utils.py
 
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
    import re
    import time
    import itertools
    import scipy.interpolate
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
 
pi=np.pi



################################################################## Supporting functions
 
def none_check(ID,LOCUST_input_type,error_message,*args):
    """
    generic function for checking if a None value appears in *args
 
    notes:
        message should be something specific to the section of code which called none_check
    """
    if all(arg is not False for arg in args):
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

def interpolate_2D(X_axis,Y_axis,Z_grid,type='RBS'):
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

    usage:
        my_interpolator=interpolate_2D(X_axis,Y_axis,data)
        interpolated_value=my_interpolator(x,y)
    """
    
    if type=='RBF':
        Y_grid,X_grid=np.meshgrid(Y_axis,X_axis) #swap since things are defined r,z 
        interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function='cubic',smooth=0)
    
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


def dump_profiles_ASCOT(filename,temperature_i,temperature_e,density_i,density_e,rotation_toroidal):
    """
    dumps collection of kinetic profiles to ASCOT input.plasma_1d format

    notes:
        profiles must be mapped to same poloidal flux axis
        currently allows single ion species - in future, pass list of profile objects and iterate
        if rotation profile added to LOCUST data in future, will need to update this function accordingly
    args:
        filename - output filename
        temperature_e - electron temperature object (eV)
        temperature_i - ion temperature object (eV)
        density_e - electron density object (#/m^3)
        density_i - ion density object (#/m^3)
        rotation_toroidal - array-like toroidal rotation (rad/s)
    """

    print("dumping profiles to ASCOT format")
 
    filepath=support.dir_input_files+filename
    with open(filepath,'w') as file:

        file.write("# Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation\n")
        file.write("# range must cover [0,1] of normalised poloidal rho. It can exceed 1. \n")
        file.write("# 18Jan08 for testing (first 3 lines are comment lines)\n")
        file.write(str(temperature_e['flux_pol_norm'].size)+"   "+"1"+" # Nrad,Nion\n")
        file.write("1           # ion Znum\n")
        file.write("1           # ion Amass\n")
        file.write("1 1         # collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons\n")
        file.write("    RHO (pol)       Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)\n")

        for RHO,Te,Ne,Vtor_I,Ti1,Ni1 in zip(temperature_e['flux_pol_norm'],temperature_e['T'],density_e['n'],rotation_toroidal,temperature_i['T'],density_i['n']): 
            line=fortran_string(RHO,16,7)+fortran_string(Te,16,7)+fortran_string(Ne,16,7)+fortran_string(Vtor_I,15,7)+fortran_string(Ti1,17,7)+fortran_string(Ni1,15,7)+"\n"
            file.write(line)

    print("finished dumping profiles to ASCOT format")

def dump_wall_ASCOT(filename,output_data):
    """
    dumps 2D wall outline to ASCOT input.wall_2d format
    
    notes:
        currently sets all divertor flags = 0 i.e. treats everything as 'wall'
    args:
        filename - output filename
        output_data - data structure holding wall with same wall variable names as GEQDSK equilibrium i.e. rlim,zlim

    """

    print("dumping wall to ASCOT format")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file:

        file.write("{number_points} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n".format(number_points=int(output_data['rlim'].size)))
        
        for r,z in zip(output_data['rlim'],output_data['zlim']):
            line=fortran_string(r,16,7)+fortran_string(z,16,7)+fortran_string(0.0,4,0,False)+"\n"
            file.write(line)

    print("finished dumping wall to ASCOT format")

def dump_rotation_MARSF(filename,output_data):
    """
    writes rotation profile to MARSF Mogui ASCII format 

    notes
        assumes structure similar number_density or temperature, with normalised poloidal flux given against rotation
        re-usable if rotation class written for LOCUST_IO
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    args:
        filename - output filename
        output_data - data structure holding 'rotation' against normalised poloidal flux
    """

    print("writing rotation to MARSF mogui")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file: #open file

        normalised_flux_sqrt=np.sqrt(np.abs(output_data['flux_pol_norm'])) #take abs sqrt()
        normalised_flux_sqrt,output_data['rotation']=utils.sort_arrays(normalised_flux_sqrt,output_data['n']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(normalised_flux_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(normalised_flux_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm_sqrt} {rotation}\n".format(flux_pol_norm_sqrt=utils.fortran_string(normalised_flux_sqrt[point],24,18),rotation=utils.fortran_string(output_data['rotation'][point],24,18)))

    print("finished writing rotation to MARSF mogui")


def dump_coil_currents_MARSF(filename,output_data):
    """
    writes coil currents to MARSF Mogui ASCII format

    notes:
        output_data structure assumed to be dict holding 'upper_coil_currents' and 'lower_coil_currents' arrays
        MARSF mogui written by David Ryan
    args:  
        filename - output filename
        output_data - data structure holding coil currents
    """

    print("writing coil currents to MARSF mogui")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file: #open file

        some_ID='34835_2000'
        file.write('{}\n'.format(some_ID)) #write headerline
        
        upper_coil_currents=''
        lower_coil_currents=''

        for coil in output_data['upper_coil_currents']:
            upper_coil_currents+=str(coil)+', '

        for coil in output_data['lower_coil_currents']:
            lower_coil_currents+=str(coil)+', '

        upper_coil_currents=upper_coil_currents[0:-2]
        lower_coil_currents=lower_coil_currents[0:-2]

        upper_coil_currents+='\n'
        lower_coil_currents+='\n'

        file.write(upper_coil_currents) #write upper coil currents
        file.write(lower_coil_currents) #write lower coil currents

    print("finished writing coil currents to MARSF mogui")






