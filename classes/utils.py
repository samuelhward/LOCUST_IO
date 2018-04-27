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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
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
    if all(arg is True for arg in args):
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
        if key is not None and value is not None: 
            data[key]=value

def safe_set(target,source):
    """
    generalised function to set a target value to source if it exists
    notes:
    """
    if source is not None:
        target=source

 
def fortran_string(number_out,length,decimals=None):
    """
    produces a string stream of specified length, padded with space at the start

    notes:
        supposed to be a quick fix to fortran format descriptors
        decimals (required for non-ints) controls how many decimal places the number is written to (assumes exponential format)
    """

    if decimals is not None:
        output_string='{:.{decimals}e}'.format(number_out,decimals=decimals)        
    else: #assume integer if decimals not specified
        output_string=str(number_out)

    number_spaces=length-len(output_string) #caclulate amount of padding needed from length of the number_out
    if number_spaces>=0:
        string_stream='{}{}'.format(' '*number_spaces,output_string) #construct the string stream
        return string_stream
    else:
        print('WARNING: fortran_string() must not take a number with more digits than desired length of output stream')


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

