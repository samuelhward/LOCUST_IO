#output_classes.py

"""
Samuel Ward
15/01/2018
----
python module for classes which hold LOCUST's output data
contains methods for reading/writing/converting/plotting/manipulating LOCUST output data
---
usage:
    see README.md for usage

notes:         
---
"""


###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from output_classes import x") but best practice to import whole output_classes module anyway
try:
    import numpy as np
    import copy
    import re
    import time
    import itertools
    import copy
except:
    raise ImportError("ERROR: initial imported modules not found!\nreturning\n")
    sys.exit(1)
try:
    import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py not found in this directory!\nreturning\n") 
    sys.exit(1)
try:
    import imas 
except:
    raise ImportError("ERROR: IMAS module not found!\nreturning\n")
    sys.exit(1)





















###################################################################################################
#Main Code









################################################################## Supporting functions

def none_check(ID,LOCUST_output_type,error_message,*args):
    """
    generic function for checking if a None value appears in *args

    notes:
        message should be something specific to the section of code which called none_check
    """
    if all(arg is not None for arg in args):
        return False
    else:
        print("WARNING: none_check returned True (LOCUST_output_type={LOCUST_output_type}, ID={ID}): {message}".format(LOCUST_output_type=LOCUST_output_type,ID=ID,message=error_message))
        return True
        














################################################################## Base class

class LOCUST_output:
    """
    base class for a generic LOCUST output object

    self.ID                     unique object identifier, good convention to fill these for error handling etc
    self.data                   holds all input data in dictionary object
    self.LOCUST_output_type     string which holds this class' output type

    notes:
    """

    LOCUST_output_type='base_output'

    def __init__(self,ID,data_format=None,input_filename=None,shot=None,run=None,properties=None): #this is common to all children (not overloaded), must have ID

        self.ID=ID #always set the ID, even if we don't invoke read_data i.e. a blank object is initialised
        if not none_check(self.ID,self.LOCUST_output_type,"read_data requires data_format, blank output initialised \n",data_format):
            self.read_data(data_format,input_filename,shot,run,properties)

    def read_data(self,data_format=None,input_filename=None,shot=None,run=None,properties=None): #bad practice to change overridden method signatures, so retain all method arguments             
        self.data=None #read_data definition is designed to be overloaded in children classes 




























################################################################## Orbits functions

def read_orbits_ASCII(input_filepath):
    """
    reads orbits stored in ASCII format -   r z phi

    notes:
    reads in a headerline for number of particles
    reads in a footerline for number of time steps
    """

    with open(input_filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_orbits_ASCII() cannot read from "+input_filepath)
    
    number_particles=int(lines[0]) #extract number of particles
    number_timesteps=int(lines[-1])-1 #extract number of time steps of each trajectory

    del(lines[0])
    del(lines[-1])

    input_data = {} #initialise the dictionary to hold the data
    input_data['orbits']=np.array([[[float(coord) for coord in time_slice.split()] for time_slice in lines[particle*number_timesteps:(particle+1)*number_timesteps]] for particle in range(number_particles)],ndmin=3)  #read in the data in one line using list comprehension again! (see dump_orbits_ASCII for a more digestable version of this)
    #XXX still need to verify if the data is nested in this order
    input_data['number_particles']=np.asarray(number_particles)
    input_data['number_timesteps']=np.asarray(number_timesteps)
   
    return input_data

def dump_orbits_ASCII(output_data,output_filepath): 
    """
    writes orbits to ASCII format -  r z phi
    
    notes:
        writes out a headerline for number of particles
        writes out a footerline for number of time steps
    """

    with open(output_filepath,'w') as file: #open file

        file.write("{}\n".format(len(output_data['number_particles']))) #re-insert line containing number of particles

        for particle in range(len(output_data['orbits'][:][0][0])): #loop over everything again
            for time_slice in range(len(output_data['orbits'][0][:][0])):

                    file.write("{r} {z} {phi}\n".format(r=output_data['orbits'][particle][time_slice][0],z=output_data['orbits'][particle][time_slice][1],phi=output_data['orbits'][particle][time_slice][2]))

        file.write("{}".format(len(output_data['number_timesteps']))) #re-insert line containing number of time steps

################################################################## Orbits class

class Orbits(LOCUST_output):
    """
    class describing orbits output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'orbits'
    class data
        self.data_format            data format of original data e.g. ASCII
        self.input_filename         name of file in output_files folder
        self.input_filepath         full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        output_filename             name of file to write to
        output_filepath             full path to output file in output_files folder

    notes:
        data is stored such that coordinate i at time t for particle p is my_orbit['orbits'][p,t,i]
        in this way, a single particle's trajectory is my_orbit['orbits'][p,:,i]
    """

    LOCUST_output_type='orbits'

    def __getitem__(self,key):

        return self.data[key]

    def __setitem__(self,key,value):

        self.data[key]=value

    def read_data(self,data_format=None,input_filename=None,shot=None,run=None,properties=None):
        """
        read orbits from file 

        notes:
        """

        if none_check(self.ID,self.LOCUST_output_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='ASCII': #here are the blocks for various file types, they all follow the same pattern
            if not none_check(self.ID,self.LOCUST_output_type,"cannot read_data from ASCII - input_filename required\n",input_filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.input_filename=input_filename
                self.input_filepath=support.dir_output_files+input_filename
                self.properties=properties
                self.data=read_orbits_ASCII(self.input_filepath) #read the file
        else:
            print("cannot read_data - please specify a compatible data_format (ASCII)\n")            

    def dump_data(self,data_format=None,output_filename=None,shot=None,run=None):
        """
        write orbits to file

        notes: 
        """
        if none_check(self.ID,self.LOCUST_output_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='ASCII':
            if not none_check(self.ID,self.LOCUST_output_type,"cannot dump_data to ASCII - output_filename required\n",output_filename):
                output_filepath=support.dir_output_files+output_filename
                dump_orbits_ASCII(self.data,output_filepath)
        else:
            print("cannot dump_data - please specify a compatible data_format (ASCII)\n")

    def copy(self,target,*keys):
        """
        copy two orbit objects 

        notes:
            if target.data is None or contains Nones then this function does nothing
            if no key supplied then copy all data over
            if key supplied then copy/append dictionary data accordingly
                        
        usage:
            my_orbits.copy(some_other_orbits) to copy all data
            my_orbits.copy(some_other_orbits,'some_arg','some_other_arg') to copy specific fields
            my_orbits.copy(some_other_orbits, *some_list_of_args) equally
        """
        if none_check(self.ID,self.LOCUST_output_type,"cannot copy() - target.data is blank\n",target.data): #return warning if any target data contains empty variables
            pass
        elif not keys: #if empty, keys will be false i.e. no key supplied --> copy everything 
            self.data=copy.deepcopy(target.data) #using = with whole dictionary results in copying by reference, so need deepcopy() here
        elif not none_check(self.ID,self.LOCUST_output_type,"cannot copy() - found key containing None\n",*keys): 
            self.set(**{key:target[key] for key in keys}) #call set function and generate the dictionary of **kwargs with list comprehension
    
    def set(self,**kwargs):
        """
        set orbit object data 

        usage:
            my_orbits.set(some_arg=5,some_other_arg=[1,2,3,4]) to set multiple values simultaneously
            my_orbits.set(**{'some_arg':100,'some_other_arg':200}) equally
        """
        keys=kwargs.keys()
        values=kwargs.values()
        allkeysvalues=keys+values #NOTE can avoid having to do this in python version 3.5
        if none_check(self.ID,self.LOCUST_output_type,"cannot set() - empty key/value pair found\n",*allkeysvalues):
            pass
        else:
            for key,value in zip(keys,values): #loop through kwargs
                self[key]=value








#################################

##################################################################

###################################################################################################
