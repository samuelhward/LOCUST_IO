#perturbation.py
 
"""
Samuel Ward
29/07/2018
----
class to handle LOCUST perturbation input data
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
try:
    import imas 
except:
    print("WARNING: IMAS module could not be imported!\nreturning\n")
try:
    from processing import utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    from classes import base_input 
except:
    raise ImportError("ERROR: base_input.py could not be imported!\nreturning\n")
    sys.exit(1) 
try:
    from classes import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
 
pi=np.pi


################################################################## Perturbation functions
 
def read_perturbation_MARSF(filepath):
    """
    reads perturbation stored in LOCUST format

    notes:
        assumes R is slower-varying dimension in file when inferring dimensions
    """

    print("reading MARSF perturbation")

    with open(filepath,'r') as file:
                 
        #initialise data
        input_data={}
        input_data['R_2D']=[]
        input_data['Z_2D']=[]
        input_data['B_field_R_real']=[]
        input_data['B_field_R_imag']=[]
        input_data['B_field_Z_real']=[]
        input_data['B_field_Z_imag']=[]
        input_data['B_field_tor_real']=[]
        input_data['B_field_tor_imag']=[]

        #read lazily
        for line in file:
            split_line=line.split()
            input_data['R_2D'].append(float(split_line[0]))
            input_data['Z_2D'].append(float(split_line[1]))
            input_data['B_field_R_real'].append(float(split_line[2]))
            input_data['B_field_R_imag'].append(float(split_line[3]))
            input_data['B_field_Z_real'].append(float(split_line[4]))
            input_data['B_field_Z_imag'].append(float(split_line[5]))
            input_data['B_field_tor_real'].append(float(split_line[6]))
            input_data['B_field_tor_imag'].append(float(split_line[7]))

        input_data['R_2D']=np.asarray(input_data['R_2D'])
        input_data['Z_2D']=np.asarray(input_data['Z_2D'])
        input_data['B_field_R_real']=np.asarray(input_data['B_field_R_real'])
        input_data['B_field_R_imag']=np.asarray(input_data['B_field_R_imag'])
        input_data['B_field_Z_real']=np.asarray(input_data['B_field_Z_real'])
        input_data['B_field_Z_imag']=np.asarray(input_data['B_field_Z_imag'])
        input_data['B_field_tor_real']=np.asarray(input_data['B_field_tor_real'])
        input_data['B_field_tor_imag']=np.asarray(input_data['B_field_tor_imag'])

        
        #infer the grid dimensions
        Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
        R_dim=int((input_data['R_2D'].size)/Z_dim)

        #reshape to grid dimensions
        input_data['R_2D']=input_data['R_2D'].reshape(R_dim,Z_dim)
        input_data['Z_2D']=input_data['Z_2D'].reshape(R_dim,Z_dim)
        input_data['B_field_R_real']=input_data['B_field_R_real'].reshape(R_dim,Z_dim)
        input_data['B_field_R_imag']=input_data['B_field_R_imag'].reshape(R_dim,Z_dim)
        input_data['B_field_Z_real']=input_data['B_field_Z_real'].reshape(R_dim,Z_dim)
        input_data['B_field_Z_imag']=input_data['B_field_Z_imag'].reshape(R_dim,Z_dim)
        input_data['B_field_tor_real']=input_data['B_field_tor_real'].reshape(R_dim,Z_dim)
        input_data['B_field_tor_imag']=input_data['B_field_tor_imag'].reshape(R_dim,Z_dim)
        

    print("finished reading MARSF perturbation")
    
    return input_data
 
def dump_perturbation_MARSF(output_data,filepath):
    """
    writes perturbation to MARSF format

    notes:
    """
 
    print("writing MARSF perturbation")

    with open(filepath,'w') as file: #open file
        pass
     
    print("finished writing MARSF perturbation")

 
################################################################## perturbation class
 
class Perturbation(base_input.LOCUST_input):
    """
    class describing magnetic field perturbation for LOCUST
 
    inheritedfrom LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'perturbation'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in input_files folder
        self.filepath               full path of file in input_files folder  
        self.shot                   shot number
        self.run                    run number
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in input_files folder
 
    notes:
    """
 
    LOCUST_input_type='perturbation'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read perturbation from file 
 
        notes:
        """
 
        if utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='MARSF': #here are the blocks for various file types, they all follow the same pattern
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_perturbation_MARSF(self.filepath) #read the file
         
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (MARSF)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write perturbation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)

        if utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='MARSF':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to MARSF - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_perturbation_MARSF(self.data,filepath)

        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (MARSF)\n")

 
#################################
 
##################################################################
 
###################################################################################################