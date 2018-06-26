#base_input.py
 
"""
Samuel Ward
02/11/2017
----
base class for LOCUST input classes
---
 
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
    from processing import utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    from classes import support
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
 
pi=np.pi


################################################################## Base class
 
class LOCUST_input:
    """
    base class for a generic LOCUST input object
 
    self.ID                     unique object identifier, good convention to fill these for error handling etc
    self.data                   holds all input data in dictionary object
    self.LOCUST_input_type      string which holds this class' input type
 
    notes:
    """
 
    LOCUST_input_type='base_input'
 
    def __init__(self,ID,data_format=None,filename=None,shot=None,run=None,**properties): #this is common to all children (not overloaded), must have ID
 
        self.ID=ID #always set the ID, even if we don't invoke read_data i.e. a blank object is initialised
        self.data={} 
        if not utils.none_check(self.ID,self.LOCUST_input_type,"read_data requires data_format, blank input initialised \n",data_format):
            self.read_data(data_format,filename,shot,run,**properties)
 
    def __getitem__(self,key):
        """
        access member data via []        
        """

        return self.data[key]
 
    def __setitem__(self,key,value):
        """
        access member data via []
        """

        self.data[key]=value

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties): #bad practice to change overridden method signatures, so retain all method arguments             
        """
        read data to be overloaded in all children classes
        """

        pass
 
    def look(self):
        """
        print class information and data
        """

        print("\n-----------------------")
        print("ID - {ID}".format(ID=self.ID))  
        print("Output Type - {LOCUST_input_type}".format(LOCUST_input_type=self.LOCUST_input_type))

        if hasattr(self,'data_format'):
            print("Data Format - {data_format}".format(data_format=self.data_format))
        
        if hasattr(self,'filename'):
            print("Input Filename - {filename}".format(filename=self.filename))
        
        if hasattr(self,'shot'):
            print("Shot - {shot}".format(shot=self.shot))
        
        if hasattr(self,'run'):
            print("Run - {run}".format(run=self.run))

        if hasattr(self,'properties') and self.properties:
            print("Properties:".format(properties=self.properties))
            for key in self.properties:
                if any(self.properties[key]): #do not print if the data is empty
                    print("{key} - {value}".format(key=key,value=self.properties[key])) 
        
        print("|")
        print("|")

        if hasattr(self,'data') and self.data:
            for key in self.data:
                if type(self.data[key])==type({}): #check for dicts since these mess things up
                    if self.data[key]: #if dict then check if dict is empty
                        print(key+":")
                        for sub_key in self.data[key]: 
                            print("     {sub_key} - {value}".format(sub_key=sub_key,value=self.data[key][sub_key]))
                elif not self.data[key].size==0: #if not a dict assume a numpy array
                    print("{key} - {value}".format(key=key,value=self.data[key]))
                
        print("-----------------------\n")
 
    def copy(self,target,*keys):
        """
        copy two input objects
 
        notes:
            if target.data is None or contains Nones then this function does nothing
            if no key supplied then copy all data over
            if key supplied then copy/append dictionary data accordingly
  
        usage:
            my_input.copy(some_other_input) to copy all data
            my_input.copy(some_other_input,'some_key','some_other_key') to copy specific fields
            my_input.copy(some_other_input, *some_list_of_args) equally
        """
        if utils.none_check(self.ID,self.LOCUST_input_type,"cannot copy() - target.data is blank\n",target.data): #return warning if any target data contains empty variables
            pass

        elif not keys: #if empty, keys will be false i.e. no key supplied --> copy everything 
            self.data=copy.deepcopy(target.data) #using = with whole dictionary results in copying by reference, so need deepcopy() here
            if hasattr(target,'properties'): #copy properties field
                self.properties=target.properties

        elif not utils.none_check(self.ID,self.LOCUST_input_type,"cannot copy() - found key containing None\n",*keys): 
            self.set(**{key:target[key] for key in keys}) #call set function and generate the dictionary of **kwargs with list comprehension
            if hasattr(target,'properties'): #copy properties field
                self.properties=target.properties
     
    def set(self,**kwargs):
        """
        set input object data 
 
        notes:
            specific to LOCUST_IO classes due to none_check call - see safe_set() for more general function

        usage:
            my_input.set(some_key=5,some_other_key=[1,2,3,4]) to set multiple values simultaneously
            my_input.set(**{'some_key':100,'some_other_key':200}) equally
        """

        keys=list(kwargs.keys())
        values=list(kwargs.values())
        allkeysvalues=keys+values

        if utils.none_check(self.ID,self.LOCUST_input_type,"cannot set() - empty key/value pair found\n",*allkeysvalues):
            pass
        else:
            for key,value in zip(keys,values): #loop through kwargs
                self[key]=value
 
    def compare(self,target,verbose=False):
        """
        compare two input objects

        notes:
            returns true if all data held by target is also held by self (self can have excess)
            verbose option prints summary of compare results
        """

        data_missing_self=[]
        data_missing_target=[]
        data_different=[]

        for key in target.data: #record fields we do not have that target does
            if not key in self.data:
                data_missing_self.append(key) 
            
            elif self[key].size!=target[key].size: #both contain data but data is different
                data_different.append(key)
            elif not np.allclose(self[key],target[key]): #must check size first before doing np.allclose()
                data_different.append(key)

        for key in self.data: #record fields we have that target does not
            if not key in target.data:
                data_missing_target.append(key) 

        if verbose is True: #if wanting to print summary
            if data_missing_self:
                print("self is missing:")
                print('\n'.join(str(key) for key in data_missing_self)) 
            if data_missing_target:
                print("target is missing:")
                print('\n'.join(str(key) for key in data_missing_target)) 
            if data_different: 
                print("shared different data:")
                print('\n'.join(str(key) for key in data_different))

        if not data_different and not data_missing_self: #if shared data is the same and self all target data 
            return True #self is same as target
        else:
            return False

    def run_check(self,verbose=False):
        """
        checks whether object is ready to run in LOCUST

        notes:
            only checks if correct data is held by object
        """

        missing_data=[]

        if not self.data:
            if verbose is True:
                print("no data")
            return False

        if self.LOCUST_input_type=='equilibrium':
            for key in support.required_equilibrium:
                if key not in self.data:
                    missing_data.append(key)

        elif self.LOCUST_input_type=='beam_deposition':
            for key in support.required_beam_deposition:            
                if key not in self.data:
                    missing_data.append(key)

        elif self.LOCUST_input_type=='temperature':
            for key in support.required_temperature:
                if key not in self.data:
                    missing_data.append(key)

        elif self.LOCUST_input_type=='number_density':
            for key in support.required_number_density:
                if key not in self.data:
                    missing_data.append(key)
 
        if missing_data: #object is not ready
            if verbose is True: #if wanting to print summary
                print("missing data:")
                print('\n'.join(str(key) for key in missing_data))
                print('\n')
            return False
        else:
            return True 



#################################
 
##################################################################
 
###################################################################################################