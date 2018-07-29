#hdf5.py

"""
Samuel Ward
27/04/2018
----
class to handle LOCUST HDF5 output data
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
    import h5py
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
    from classes import base_output 
except:
    raise ImportError("ERROR: base_output.py could not be imported!\nreturning\n")
    sys.exit(1) 
try:
    from classes import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)


np.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays

################################################################## Orbits functions

def read_HDF5_LOCUST(filepath):
    """
    reads generic HDF5 data output from LOCUST

    notes:

    """

    print("reading HDF5 from LOCUST")

    input_data={} #initialise data dictionary
    file = h5py.File(filepath, 'r')
    
    input_data['sqrt(PSIn)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['sqrt(PSIn)']) 
    input_data['density']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Density']['data'])
    input_data['density_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Density']['error'])
    input_data['energy']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Energy']['data'])
    input_data['energy_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Energy']['error'])
    input_data['J(NBCD)-raw']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['J(NBCD)-raw']['data'])
    input_data['J(NBCD)-raw_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['J(NBCD)-raw']['error'])
    input_data['NBI-heating-power(TOT)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(TOT)']['data'])
    input_data['NBI-heating-power(TOT)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(TOT)']['error'])
    input_data['NBI-heating-power(e-)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(e-)']['data'])
    input_data['NBI-heating-power(e-)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(e-)']['error'])
    input_data['NBI-heating-power(i1)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i1)']['data'])
    input_data['NBI-heating-power(i1)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i1)']['error'])
    input_data['NBI-heating-power(i1)A']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i1)']['A'])
    input_data['NBI-heating-power(i1)Z']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i1)']['Z'])
    input_data['NBI-heating-power(i2)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i2)']['data'])
    input_data['NBI-heating-power(i2)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i2)']['error'])
    input_data['NBI-heating-power(i2)Z']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i2)']['Z'])
    input_data['NBI-heating-power(i2)A']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['NBI-heating-power(i2)']['A'])
    input_data['P(para)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['P(para)']['data'])
    input_data['P(para)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['P(para)']['error'])
    input_data['P(perp)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['P(perp)']['data'])
    input_data['P(perp)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['P(perp)']['error'])
    input_data['residual-angular-momentum-density']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Residual-angular-momentum-density']['data'])
    input_data['residual-angular-momentum-density_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Residual-angular-momentum-density']['error'])
    input_data['torque-density(JxB-inst)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(JxB-inst)']['data'])
    input_data['torque-density(JxB-inst)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(JxB-inst)']['error'])
    input_data['torque-density(JxB-sweep)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(JxB-sweep)']['data'])
    input_data['torque-density(JxB-sweep)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(JxB-sweep)']['error'])
    input_data['torque-density(coll)']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(coll)']['data'])
    input_data['torque-density(coll)_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['Torque-density(coll)']['error'])
    input_data['e-source']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['e-Source']['data'])
    input_data['e-source_err']=np.array(file['Output Data']['Fast Iones']['Profiles (1D)']['e-Source']['error'])
    
    print("finished reading HDF5 from LOCUST")

    return input_data



################################################################## Orbits class

class HDF5(base_output.LOCUST_output):
    """
    class describing a generic HDF5 output for LOCUST
    
    inheritedfrom LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'HDF5'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in output_files folder
        self.filepath               full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in output_files folder

    notes:
        data is stored such that...
    """

    LOCUST_output_type='HDF5'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None):
        """
        read HDF5 from file 

        notes:
        """

        if utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties=properties
                self.data=read_HDF5_LOCUST(self.filepath) #read the file
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write HDF5 to file

        notes: 
        """
        if utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_HDF5_LOCUST(self.data,filepath)
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST)\n")



#################################

##################################################################

###################################################################################################