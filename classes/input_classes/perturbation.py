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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_input 
except:
    raise ImportError("ERROR: LOCUST_IO/classes/base_input.py could not be imported!\nreturning\n")
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

################################################################## Perturbation read functions
 
def read_perturbation_LOCUST(filepath):
    """
    reads perturbation stored in LOCUST format

    notes:
        assumes R is slower-varying dimension in file when inferring dimensions
    """

    print("reading LOCUST perturbation")

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
        

    print("finished reading LOCUST perturbation")
    
    return input_data

def read_perturbation_IDS(shot,run):
    """
    notes:

    """ 

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_perturbation_IDS could not import IMAS module!\nreturning\n")
        return

def read_perturbation_JOREK(filepath):
    """
    notes:
    """

    try:
        import h5py
    except:
        print("ERROR: read_perturbation_JOREK could not import h5py!\n")
        return

    input_data={} #initialise data dictionary
    file = h5py.File(filepath, 'r')

    input_data['']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['sqrt(PSIn)']) 

    return input_data
    
################################################################## Perturbation write functions
 
def dump_perturbation_LOCUST(output_data,filepath):
    """
    writes perturbation to LOCUST format

    notes:
    """
 
    print("writing LOCUST perturbation")

    with open(filepath,'w') as file: #open file
        pass
     
    print("finished writing LOCUST perturbation")

################################################################## perturbation class
 
class Perturbation(classes.base_input.LOCUST_input):
    """
    class describing magnetic field perturbation for LOCUST
 
    inherited from LOCUST_input:
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
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_perturbation_LOCUST(self.filepath) #read the file

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from magnetics IDS - shot and run required\n",shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_perturbation_IDS(self.shot,self.run)

        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (MARSF/IDS)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write perturbation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)

        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_perturbation_LOCUST(self.data,filepath)

        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (MARSF)\n")

    def plot(self,some_equilibrium=False,key='Br',axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,number_bins=20,fill=True,ax=False,fig=False):
        """
        notes:
        args:
            some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
            key - select which data to plot
            axes - define plot axes in x,y order or as full list of indices/slices (see dfn_transform())
            LCFS - show plasma boundary outline (requires equilibrium arguement)
            limiters - toggles limiters on/off in 2D plots
            real_scale - plot to Tokamak scale
            colmap - set the colour map (use get_cmap names)
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """
        
        import scipy
        import numpy as np
        import matplotlib
        from matplotlib import cm
        import matplotlib.pyplot as plt
        from mpl_toolkits import mplot3d #import 3D plotting axes
        from mpl_toolkits.mplot3d import Axes3D

        pass
#################################
 
##################################################################
 
###################################################################################################