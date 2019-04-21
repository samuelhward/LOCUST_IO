#rotation.py
 
"""
Samuel Ward
21/04/2019
----
class to handle LOCUST rotation input data
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
    import pathlib
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


################################################################## rotation read functions
 
def read_rotation_LOCUST(filepath,**properties):
    """
    reads rotation profile stored in LOCUST format
    """

    pass #XXX not yet implemented

################################################################## rotation write functions

def dump_rotation_LOCUST(output_data,filepath,**properties):
    """
    writes rotation profile stored in LOCUST format
    """

    pass #XXX not yet implemented

def dump_rotation_MARSF(output_data,filepath,**properties):
    """
    writes rotation profile to MARSF Mogui ASCII format 

    notes
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    """

    print("writing rotation to MARSF mogui")

    filepath=support.dir_input_files / filename

    with open(filepath,'w') as file: #open file

        flux_pol_norm_sqrt=np.sqrt(np.abs(some_rotation['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,rotation=processing.utils.sort_arrays(flux_pol_norm_sqrt,some_rotation['rotation']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(flux_pol_norm_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(flux_pol_norm_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm_sqrt}{rotation}\n".format(flux_pol_norm_sqrt=processing.utils.fortran_string(flux_pol_norm_sqrt[point],24,18),rotation=processing.utils.fortran_string(rotation[point],24,18)))

    print("finished writing rotation to MARSF mogui")
        
################################################################## rotation class
 
class rotation(classes.base_input.LOCUST_input):
    """
    class describing rotation profile input for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'rotation'
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
 
    LOCUST_input_type='rotation'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,time=None,**properties):
        """
        read rotation from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_rotation_LOCUST(self.filepath,**properties)
    
        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write rotation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID={})".format(self.ID))
            
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_rotation_LOCUST(self.data,filepath)
 
        elif data_format=='MARSF':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to MARSF - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_rotation_MARSF(self.data,filepath)

        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/MARSF)\n".format(self.ID))
 
    def plot(self,axis='flux_pol_norm',colmap='blue',ax=False,fig=False):
        """
        plots rotation

        notes:
            axis - selects x axis of plot
            colmap - set the colour map (use get_cmap names)
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """
        
        import matplotlib
        from matplotlib import cm
        import matplotlib.pyplot as plt

        if ax is False:
            ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
        else:
            ax_flag=True

        if fig is False:
            fig_flag=False
        else:
            fig_flag=True

        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)
        ax.set_title(self.ID)
       
        ax.plot(self[axis],self['rotation'],color=colmap)
        ax.set_xlabel(axis)
        ax.set_ylabel('rotation')

        if ax_flag is False and fig_flag is False:
            plt.show()
 
#################################
 
##################################################################
 
###################################################################################################