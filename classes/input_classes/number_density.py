#number_density.py
 
"""
Samuel Ward
02/11/2017
----
class to handle LOCUST number_density input data
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
    import scipy
    import numpy as np
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d #import 3D plotting axes
    from mpl_toolkits.mplot3d import Axes3D
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

################################################################## Number_Density read functions
 
def read_number_density_LOCUST(filepath):
    """
    reads number density profile stored in LOCUST format - normalised_poloidal_flux n [m^-3]
 
    notes:
        reads in a headerline for length of file
    """

    print("reading number density from LOCUST")

    with open(filepath,'r') as file:
         
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_number_density_LOCUST() cannot read from "+filepath)
     
        del(lines[0]) #first line contains the number of points
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['flux_pol_norm']=[] #initialise the arrays 
        input_data['n']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['flux_pol_norm'].append(float(split_line[0]))
            input_data['n'].append(float(split_line[1]))
     
        input_data['flux_pol_norm']=np.asarray(input_data['flux_pol_norm']) #convert to arrays
        input_data['n']=np.asarray(input_data['n'])

    print("finished reading number density from LOCUST")
    
    return input_data
 
def read_number_density_IDS(shot,run,**properties):
    """
    reads LOCUST number density data from a core_profiles IDS and returns as a dictionary
 
    notes:
        
    """
 
    print("reading number density from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_number_density_IDS could not import IMAS module!\nreturning\n")
        return

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.open()
    input_IDS.core_profiles.get() #open the file and get all the data from it
 
    input_data = {} #initialise blank dictionary to hold the data
    
    #read in number density depending on species
    if properties['species']=='electrons':
        input_data['n']=np.asarray(input_IDS.core_profiles.profiles_1d[0].electrons.density)
    elif properties['species']=='ions':
        input_data['n']=np.asarray(input_IDS.core_profiles.profiles_1d[0].ion[0].density)
    else:
        print("ERROR: cannot read_number_density_IDS - Number_Density.properties['species'] must be set to 'electrons' or 'ions'\n")
        return 

    #read in axes
    processing.utils.dict_set(input_data,flux_pol=np.asarray(input_IDS.core_profiles.profiles_1d[0].grid.psi)/(2.0*pi)) #convert to Wb/rad
    processing.utils.dict_set(input_data,flux_tor_coord=np.asarray(input_IDS.core_profiles.profiles_1d[0].grid.rho_tor))
    processing.utils.dict_set(input_data,q=np.asarray(input_IDS.core_profiles.profiles_1d[0].q))
    if input_IDS.core_profiles.vacuum_toroidal_field.b0: #if we are supplied a vacuum toroidal field to derive toroidal flux, then derive it
        processing.utils.dict_set(input_data,flux_tor=np.asarray(input_IDS.core_profiles.vacuum_toroidal_field.b0*(input_data['flux_tor_coord']**2)/2.)) #in Wb/rad

    input_IDS.close()
    print("finished reading number density from IDS")
 
    return input_data
 
################################################################## Number_Density write functions

def dump_number_density_LOCUST(output_data,filepath):
    """
    writes number density profile to LOCUST format - normalised_poloidal_flux n [m^-3]
     
    notes:
        writes out a headerline for length of file
    """
 
    print("writing number density to LOCUST")

    with open(filepath,'w') as file: #open file

        normalised_flux=np.abs(output_data['flux_pol_norm']) #take abs
        normalised_flux,output_n=processing.utils.sort_arrays(normalised_flux,output_data['n']) #check order
 
        file.write("{}\n".format(processing.utils.fortran_string(output_data['flux_pol_norm'].size,12))) #re-insert line containing length
        for point in range(normalised_flux.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm}{n}\n".format(flux_pol_norm=processing.utils.fortran_string(normalised_flux[point],16,8),n=processing.utils.fortran_string(output_n[point],16,8)))
 
    print("finished writing number density to LOCUST")

def dump_number_density_IDS(ID,output_data,shot,run,**properties):
    """
    writes LOCUST number density data to a core_profiles IDS
    """
    
    print("writing number density to IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_number_density_IDS could not import IMAS module!\nreturning\n")
        return

    output_IDS=imas.ids(shot,run) 
    output_IDS.create() #this will overwrite any existing IDS for this shot/run
 
    #write out code properties
    output_IDS.core_profiles.ids_properties.comment=ID #write out identification
    output_IDS.core_profiles.code.name="LOCUST_IO"
    output_IDS.core_profiles.code.version=support.LOCUST_IO_version
    output_IDS.core_profiles.ids_properties.homoegeneous_time=0   #must set homogeneous_time variable
     
    #add a time_slice and set the time
    output_IDS.core_profiles.profiles_1d.resize(1) #add a time_slice
    output_IDS.core_profiles.profiles_1d[0].time=0.0 #set the time of the time_slice
 
    #write out number density depending on species
    if properties['species']=='electrons':
        output_IDS.core_profiles.profiles_1d[0].electrons.density=output_data['n']
    elif properties['species']=='ions':
        output_IDS.core_profiles.profiles_1d[0].ion.resize(1) #add an ion species 
        #TODO need to add additional species data here e.g. mass, charge
        output_IDS.core_profiles.profiles_1d[0].ion[0].density=output_data['n']
    else:
        print("cannot dump_number_density_IDS - Number_Density.properties['species'] must be set to 'electrons' or 'ions'\n")

    #write out the axes
    processing.utils.safe_set(output_IDS.core_profiles.profiles_1d[0].grid.psi,output_data['flux_pol'])
    processing.utils.safe_set(output_IDS.core_profiles.profiles_1d[0].grid.rho_tor,output_data['flux_tor_coord'])
    processing.utils.safe_set(output_IDS.core_profiles.profiles_1d[0].q,output_data['q'])

    #'put' all the output_data into the file and close
    output_IDS.core_profiles.put()
    output_IDS.close()

    print("finished writing number density to IDS")

def dump_number_density_MARSF(output_data,filepath):
    """
    writes number density data to MARSF Mogui ASCII format 

    note:
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    """

    print("writing number density to MARSF mogui")

    with open(filepath,'w') as file: #open file

        flux_pol_norm_sqrt=np.sqrt(np.abs(output_data['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,output_n=processing.utils.sort_arrays(flux_pol_norm_sqrt,output_data['n']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(flux_pol_norm_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(flux_pol_norm_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm_sqrt} {n}\n".format(flux_pol_norm_sqrt=processing.utils.fortran_string(flux_pol_norm_sqrt[point],24,18),n=processing.utils.fortran_string(output_n[point],24,18)))

    print("finished writing number density to MARSF mogui")

################################################################## Number_Density class
 
class Number_Density(classes.base_input.LOCUST_input):
    """
    class describing number density profile input for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'number density'
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
        data is stored such that a reading of number density at coordinate 'coord' is:
            my_number_density['flux_tor'][coord], my_number_density['n'][coord]
    """
 
    LOCUST_input_type='number_density'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read number density from file 
 
        notes:
        """
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"Number_Density.properties['species'] not specified - set to 'electrons' or 'ions' for IDS functionality\n",properties):
            pass
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_number_density_LOCUST(self.filepath) #read the file
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from core_profiles IDS - shot, run and ion species property required\n",shot,run,properties):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_number_density_IDS(self.shot,self.run,**self.properties)
 
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST/IDS)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write number density to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)

        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_number_density_LOCUST(self.data,filepath)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to core_profiles IDS - shot, run and ion species property required\n",shot,run,self.properties):
                dump_number_density_IDS(self.ID,self.data,shot,run,**properties)

        elif data_format=='MARSF':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to MARSF - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_number_density_MARSF(self.data,filepath)
 
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST/IDS/MARSF)\n")

    def plot(self,axis='flux_pol_norm',colmap='blue',ax=False,fig=False):
        """
        plots number density
         
        notes:
            axis - selects x axis of plot
            colmap - set the colour map (use get_cmap names)
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

        ax.plot(self[axis],self['n'],color=colmap)
        ax.set_xlabel(axis)
        ax.set_ylabel('number density [m^-3]')
        ax.set_title(self.ID)
        
        if ax_flag is False and fig_flag is False:
            plt.show()
 
#################################
 
##################################################################
 
###################################################################################################