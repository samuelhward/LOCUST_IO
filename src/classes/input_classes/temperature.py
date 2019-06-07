#temperature.py
 
"""
Samuel Ward
02/11/2017
----
class to handle LOCUST temperature input data
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
    import classes.input_classes.equilibrium 
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/equilibrium.py could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


################################################################## Temperature read functions
 
def read_temperature_LOCUST(filepath,**properties):
    """
    reads temperature profile stored in LOCUST format - normalised_poloidal_flux T(ev)
    """
    
    print("reading temperature from LOCUST")

    with open(filepath,'r') as file:
         
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_temperature_LOCUST() cannot read from "+str(filepath))
     
        del(lines[0]) #first line contains the number of points
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['flux_pol_norm']=[] #initialise the arrays 
        input_data['T']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['flux_pol_norm'].append(float(split_line[0]))
            input_data['T'].append(float(split_line[1]))
     
        input_data['flux_pol_norm']=np.asarray(input_data['flux_pol_norm']) #convert to arrays
        input_data['T']=np.asarray(input_data['T'])
        
        print("finished reading temperature from LOCUST")

    return input_data
 
def read_temperature_LOCUST_h5(filepath,**properties):
    """

    notes:
    """

    print("reading temperature from LOCUST_h5")
    
    try:
        import h5py
    except:
        raise ImportError("ERROR: read_temperature_LOCUST_h5 could not import h5py module!\n") 
        return

    with h5py.File(filepath,'r') as file:
        
        input_data={}

        input_data['flux_pol_norm']=file['Input Data/Kinetic Data/Profiles (1D)/PSIn'].value

        if properties['species']=='electrons':
            input_data['T']=file['Input Data/Kinetic Data/Profiles (1D)/Te/data'].value
        elif properties['species']=='ions':
            input_data['T']=file['Input Data/Kinetic Data/Profiles (1D)/Ti/data'].value
        else:
            print("ERROR: cannot read_temperature_LOCUST_h5 - properties['species'] must be set to 'electrons' or 'ions'\n") 
            return

    print("finished reading temperature from LOCUST_h5")

    return input_data

def read_temperature_IDS(shot,run,**properties):
    """
    reads relevant LOCUST temperature data from a core_profiles IDS and returns as a dictionary
 
    notes:
        
    """
 
    print("reading temperature from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_temperature_IDS could not import IMAS module!\nreturning\n")
        return

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.open_env(username,imasdb,'3')
    input_IDS.core_profiles.get() #open the file and get all the data from it
 
    input_data = {} #initialise blank dictionary to hold the data
    
    #read in temperature depending on species
    if properties['species']=='electrons':
        input_data['T']=np.asarray(input_IDS.core_profiles.profiles_1d[0].electrons.temperature)
    elif properties['species']=='ions':
        input_data['T']=np.asarray(input_IDS.core_profiles.profiles_1d[0].ion[0].temperature)
    else:
        print("ERROR: cannot read_temperature_IDS - properties['species'] must be set to 'electrons' or 'ions'\n")
        return
     
    #read in axes
    processing.utils.dict_set(input_data,flux_pol=np.asarray(input_IDS.core_profiles.profiles_1d[0].grid.psi)/(2.0*constants.pi)) #convert to Wb/rad
    processing.utils.dict_set(input_data,flux_tor_coord=np.asarray(input_IDS.core_profiles.profiles_1d[0].grid.rho_tor))
    processing.utils.dict_set(input_data,q=np.asarray(input_IDS.core_profiles.profiles_1d[0].q))
    if input_IDS.core_profiles.vacuum_toroidal_field.b0: #if we are supplied a vacuum toroidal field to derive toroidal flux, then derive it
        processing.utils.dict_set(input_data,flux_tor=np.asarray(input_IDS.core_profiles.vacuum_toroidal_field.b0*(input_data['flux_tor_coord']**2)/2.)) #in Wb/rad

    input_IDS.close()

    print("finished reading temperature from IDS")
 
    return input_data
 
def read_temperature_UDA(shot,time,**properties):
    """
    reads ion or electron temperature profiles from CCFE UDA database

    notes:
        this script heavily relies on the structure of UDA - please contant Stuart Henderson for updates
    """

    print("reading temperature from UDA")

    try:
        import pyuda
        udaClient=pyuda.Client()
        getdata=udaClient.get
    except:
        raise ImportError("ERROR: read_temperature_UDA could not import pyuda!\nreturning\n")

    def replace_nan(quantity,time_index):
        """
        quick little function to replace NaNs in data with 1D interpolated values
        """

        def nan_helper(y): 
            return np.isfinite(y), lambda z:z.nonzero()[0]

        nansid, xdim = nan_helper(quantity.data[time_index,:])
        try:
            quantity.data[time_index,~nansid] = np.interp(xdim(~nansid),xdim(nansid),quantity.data[time_index,nansid])
        except:
            quantity.data[time_index,~nansid] = np.full(len(~nansid),-1.)

        return

    equilibrium=classes.input_classes.equilibrium.Equilibrium(ID='',data_format='UDA',shot=shot,time=time) #needs to read corresponding equilibrium
    input_data={}

    if properties['species']=='electrons':

        if shot<23000:
            signal_T='atm_te'
            signal_R='atm_R'
        else:
            signal_T='ayc_te'
            signal_R='ayc_r'

        T=getdata(signal_T,shot)
        R_data=getdata(signal_R,shot) #electron and ion temperature signals are stored against different radial bases - must access another dataset for electrons
        time_grid=T.dims[0].data
        time_index=np.abs(time_grid-time).argmin() #figure out what time we are wanting to output (pick closest)
        R_grid=R_data.data[time_index]
    
    elif properties['species']=='ions':

        T=getdata('act_ss_temperature',shot)
        time_grid=T.dims[0].data
        time_index=np.abs(time_grid-time).argmin() #figure out what time we are wanting to output (pick closest)
        R_grid=T.dims[1].data

    replace_nan(T,time_index)
    T=T.data[time_index,:]

    flux_pol_grid=processing.utils.RZ_to_Psi(R_grid,np.full(len(R_grid),0.0),equilibrium) #assume measurements are along Z=0
    flux_pol_grid_norm=(flux_pol_grid-equilibrium['simag'])/(equilibrium['sibry']-equilibrium['simag'])
    flux_pol_grid_norm,T=processing.utils.sort_arrays(flux_pol_grid_norm,T) #take all data points along the line of sight

    input_data['time']=time_grid[time_index]
    input_data['flux_pol_norm']=flux_pol_grid_norm #or could interpolate here onto new grid
    input_data['T']=T

    print("finished reading temperature from UDA")

    return input_data

def read_temperature_UFILE(filepath,**properties):
    """
    reads number density profiles from UFILE 

    notes:
        if more than one time point then temperature field is multidimensional-[time_point,flux_pol_norm]        
    """

    print("reading temperature from UFILE")

    with open(filepath,'r') as file:

        input_data={}
        lines=file.readlines()

        for line_number,line in enumerate(lines): #extract dimensionality of the data
            split_line=line.split()

            if 'TemperatureKEV' in split_line or 'KEV' in split_line:
                in_kev=True        

            if 'X0' in split_line:
                length_temperature=int(split_line[0])
            elif 'X1' in split_line:
                length_time=int(split_line[0])
                del(lines[:line_number+1]) #delete all lines up to this point
                del(lines[-1]) #remove the last two lines with Ufile authorship
                del(lines[-1])
                break    

        data=[]
        for line_number,line in enumerate(lines): #extract data now
            split_line=line.split()   
            for number in split_line:
                data.append(float(number))

        input_data['flux_pol_norm']=np.array(data[:length_temperature])
        del(data[:length_temperature])
        input_data['time']=np.array(data[:length_time])
        del(data[:length_time])
        input_data['T']=np.array(data).reshape(length_time,length_temperature)
        input_data['T']=np.squeeze(input_data['T']) #get rid of redundant axis if time is only 1 element long

        if in_kev:
            input_data['T']*=1e3

    print("finished reading temperature from UFILE")

    return input_data

def read_temperature_ASCOT(filepath,**properties):
    """
    notes:
        must include 'species' in properties - either 'electrons' or 'ions'
        include integer species_number in properties to select species (1 by default)
    """

    with open(filepath,'r') as file:

        for line in file:
            if 'collision mode' in line:
                break
        line=file.readline()
        split_line=line.split()

        fields=[]

        if properties['species']=='electrons':
            desired_field='Te'
        elif properties['species']=='ions':
            if 'species_number' in properties.keys():
                desired_field='Ti'+str(properties['species_number'])    
            else:
                desired_field='Ti1'
        else:
            print("ERROR: cannot read_number_density_ASCOT - properties['species'] must be set to 'electrons' or 'ions' ")
            return 

        for counter,field in enumerate(split_line[0::2]):
            if field==desired_field:
                desired_column_T=counter
            if field=='RHO':
                desired_column_flux_pol_norm_sqrt=counter

        input_data={}
        input_data['T']=[]
        input_data['flux_pol_norm_sqrt']=[]
        
        for line in file:
            split_line=line.split()
            input_data['T'].extend([float(line.split()[desired_column_T])])
            input_data['flux_pol_norm_sqrt'].extend([float(line.split()[desired_column_flux_pol_norm_sqrt])])

        input_data['T']=np.asarray(input_data['T'])
        input_data['flux_pol_norm_sqrt']=np.asarray(input_data['flux_pol_norm_sqrt'])
        input_data['flux_pol_norm']=input_data['flux_pol_norm_sqrt']**2

    return input_data

################################################################## Temperature write functions

def dump_temperature_LOCUST(output_data,filepath,**properties):
    """
    writes temperature profile to LOCUST format - normalised_poloidal_flux T(ev)    
     
    notes:
        writes out a headerline for length of file
    """

    print("writing temperature to LOCUST")
 
    with open(filepath,'w') as file: #open file

        normalised_flux=np.abs(output_data['flux_pol_norm']) #take abs
        normalised_flux,output_T=processing.utils.sort_arrays(normalised_flux,output_data['T']) #check order
 
        file.write("{}\n".format(processing.utils.fortran_string(output_T.size,8))) #re-insert line containing length
        
        for point in range(output_T.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm}{T}\n".format(flux_pol_norm=processing.utils.fortran_string(normalised_flux[point],16,8),T=processing.utils.fortran_string(output_T[point],16,8)))
 
    print("finished writing temperature to LOCUST")

def dump_temperature_IDS(ID,output_data,shot,run,**properties):   
    """
    writes relevant LOCUST temperature data to a core_profiles IDS
    """

    print("writing temperature to IDS")
 
    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_temperature_IDS could not import IMAS module!\nreturning\n")
        return

    output_IDS=imas.ids(shot,run) 
    output_IDS.create_env(username,imasdb,'3') #this will overwrite any existing IDS for this shot/run
 
    #write out code properties
    output_IDS.core_profiles.ids_properties.comment=ID #write out identification
    output_IDS.core_profiles.code.name="LOCUST_IO"
    output_IDS.core_profiles.code.version=support.LOCUST_IO_version
    output_IDS.core_profiles.ids_properties.homoegeneous_time=0   #must set homogeneous_time variable
     
    #add a time_slice and set the time
    output_IDS.core_profiles.profiles_1d.resize(1) #add a time_slice
    output_IDS.core_profiles.profiles_1d[0].time=0.0 #set the time of the time_slice

    #write out temperature depending on species
    if properties['species']=='electrons':
        output_IDS.core_profiles.profiles_1d[0].electrons.temperature=output_data['T']
    elif properties['species']=='ions':
        output_IDS.core_profiles.profiles_1d[0].ion.resize(1) #add an ion species 
        #TODO need to add additional species data here e.g. mass, charge
        output_IDS.core_profiles.profiles_1d[0].ion[0].temperature=output_data['T']
    else:
        print("ERROR: cannot dump_temperature_IDS - properties['species'] must be set to 'electrons' or 'ions'\n")
        return
        
    #write out the axes
    processing.utils.safe_set(output_IDS.core_profiles.profiles_1d[0].grid.psi,output_data['flux_pol'])


    #'put' all the output_data into the file and close
    output_IDS.core_profiles.put()
    output_IDS.close()

    print("finished writing temperature to IDS")

def dump_temperature_MARSF(output_data,filepath,**properties):
    """
    writes temperature profile to MARSF Mogui ASCII format
     
    notes:
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    """

    print("writing temperature to MARSF mogui")
 
    with open(filepath,'w') as file: #open file

        flux_pol_norm_sqrt=np.sqrt(np.abs(output_data['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,output_T=processing.utils.sort_arrays(flux_pol_norm_sqrt,output_data['T']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(flux_pol_norm_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(flux_pol_norm_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays 
            file.write("{flux_pol_norm}{T}\n".format(flux_pol_norm=processing.utils.fortran_string(flux_pol_norm_sqrt[point],24,18),T=processing.utils.fortran_string(output_T[point],24,18)))
 
    print("finished writing temperature to MARSF mogui")

################################################################## Temperature class
 
class Temperature(classes.base_input.LOCUST_input):
    """
    class describing temperature profile input for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'temperature'
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
 
    LOCUST_input_type='temperature'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,time=None,**properties):
        """
        read temperature from file 
 
        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"Temperature.properties['species'] not specified - set to 'electrons' or 'ions' for IDS/LOCUST_h5 functionality\n",properties['species']):
            pass
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_temperature_LOCUST(self.filepath,**properties)
    
        elif data_format=='LOCUST_h5': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_h5 - filename and ion species property required\n".format(self.ID),filename,properties['species']): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename #remember this is in output files director!
                self.properties={**properties}
                self.data=read_temperature_LOCUST_h5(self.filepath,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from core_profiles IDS - shot, run and ion species property required\n".format(self.ID),shot,run,properties['species']):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_temperature_IDS(self.shot,self.run,**properties)
 
        elif data_format=='UDA':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from UDA - shot and time required\n".format(self.ID),shot,time):
                self.data_format=data_format
                self.shot=shot
                self.time=time
                self.properties={**properties}
                self.data=read_temperature_UDA(self.shot,self.time,**properties)

        elif data_format=='UFILE': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from UFILE - filename required\n".format(self.ID),filename): 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_temperature_UFILE(self.filepath,**properties)

        elif data_format=='ASCOT': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT - filename required\n".format(self.ID),filename): 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_temperature_ASCOT(self.filepath,**properties)              

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/LOCUST_h5/IDS/UDA/UFILE/ASCOT)\n".format(self.ID))            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write temperature to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID={})".format(self.ID))
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_temperature_LOCUST(self.data,filepath)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to core_profiles IDS - shot, run and ion species property required\n".format(self.ID),shot,run,properties['species']):
                dump_temperature_IDS(self.ID,self.data,shot,run,**properties)

        elif data_format=='MARSF':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to MARSF - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_temperature_MARSF(self.data,filepath,**properties)
 
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/IDS/MARSF)\n".format(self.ID))
 
    def plot(self,axis='flux_pol_norm',colmap=cmap_default,ax=False,fig=False):
        """
        plots temperature

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
       
        ax.plot(self[axis],self['T'],color=colmap(np.random.uniform()))
        ax.set_xlabel(axis)
        ax.set_ylabel('temperature [eV]')

        if ax_flag is False and fig_flag is False:
            plt.show()
 
#################################
 
##################################################################
 
###################################################################################################