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
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_input 
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/base_input.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


################################################################## rotation read functions
 
def read_rotation_LOCUST(filepath,**properties):
    """
    reads rotation profile stored in LOCUST format - normalised_poloidal_flux rotation_ang [rad/s]
    """
    
    print("reading rotation from LOCUST")

    with open(filepath,'r') as file:
         
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_rotation_LOCUST() cannot read from "+str(filepath))
     
        del(lines[0]) #first line contains the number of points
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['flux_pol_norm']=[] #initialise the arrays 
        input_data['rotation_ang']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['flux_pol_norm'].append(float(split_line[0]))
            input_data['rotation_ang'].append(float(split_line[1]))
     
        input_data['flux_pol_norm']=np.asarray(input_data['flux_pol_norm']) #convert to arrays
        input_data['rotation_ang']=np.asarray(input_data['rotation_ang'])
        
        print("finished reading rotation from LOCUST")

    return input_data

def read_rotation_UFILE(filepath,**properties):
    """
    reads rotation profiles from UFILE 

    notes:
        XXX untested
        currently no way of converting to rotation_ang
        if more than one time point then rotation field is multidimensional-[time_point,flux_pol_norm]        
    """

    print("reading rotation from UFILE")

    with open(filepath,'r') as file:

        input_data={}
        lines=file.readlines()

        for line_number,line in enumerate(lines): #extract dimensionality of the data
            split_line=line.split()

            if 'M/SEC' in split_line:
                in_ms=True        

            if 'X0' in split_line:
                length_rotation=int(split_line[0])
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

        input_data['flux_pol_norm']=np.array(data[:length_rotation])
        del(data[:length_rotation])
        input_data['time']=np.array(data[:length_time])
        del(data[:length_time])
        input_data['rotation_vel']=np.array(data).reshape(length_time,length_rotation)
        input_data['rotation_vel']=np.squeeze(input_data['rotation_vel']) #get rid of redundant axis if time is only 1 element long

    print("finished reading rotation from UFILE")

    return input_data

def read_rotation_ASCOT(filepath,**properties):
    """
    notes:
        XXX untested
    """

    with open(filepath,'r') as file:

        for line in file:
            if 'collision mode' in line:
                break
        line=file.readline()
        split_line=line.split()

        fields=[]
        desired_field='Vtor_I'

        for counter,field in enumerate(split_line[0::2]):
            if field==desired_field:
                desired_column_rotation=counter
            if field=='RHO':
                desired_column_flux_pol_norm_sqrt=counter

        input_data={}
        input_data['rotation_ang']=[]
        input_data['flux_pol_norm_sqrt']=[]

        for line in file:
            split_line=line.split()
            input_data['rotation_ang'].extend([float(line.split()[desired_column_rotation])])
            input_data['flux_pol_norm_sqrt'].extend([float(line.split()[desired_column_flux_pol_norm_sqrt])])

        input_data['rotation_ang']=np.asarray(input_data['rotation_ang'])
        input_data['flux_pol_norm_sqrt']=np.asarray(input_data['flux_pol_norm_sqrt'])
        input_data['flux_pol_norm']=input_data['flux_pol_norm_sqrt']**2

    return input_data

def read_rotation_excel_1(filepath,**properties):
    """
    reads rotation from excel spreadsheets supplied by Yueqiang Liu for ITER RMP study

    notes:  
        must include spreadsheet name holding minor radius in properties['sheet_name']
        must include spreadsheet name holding rotation in properties['sheet_name_rotation']
        must include name of rotation variable e.g. Vt(tF/tE=2) in properties['rotation_name']
        R_axis value is hardcoded here, please update accordingly
    """

    if 'sheet_name' not in properties: #must supply some sort of sheet_name
        print("ERROR: cannot read_rotation_excel_1 - properties['sheet_name'] must be set!\nreturning\n")
        return

    if 'rotation_name' not in properties: #must supply rotation_name
        print("ERROR: cannot read_rotation_excel_1 - properties['rotation_name'] must be set!\nreturning\n")
        return

    input_data={}
    input_data['flux_pol_norm'],radius_minor=run_scripts.utils.read_kinetic_profile_data_excel_1(filepath=filepath,x='Fp',y='a',sheet_name=properties['sheet_name'])
    input_data['rotation_vel']=run_scripts.utils.read_kinetic_profile_data_excel_1(filepath=filepath,y=properties['rotation_name'],sheet_name=properties['sheet_name_rotation'])
    input_data['flux_pol_norm_sqrt']=np.sqrt(input_data['flux_pol_norm'])
    input_data['rotation_vel']*=1000. #convert from km/s

    R_axis=6.2 #XXX warning this is hardcoded
    input_data['rmaj']=radius_minor+R_axis 
    input_data['rotation_ang']=input_data['rotation_vel']/input_data['rmaj']

    return input_data

################################################################## rotation write functions

def dump_rotation_LOCUST(output_data,filepath,**properties):
    """
    writes rotation profile to LOCUST format - normalised_poloidal_flux rotation_ang [rad/s]    
     
    notes:
        writes out a headerline for length of file
    """

    print("writing rotation to LOCUST")
 
    with open(filepath,'w') as file: #open file

        normalised_flux=np.abs(output_data['flux_pol_norm']) #take abs
        normalised_flux,output_rot=processing.utils.sort_arrays(normalised_flux,output_data['rotation_ang']) #check order
 
        file.write("{}\n".format(processing.utils.fortran_string(output_rot.size,12))) #re-insert line containing length
        
        for point in range(output_rot.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm}{rot}\n".format(flux_pol_norm=processing.utils.fortran_string(normalised_flux[point],16,8),rot=processing.utils.fortran_string(output_rot[point],16,8)))
 
    print("finished writing rotation to LOCUST")

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

def dump_rotation_IDS(ID,output_data,shot,run,**properties):
    """
    notes:
        performing operations as per NEMO
        assumes all species are comprised of a single element
    """
    
    print("writing rotation to IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_rotation_IDS could not import IMAS module!\nreturning\n")
        return

    output_IDS=imas.ids(int(shot),int(run)) 
    output_IDS.open_env(username,imasdb,'3') #open the IDS
    output_IDS.core_profiles.get()
 
    #write out code properties
    output_IDS.core_profiles.ids_properties.comment=ID #write out identification
    output_IDS.core_profiles.code.name="LOCUST_IO"
    if settings.commit_hash_default_LOCUST_IO: output_IDS.core_profiles.code.commit=str(settings.commit_hash_default_LOCUST_IO)
    output_IDS.core_profiles.code.version=support.LOCUST_IO_version
    output_IDS.core_profiles.ids_properties.homogeneous_time=1   #must set homogeneous_time variable
    output_IDS.core_profiles.time=np.array([0.0])
     
    #add a time_slice and set the time
    if len(output_IDS.core_profiles.profiles_1d)==0: output_IDS.core_profiles.profiles_1d.resize(1) #add a time_slice
    output_IDS.core_profiles.profiles_1d[0].time=0.0 #set the time of the time_slice
 
    #write out rotation
    if 'A' not in properties or 'Z' not in properties:
        print("ERROR: could not dump_rotation_IDS - properties['A'] and properties['Z'] required!\nreturning\n")
        return
    else:

        species_number=[counter for counter,ion in enumerate(
                        output_IDS.core_profiles.profiles_1d[0].ion)   
                        if len(ion.element)!=0
                        if ion.element[0].a==properties['A']
                        if ion.element[0].z_n==properties['Z']]
        if not species_number:
            species_number=[-1] #if matching ion found in IDS, its number now held in species number
            output_IDS.core_profiles.profiles_1d[0].ion.resize(len(output_IDS.core_profiles.profiles_1d[0].ion)+1) #add an ion species if desired species does not already exist in IDS
            output_IDS.core_profiles.profiles_1d[0].ion[species_number[0]].element.resize(1)
            output_IDS.core_profiles.profiles_1d[0].ion[species_number[0]].element[0].a=properties['A']
            output_IDS.core_profiles.profiles_1d[0].ion[species_number[0]].element[0].z_n=properties['Z']

        output_IDS.core_profiles.profiles_1d[0].ion[species_number[0]].rotation_frequency_tor=output_data['rotation_ang']
            
    #write out the axes
    output_IDS.core_profiles.profiles_1d[0].grid.psi=output_data['flux_pol']

    #'put' all the output_data into the file and close
    output_IDS.core_profiles.put()
    output_IDS.close()

    print("finished writing rotation to IDS")

################################################################## rotation class
 
class Rotation(classes.base_input.LOCUST_input):
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

        elif data_format=='EXCEL1': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from EXCEL1 - filename required\n".format(self.ID),filename): 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_rotation_excel_1(self.filepath,**properties)    

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/EXCEL1)\n".format(self.ID))            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write rotation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID = {})".format(self.ID))
            
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

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to core_profiles IDS - shot, run and ion species property required\n".format(self.ID),shot,run,self.properties):
                dump_rotation_IDS(self.ID,self.data,shot,run,**properties)
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/MARSF/IDS)\n".format(self.ID))
 
    def plot(self,axis='flux_pol_norm',colmap=cmap_default,colmap_val=np.random.uniform(),ax=False,fig=False):
        """
        plots rotation

        notes:
            axis - selects x axis of plot
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
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
       
        ax.plot(self[axis],self['rotation_ang'],color=colmap(colmap_val))
        ax.set_xlabel(axis)
        ax.set_ylabel('rotation_ang')

        if ax_flag is False and fig_flag is False:
            plt.show()
 
#################################
 
##################################################################
 
###################################################################################################