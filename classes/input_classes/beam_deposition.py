#beam_deposition.py
 
"""
Samuel Ward
02/11/2017
----
class to handle LOCUST beam deposition input data
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
    raise ImportError("ERROR: IMAS module could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes import utils
except:
    raise ImportError("ERROR: utils.py could not be imported!\nreturning\n")
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


################################################################## Beam_Deposition functions
 
def read_beam_depo_LOCUST(filepath):
    """
    reads birth profile stored in LOCUST format - r phi z v_r v_tor v_z

    notes:
    """
 
 
    print("reading beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+filepath)
     
    del(lines[0]) #first two lines are junk
    del(lines[0])
 
    input_data = {} #initialise the dictionary to hold the data
    input_data['r']=[] #initialise the arrays 
    input_data['phi']=[]
    input_data['z']=[]
    input_data['v_r']=[]
    input_data['v_tor']=[]
    input_data['v_z']=[]
 
    for line in lines:
 
        split_line=line.split()
        input_data['r'].append(float(split_line[0]))
        input_data['phi'].append(float(split_line[1]))
        input_data['z'].append(float(split_line[2]))
        input_data['v_r'].append(float(split_line[3]))
        input_data['v_tor'].append(float(split_line[4]))
        input_data['v_z'].append(float(split_line[5]))
 
    input_data['r']=np.asarray(input_data['r']) #convert to arrays
    input_data['phi']=np.asarray(input_data['phi'])
    input_data['z']=np.asarray(input_data['z'])
    input_data['v_r']=np.asarray(input_data['v_r'])
    input_data['v_tor']=np.asarray(input_data['v_tor'])
    input_data['v_z']=np.asarray(input_data['v_z'])

    print("finished reading beam deposition from LOCUST")
 
    return input_data
 
def dump_beam_depo_LOCUST(output_data,filepath):
    """
    writes birth profile to LOCUST format - r z phi v_r v_z v_tor
     
    notes:
        writes out two headerlines
    """
 
    print("writing beam deposition to LOCUST")

    with open(filepath,'w') as file: #open file
 
        file.write("{}\n".format(utils.fortran_string(1.0,13))) #re-insert junk lines
        file.write("{}\n".format(utils.fortran_string(1.0,13)))
 
        for this_particle in range(output_data['r'].size): #iterate through all particles i.e. length of our dictionary's arrays
 
            r_out=output_data['r'][this_particle] #briefly set to a temporary variable to improve readability
            phi_out=output_data['phi'][this_particle]
            z_out=output_data['z'][this_particle]
            v_r_out=output_data['v_r'][this_particle]
            v_tor_out=output_data['v_tor'][this_particle]
            v_z_out=output_data['v_z'][this_particle]
 
            file.write("{r}{phi}{z}{v_r}{v_tor}{v_z}\n".format(r=utils.fortran_string(r_out,14,6),phi=utils.fortran_string(phi_out,14,6),z=utils.fortran_string(z_out,14,6),v_r=utils.fortran_string(v_r_out,14,6),v_tor=utils.fortran_string(v_tor_out,14,6),v_z=utils.fortran_string(v_z_out,14,6)))
    
    print("finished writing beam deposition to LOCUST") 
 
def read_beam_depo_IDS(shot,run):
    """
    reads birth profile from a distribution_sources IDS and returns as a dictionary
 
    notes:
        reads in an arbitrary number of coordinates and injectors for each source        
        assumes that all sources have the same coordinate structure
        assumes markers hold only one time slice, at ...markers[0]
    """
    print("reading beam deposition from IDS")

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.open()
    input_IDS.distribution_sources.get() #open the file and get all the data from it

    input_data = {} #initialise blank dictionary to hold the data
    for identifier in input_IDS.distribution_sources.source[0].markers[0].coordinate_identifier: #generate keys for input_data by looking at the coordinates of the particle markers
        input_data[identifier.name.replace('\x00','').strip()]=[] #need to remove the unicode bits

    for source in input_IDS.distribution_sources.source: #cycle through all possible sources
        if len(source.markers[0].positions)>0:

            for coordinate_index in range(len(source.markers[0].positions[0,:])): #loop over the possible coordinate types e.g. r, phi, z
                coordinate_name=source.markers[0].coordinate_identifier[coordinate_index].name.replace('\x00','').strip()

                for marker in source.markers[0].positions[:,coordinate_index]: #this range should/must be the same for all values of coordinate_index

                    input_data[coordinate_name].extend([marker])    


    for key in input_data: #convert to numpy arrays
        input_data[key]=np.asarray(input_data[key])
 
    input_IDS.close()

    print("finished reading beam deposition from IDS")
 
    return input_data
 
def dump_beam_depo_IDS(ID,output_data,shot,run):
    """
    writes birth profile to a distribution_sources IDS
 
    notes:
    """
    
    print("writing beam deposition to IDS")

    output_IDS=imas.ids(shot,run) 
    output_IDS.create() #this will overwrite any existing IDS for this shot/run
 
    #write out code properties
    output_IDS.distribution_sources.ids_properties.comment=ID #write out identification
    output_IDS.distribution_sources.code.name="LOCUST_IO"
    output_IDS.distribution_sources.code.version=support.LOCUST_IO_version
    output_IDS.distribution_sources.ids_properties.homoegeneous_time=0   #must set homogeneous_time variable
     
    #add a type of source and add a time_slice for this source
    output_IDS.distribution_sources.source.resize(1) #adds a type of source here
    output_IDS.distribution_sources.source[0].markers.resize(1) #adds a time_slice here    
    output_IDS.distribution_sources.source[0].markers[0].time=0.0 #set the time of this time_slice
 
    #add definition of our coordinate basis - r,z,phi,v_r,v_z,v_tor in this case
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier.resize(1)
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].name="r" #name of coordinate
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].index=0 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].description="major radius [m]]" #description of coordinate

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].name="phi" 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].index=1 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].description="toroidal angle [rad]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=2 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="vertical position [m]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="v_r"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=3
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="radial velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="v_tor"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=4
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="toroidal velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="v_z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=5
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="vertical velocity [m/s]"

    #start storing particle data
    output_IDS.distribution_sources.source[0].markers[0].weights=np.ones(output_data['r'].size) #define the weights, i.e. number of particles per marker 
    positions=np.array([output_data['r'],output_data['phi'],output_data['z'],output_data['v_r'],output_data['v_tor'],output_data['v_z']]) #create 2D array of positions
    output_IDS.distribution_sources.source[0].markers[0].positions=np.transpose(positions) #swap the indices due to data dictionary convention
 
    #'put' all the output_data into the file and close
    output_IDS.distribution_sources.put()
    output_IDS.close()

    print("finished writing beam deposition to IDS")

def read_beam_depo_TRANSP(filepath):
    """
    reads birth profile from TRANSP ASCII file
 
    notes:
    """

    print("reading beam deposition from TRANSP format")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_TRANSP() cannot read from "+filepath)

    for counter,line in enumerate(lines): #look for the start of the data, marked by a certain string
        if str(line.split()[0])=='<start-of-data>':
            del(lines[0:counter+1])
            break

    input_data = {} #initialise the dictionary to hold the data
    input_data['x']=[] #initialise the arrays
    input_data['y']=[]
    input_data['z']=[]
    input_data['r']=[]  
    input_data['phi']=[]
    input_data['v_x']=[]
    input_data['v_y']=[]
    input_data['v_z']=[]
    input_data['v_r']=[]
    input_data['v_tor']=[]

    for line in lines:
        split_line=line.split()
        split_line[0]=float(split_line[0])
        split_line[1]=float(split_line[1])
        split_line[2]=float(split_line[2])
        split_line[3]=float(split_line[3])
        split_line[4]=float(split_line[4])
        split_line[5]=float(split_line[5])

        input_data['x'].append(split_line[0]) #only read in x,y,z with append for speed
        input_data['y'].append(split_line[1])
        input_data['z'].append(split_line[2])
        input_data['v_x'].append(split_line[3])
        input_data['v_y'].append(split_line[4])
        input_data['v_z'].append(split_line[5])

    input_data['x']=0.01*np.asarray(input_data['x']) #convert to arrays and from cm to m
    input_data['y']=0.01*np.asarray(input_data['y'])
    input_data['z']=0.01*np.asarray(input_data['z'])
    input_data['v_x']=0.01*np.asarray(input_data['v_x'])
    input_data['v_y']=0.01*np.asarray(input_data['v_y'])
    input_data['v_z']=0.01*np.asarray(input_data['v_z'])
    
    input_data['r']=np.asarray(np.sqrt(input_data['x']**2+input_data['y']**2)) #need to convert from x,y,z 
    input_data['v_r']=np.asarray(np.sqrt(input_data['v_x']**2+input_data['v_y']**2))
    input_data['phi']=np.asarray(np.arctan2(input_data['y'],input_data['x']))
    input_data['v_tor']=np.asarray(input_data['v_x']*np.sin(input_data['phi']))

    print("finished reading beam deposition from TRANSP format")

    return input_data


################################################################## Beam_Deposition class
 
class Beam_Deposition(base_input.LOCUST_input):
    """
    class describing neutral beam deposition profile input for LOCUST
 
    inZ_1Derited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'beam_deposition'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in input_files folder
        self.filepath               full path of file in input_files folder  
        self.shot                   shot number
        self.run                    run number
        self.properties             data to hold additional class-specific information e.g. ion species in Temperature
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in input_files folder
 
    notes:
        data is stored such that the coordinate 'r' for all particles is stored in my_beam_deposition['r']
        therefore the phase space position of particle p is:
            (my_beam_deposition['r'][p], my_beam_deposition['phi'][p], my_beam_deposition['z'][p], my_beam_deposition['v_r'][p], my_beam_deposition['v_tor'][p], my_beam_deposition['v_z'][p])
    """
 
    LOCUST_input_type='beam_deposition'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None): 
        """
        read beam_deposition from file 
 
        notes:
        """
 
        if utils.none_check(self.ID,self.LOCUST_input_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not utils.none_check(self.ID,self.LOCUST_input_type,"cannot read_data from LOCUST - filename required\n",filename): #must check we have all info required for reading GEQDSKs
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties=properties
                self.data=read_beam_depo_LOCUST(self.filepath) #read the file
         
        elif data_format=='IDS':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"cannot read_data from distribution_sources IDS - shot and run data required\n",shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties=properties
                self.data=read_beam_depo_IDS(self.shot,self.run)

        elif data_format=='TRANSP':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"cannot read_data from TRANSP - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties=properties
                self.data=read_beam_depo_TRANSP(self.filepath) #read the file
 
        else:
            print("cannot read_data - please specify a compatible data_format (LOCUST/IDS/TRANSP)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write beam_deposition to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run")
 
        if utils.none_check(self.ID,self.LOCUST_input_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"cannot dump_data to LOCUST - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_LOCUST(self.data,filepath)
         
        elif data_format=='IDS':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"cannot dump_data to distribution_sources IDS - shot and run required\n",shot,run):
                dump_beam_depo_IDS(self.ID,self.data,shot,run)
 
        else:
            print("cannot dump_data - please specify a compatible data_format (LOCUST/IDS)\n")

 
 
#################################
 
##################################################################
 
###################################################################################################