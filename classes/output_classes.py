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


np.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays



















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
        
def dict_set(data,**kwargs):
    """
    generalised function (based upon safe_set) to set values in dictionary 'data' after checking source data exists 

    usage:
        dict_set(data,some_key=some_source,some_other_key=[1,2,3,4]) to set multiple values simultaneously
        dict_set(data,**{'some_key':100,'some_other_key':200}) equally
    """
    keys=kwargs.keys()
    values=kwargs.values()
    
    for key,value in zip(keys,values): #loop through kwargs
        if key is not None and value is not None: 
            data[key]=value

def safe_set(target,source):
    """
    generalised function to set a target value to source if it exists 
    """
    if source is not None:
        target=source













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

    def __init__(self,ID,data_format=None,filename=None,shot=None,run=None,properties=None): #this is common to all children (not overloaded), must have ID

        self.ID=ID #always set the ID, even if we don't invoke read_data i.e. a blank object is initialised
        if not none_check(self.ID,self.LOCUST_output_type,"read_data requires data_format, blank output initialised \n",data_format):
            self.read_data(data_format,filename,shot,run,properties)

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

    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None): #bad practice to change overridden method signatures, so retain all method arguments             
        """
        read data to be overloaded in all children classes
        """

        self.data={} #read_data definition is designed to be overloaded in children classes 

    def look(self):
        """
        print class information and data
        """

        print("\n-----------------------")
        print("ID - {ID}".format(ID=self.ID))  
        print("Output Type - {LOCUST_output_type}".format(LOCUST_output_type=self.LOCUST_output_type))

        if hasattr(self,'properties'): 
            print("Properties - {properties}".format(properties=self.properties))

        if hasattr(self,'data_format'):
            print("Data Format - {data_format}".format(data_format=self.data_format))
        
        if hasattr(self,'filename'):
            print("Input Filename - {filename}".format(filename=self.filename))
        
        if hasattr(self,'shot'):
            print("Shot - {shot}".format(shot=self.shot))
        
        if hasattr(self,'run'):
            print("Run - {run}".format(run=self.run))
        
        print("|")
        print("|")

        if hasattr(self,'data') and self.data:
            for key in self.data:
                if not self.data[key].size==0: #do not print if the data is empty
                    print("{key} - {value}".format(key=key,value=self.data[key]))
        
        print("-----------------------\n")

    def copy(self,target,*keys):
        """
        copy two output objects
 
        notes:
            if target.data is None or contains Nones then this function does nothing
            if no key supplied then copy all data over
            if key supplied then copy/append dictionary data accordingly
  
        usage:
            my_output.copy(some_other_output) to copy all data
            my_output.copy(some_other_output,'some_key','some_other_key') to copy specific fields
            my_output.copy(some_other_output, *some_list_of_args) equally
        """
        if none_check(self.ID,self.LOCUST_input_type,"cannot copy() - target.data is blank\n",target.data): #return warning if any target data contains empty variables
            pass

        elif not keys: #if empty, keys will be false i.e. no key supplied --> copy everything 
            self.data=copy.deepcopy(target.data) #using = with whole dictionary results in copying by reference, so need deepcopy() here
            if hasattr(target,'properties'): #copy properties field
                self.properties=target.properties

        elif not none_check(self.ID,self.LOCUST_input_type,"cannot copy() - found key containing None\n",*keys): 
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

        if none_check(self.ID,self.LOCUST_input_type,"cannot set() - empty key/value pair found\n",*allkeysvalues):
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






















################################################################## Orbits functions

def read_orbits_ASCII(filepath):
    """
    reads orbits stored in ASCII format - r phi z

    notes:
        reads in a headerline for number of particles
        reads in a footerline for number of time steps
    """

    print("reading orbits from ASCII")

    with open(filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_orbits_ASCII() cannot read from "+filepath)
    
    number_particles=int(lines[0]) #extract number of particles
    number_timesteps=int(lines[-1])-1 #extract number of time steps of each trajectory
    number_coords=int(3)

    del(lines[0])
    del(lines[-1])

    input_data = {} #initialise the dictionary to hold the dat
    input_data['orbits']=np.array([[float(value) for value in line.split()] for line in lines]).reshape((number_timesteps+1,number_coords,number_particles)) #read all data and reshape accordingly
    input_data['number_particles']=np.asarray(number_particles)
    input_data['number_timesteps']=np.asarray(number_timesteps)
   
    print("finished reading orbits from ASCII")

    return input_data

def dump_orbits_ASCII(output_data,filepath): 
    """
    writes orbits to ASCII format - r phi z
    
    notes:
        writes out a headerline for number of particles
        writes out a footerline for number of time steps
    """

    print("writing orbits to ASCII")

    with open(filepath,'w') as file: #open file

        file.write("{}\n".format(output_data['number_particles'].size)) #re-insert line containing number of particles

        for time_slice in range(output_data['number_timesteps']): #loop over everything again 
            for particle in range(output_data['number_particles']):

                    file.write("{r} {phi} {z}\n".format(r=output_data['orbits'][time_slice][particle][0],phi=output_data['orbits'][time_slice][particle][1],z=output_data['orbits'][time_slice][particle][2]))

        file.write("{}".format(output_data['number_timesteps'].size)) #re-insert line containing number of time steps

    print("finished writing orbits to ASCII")

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
        self.filename               name of file in output_files folder
        self.filepath               full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in output_files folder

    notes:
        data is stored such that coordinate i at time t for particle p is my_orbit['orbits'][t,i,p]
        in this way, a single particle's trajectory is my_orbit['orbits'][:,i,p] where i=0,1,2=r,phi,z
    """

    LOCUST_output_type='orbits'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None):
        """
        read orbits from file 

        notes:
        """

        if none_check(self.ID,self.LOCUST_output_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='ASCII': #here are the blocks for various file types, they all follow the same pattern
            if not none_check(self.ID,self.LOCUST_output_type,"cannot read_data from ASCII - filename required\n",filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties=properties
                self.data=read_orbits_ASCII(self.filepath) #read the file
        else:
            print("cannot read_data - please specify a compatible data_format (ASCII)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write orbits to file

        notes: 
        """
        if none_check(self.ID,self.LOCUST_output_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='ASCII':
            if not none_check(self.ID,self.LOCUST_output_type,"cannot dump_data to ASCII - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_orbits_ASCII(self.data,filepath)
        else:
            print("cannot dump_data - please specify a compatible data_format (ASCII)\n")






















################################################################## Final_Particle_List functions

def read_final_particle_list_ASCII(filepath='ptcl_cache.dat'):
    """
    reads final particle list stored in ASCII format

    notes:
        contains lots of references to process_ptcles.pro, written by Rob Akers
        status_flag describes each particle's final status (guide stored in status_flags, verbose guide in LOCUST/ctrk_mod.f90/ctrk_kernel) 
    """

    print("reading final particle list from ASCII")

    with open(filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_final_particle_list_ASCII() cannot read from "+filepath)

    #read in headerlines
    header=lines[0].split()
    n=int(header[0]) 
    ngpu=int(header[1]) #n*ngpu=number of particles
    niter=int(header[2]) #time iterations
    npt_=int(header[3]) #info slots
    nphc=int(header[4]) #levels in -DSPLIT split cache (always use the first)
    ntri=int(header[5]) #triangle grid "dimension" 

    #initialise particle_list and data dictionary
    input_data={}
    input_data['r']=np.array([])
    input_data['phi']=np.array([])
    input_data['z']=np.array([])
    input_data['v_r']=np.array([])
    input_data['v_phi']=np.array([])
    input_data['v_z']=np.array([])
    input_data['t']=np.array([])
    input_data['status_flag']=np.array([])
    input_data['additional_flag1']=np.array([])
    input_data['additional_flag2']=np.array([])
    input_data['additional_flag3']=np.array([])
    input_data['additional_flag4']=np.array([])
    input_data['additional_flag5']=np.array([])
    input_data['additional_flag6']=np.array([])
    input_data['additional_flag7']=np.array([])
    input_data['additional_flag8']=np.array([])
    input_data['additional_flag9']=np.array([])
    input_data['n']=n
    input_data['ngpu']=ngpu
    input_data['niter']=nniter
    input_data['npt_']=npt_
    input_data['nphc']=nphc
    input_data['ntri']=ntri
    input_data['number_particles']=n*ngpu

    #get rid of white space and completely flatten IDL/FORTRAN-style
    lines=[[float(number) for number in line.split()] for line in lines]
    lines=[number for line in lines for number  in line]
    del(lines[0:6])
    input_data['f']=lines[-1]

    for i in range(0,niter):
        
        #transfer chunk from lines to file_buffer and assimilate into dictionary
        file_buffer=np.array([lines[0:(n*ngpu)*npt_*nphc]]).reshape((n*ngpu),npt_,nphc,order='F')
        del(lines[0:(n*ngpu)*npt_*nphc])
        
        input_data['r']=np.append(input_data['r'],file_buffer[:,0,0])
        input_data['phi']=np.append(input_data['phi'],file_buffer[:,1,0])
        input_data['z']=np.append(input_data['z'],file_buffer[:,2,0])
        input_data['v_r']=np.append(input_data['v_r'],file_buffer[:,3,0])
        input_data['v_phi']=np.append(input_data['v_phi'],file_buffer[:,4,0])
        input_data['v_z']=np.append(input_data['v_z'],file_buffer[:,5,0])
        input_data['t']=np.append(input_data['t'],file_buffer[:,6,0])
        input_data['status_flag']=np.append(input_data['status_flag'],file_buffer[:,7,0])
        input_data['additional_flag1']=np.append(input_data['additional_flag1'],file_buffer[:,8,0])
        input_data['additional_flag2']=np.append(input_data['additional_flag2'],file_buffer[:,9,0])
        input_data['additional_flag3']=np.append(input_data['additional_flag3'],file_buffer[:,10,0])
        input_data['additional_flag4']=np.append(input_data['additional_flag4'],file_buffer[:,11,0])
        input_data['additional_flag5']=np.append(input_data['additional_flag5'],file_buffer[:,12,0])
        input_data['additional_flag6']=np.append(input_data['additional_flag6'],file_buffer[:,13,0])
        input_data['additional_flag7']=np.append(input_data['additional_flag7'],file_buffer[:,14,0])
        input_data['additional_flag8']=np.append(input_data['additional_flag8'],file_buffer[:,15,0])
        input_data['additional_flag9']=np.append(input_data['additional_flag9'],file_buffer[:,16,0])   

    input_data['status_flags']={} #nested dictionary to hold possible status flags for the particle list
    input_data['status_flags']['ok_if_greater']=0.0 
    input_data['status_flags']['undefined']=0.0
    input_data['status_flags']['left_space_grid']=-1.0
    input_data['status_flags']['not_poss_on_1st_call']=-1000.0 
    input_data['status_flags']['track_failure']=-2000.0
    input_data['status_flags']['unresolved_hit']=-3.0
    input_data['status_flags']['left_mesh']=-3000.0
    input_data['status_flags']['track_problem']=-4000.0
    input_data['status_flags']['ptcl_disconnect']=-5000.0
    input_data['status_flags']['PFC_intercept']=-5.0
    input_data['status_flags']['left_field_grid']=-6.0
    input_data['status_flags']['goose_fail']=-7000.0
    input_data['status_flags']['left_plasma']=-8.0
    input_data['status_flags']['thermalised']=-9.0
    input_data['status_flags']['coll_op_fail']=-10000.0
    input_data['status_flags']['GC_calc_fail']=-10.0 
    input_data['status_flags']['CX_loss']=-11.0
    input_data['status_flags']['gc_init_fail']=-11000.0
    input_data['status_flags']['bin_fail_soft']=-12.0
    input_data['status_flags']['bin_fail_hard_1']=-13000.0
    input_data['status_flags']['time_limit_reached']=-14.0
    input_data['status_flags']['cross_open_face']=-15.0
    input_data['status_flags']['bin_fail_hard_2']=-16000.0
    input_data['status_flags']['generic_fail_hard']=-99999.0
  
    print("finished reading final particle list from ASCII")

    return input_data


def dump_final_particle_list_ASCII(output_data,filepath): 
    """
    writes final particle list to ASCII format
    
    notes:

    """

    pass

################################################################## Final_Particle_List class

class Final_Particle_List(LOCUST_output):
    """
    class describing final particle list output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'final particle list'
    class data
        self.data_format            data format of original data e.g. ASCII
        self.filename               name of file in output_files folder
        self.filepath               full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in output_files folder

    notes:
        my_final_particle_list['status_flags'] contains a guide to the values a particle's status flag may contain 
    """

    LOCUST_output_type='final particle list'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None):
        """
        read final particle list from file 

        notes:
        """

        if none_check(self.ID,self.LOCUST_output_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='ASCII': #here are the blocks for various file types, they all follow the same pattern
            if not none_check(self.ID,self.LOCUST_output_type,"cannot read_data from ASCII - filename required\n",filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties=properties
                self.data=read_final_particle_list_ASCII(self.filepath) #read the file
        else:
            print("cannot read_data - please specify a compatible data_format (ASCII)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write final particle list to file

        notes: 
        """
        if none_check(self.ID,self.LOCUST_output_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='ASCII':
            if not none_check(self.ID,self.LOCUST_output_type,"cannot dump_data to ASCII - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_final_particle_list_ASCII(self.data,filepath)
        else:
            print("cannot dump_data - please specify a compatible data_format (ASCII)\n")















'''

################################################################## Distribution_Function functions

def read_distribution_function_bin(filepath,DTYPE=3,ITER=False,wtot=False,WIPE=False,test=False):
    """
    reads distribution function stored in unformatted fortran file

    notes:
        IDFTYP==1:F
                2:H
                3:X
                4:I
            else: #

        all n* variables are really n*-1 e.g. nVF is nVF-1 in LOCUST
    """

    file=FortranFile(filepath,'r')

    DFN_HEAD=filepath[0] #infer IDFTYP from first character of file name
    if DFN_HEAD=='F': 
        IDFTYP=1
    elif DFN_HEAD=='H':    
        IDFTYP=2
    elif DFN_HEAD=='X':
        IDFTYP=3
    elif DFN_HEAD=='I':
        IDFTYP=4
    else:
        pass
    
    input_data={} #initialise blank dictionary

    if DTYPE==3:

        #32-char checksum (32 bytes in Fortran, so essentially reading 8x32 bit floats)
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum

        #nVF-1
        input_data['nVF']=file.read_ints() #V

        #nPP-1
        input_data['nPP']=file.read_ints() #PPhi
        
        #nMU-1
        input_data['nMU']=file.read_ints() #MU
        
        #dVFh
        input_data['dVFh']=file.read_reals(dtype=np.float32) #XXX ARE THESE INTS?
        
        #dPPh
        input_data['dPPh']=file.read_reals(dtype=np.float32) #dPPh    = P_phi1h - P_phi0h
        
        #dMUh
        input_data['dMUh']=file.read_reals(dtype=np.float32)
        
        #V+dV/2 (nVF long)
        input_data['V+dV']=file.read_reals(dtype=np.float32)
        
        #PP+dPP/2 (nPP long)
        input_data['PP+dPP']=file.read_reals(dtype=np.float32) #PPhi
        
        #MU+MU/2 (nMU long)
        input_data['MU+dMU']=file.read_reals(dtype=np.float32)
        
        if wtot: #cumulative energy inventory (total energy injected so far)
            #Fh_norm (nVF by nPP by nMU)
            input_data['Fh_norm']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
        else:
            #Fh (NOTE check dimensions e.g. nVF by nPP by nMU)
            input_data['Fh']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
    
        #Fh_s (nVF by nPP by nMU) if Fh_s present
        input_data['Fh_s']=file.read_reals(dtype=np.float32) #Dfn. M.C. error
        
        #Jh
        input_data['Jh']=file.read_reals(dtype=np.float64) #Jacobian for IDFTYP=3 
        
        #Jh_s
        input_data['Jh_s']=file.read_ints(dtype=np.float64) #Jacobian error for IDFTYP=3
        
        #cpuTime
        input_data['cpu_time']=file.read_reals(dtype=np.float64) #this may not be present in the file
    

    #THIS IS THE TEST FILE - NOTE CAN INFER IF DTYPE=3 OR NOT DUE TO THIS VALUE (unless nvf is 0?)
    else:

        #32-char checksum
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum

        #just contains zero
        input_data['0_1']=file.read_ints()



        if IDFTYP==4:
            #nPSIF
            input_data['nPSIF']=file.read_ints() #number of surface contours
            #nPOLF
            input_data['nPOLF']=file.read_ints() #number of poloidal cells
        else:
            #nF
            input_data['nF']=file.read_ints() #R and Z dimension of the distribution function grid
            input_data['nF']=file.read_ints()

        if IDFTYP==1 or IDFTYP==4:
            #nLF
            input_data['nLF']=file.read_ints() #Vphi/V cell boundaries
        else:
            #nPP
            input_data['nPP']=file.read_ints() #PPhi          

        #nVF
        input_data['nVF']=file.read_ints() #V
        
        #nPF
        input_data['nPF']=file.read_ints() #poloidal gyro-phase cell boundaries

        if ITER:

            if WIPE:
                #Fh
                input_data['Fh']=[] ##Final combined DFn. grid
                for line in range(input_data['nPF']):
                    input_data['Fh'].append(file.read_reals(dtype=np.float32))
                
                input_data['nc']=len(input_data['Fh'])/input_data['nPF']

            else:
                input_data['Fh/iter']=[] ##Final combined DFn. grid  
                for line in range(input_data['nPF']):
                    input_data['Fh/iter'].append(file.read_reals(dtype=np.float32))
                
                input_data['nc']=len(input_data['Fh/iter'])/input_data['nPF']
        
        else:
            input_data['Fh']=[] ##Final combined DFn. grid   
            for line in range(input_data['nPF']):
                input_data['Fh'].append(file.read_reals(dtype=np.float32))
            
            input_data['nc']=len(input_data['Fh'])/input_data['nPF']
        
        #Fh_s
        input_data['Fh_s']=[] #Dfn. M.C. error
        for line in range(input_data['nPF']):
            input_data['Fh_s'].append(file.read_reals(dtype=np.float32))



        if IDFTYP==4: #XXX check the real() fortran function but I think real(some_number,single) can still be an int!

            #dVOL (nPSIF by nPOLF)
            input_data['dVOL']=file.read_reals(dtype=np.float32) 

            #npn
            input_data['npn']=file.read_ints() #cell granularity for IDFTYP=4          

            #csb (nPSIF-1 by nPOLF by 2*npn by 2)
            input_data['csb']=file.read_reals(dtype=np.float32)

            #npolh (nPSIF long)
            input_data['npolh']=file.read_ints()            

            #1.5 
            input_data['1.5 ']=file.read_reals(dtype=np.float32)

        else:

            #R+dR/2 (nF long)
            input_data['R+dR']=file.read_reals(dtype=np.float32)

            #Z+dZ/2 (nF long)
            input_data['Z+dZ']=file.read_reals(dtype=np.float32)

        if IDFTYP==1 or IDFTYP==4:

            #L+dL/2 (nLF long)
            input_data['L+dL']=file.read_reals(dtype=np.float32)

        else:

            #PP+dPP/2 (nPP long)
            input_data['PP+dPP']=file.read_reals(dtype=np.float32)

        #V+dV/2 (nVF long)
        input_data['V+dV']=file.read_reals(dtype=np.float32)

        #PG+dPG/2 (nPF long)
        input_data['PG+dPG']=file.read_reals(dtype=np.float32)

        input_data['Ab']=file.read_reals(dtype=np.float32) #fast ion masses
        input_data['Ai_1']=file.read_reals(dtype=np.float32) #first value of Ab
        input_data['Zb']=file.read_reals(dtype=np.float32) #trace particle Z
        input_data['Zi_1']=file.read_reals(dtype=np.float32) #first value of Zb
        input_data['Vsclh']=file.read_reals(dtype=np.float32) #vgrid upper bound
        input_data['Vnrm']=file.read_reals(dtype=np.float32)
        input_data['icoll']=file.read_ints() #collisions (1=on 0=off)
        input_data['iscat']=file.read_ints() #scattering (1=on 0=off)
        input_data['idiff']=file.read_ints() #energy diffusion (1=on 0=off)
        input_data['iloss']=file.read_ints() #charge exchange losses (1=on 0=off)
        input_data['iterm']=file.read_ints() #terminate if ptcl. leaves plasma
        input_data['niter']=file.read_ints() #number of iterations for isym=1 simulation
        input_data['LEIID']=file.read_ints() #integrator type(?)
        input_data['npnt']=file.read_ints() #points per gyration
        input_data['one_1']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['one_2']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['one_3']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['999']=file.read_reals(dtype=np.float32) #the number 999.0
        input_data['0_2']=file.read_reals(dtype=np.float32) #the maximum integrator step size
        input_data['dt0']=file.read_reals(dtype=np.float32) #the number 0.0
        input_data['tstp']=file.read_reals(dtype=np.float32) #timestep requested
        input_data['threadsPerBlock']=file.read_ints() #gpu threads per block 
        input_data['blocksPerGrid']=file.read_ints() #gpu blockers per grid

        if test:

            #Pdep/E0
            input_data['Pdep/E0']=file.read_reals(dtype=np.float32) #pdep is injected power

            #tau_s
            input_data['tau_s']=file.read_reals(dtype=np.float32) #zeroth order slowing down time

            #E0
            input_data['E0']=file.read_reals(dtype=np.float32) #energy (plasma frame)

            #EC
            input_data['EC']=file.read_reals(dtype=np.float32)

            #rho=Ai_1/(2*Ab)
            input_data['rho']=file.read_reals(dtype=np.float32)

            #siglg
            input_data['siglg']=file.read_reals(dtype=np.float32) #r.m.s. width of test src

    file.close()

    return input_data



def dump_distribution_function_bin(output_data,output_filepath): 
    """
    writes distribution function to binary
    
    notes:
        
    """

    pass

################################################################## Distribution_Function class

class Distribution_Function(LOCUST_output):
    """
    class describing orbits output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'orbits'
    class data
        self.data_format            data format of original data e.g. ASCII
        self.filename               name of file in output_files folder
        self.filepath               full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in output_files folder

    notes:
        data is stored such that
    """

    LOCUST_output_type='distribution_function'

    def read_data(self,data_format=None,input_filename=None,shot=None,run=None,properties=None):
        """
        read orbits from file 

        notes:
        """

        if none_check(self.ID,self.LOCUST_output_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='bin': #here are the blocks for various file types, they all follow the same pattern
            if not none_check(self.ID,self.LOCUST_output_type,"cannot read_data from bin - input_filename required\n",input_filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.input_filename=input_filename
                self.input_filepath=support.dir_output_files+input_filename
                self.properties=properties
                self.data=read_distribution_function_bin(self.input_filepath) #read the file
        else:
            print("cannot read_data - please specify a compatible data_format (bin)\n")            

    def dump_data(self,data_format=None,output_filename=None,shot=None,run=None):
        """
        write orbits to file

        notes: 
        """
        if none_check(self.ID,self.LOCUST_output_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='bin':
            if not none_check(self.ID,self.LOCUST_output_type,"cannot dump_data to bin - output_filename required\n",output_filename):
                output_filepath=support.dir_output_files+output_filename
                dump_distribution_function_bin(self.data,output_filepath)
        else:
            print("cannot dump_data - please specify a compatible data_format (bin)\n")



'''


#################################

##################################################################

###################################################################################################

