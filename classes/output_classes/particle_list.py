#particle_list.py

"""
Samuel Ward
15/01/2018
----
class to handle LOCUST particle list output data
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
    from scipy.io import  FortranFile
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


################################################################## Final_Particle_List functions

def read_final_particle_list_LOCUST(filepath='ptcl_cache.dat'):
    """
    reads final particle list stored in LOCUST format

    notes:
        contains lots of references to process_ptcles.pro, written by Rob Akers
        status_flag describes each particle's final status (guide stored in status_flags, verbose guide in LOCUST/ctrk_mod.f90/ctrk_kernel) 
    """

    print("reading final particle list from LOCUST")

    with open(filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_final_particle_list_LOCUST() cannot read from "+filepath)

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
    input_data['v_tor']=np.array([])
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
    input_data['n']=np.array(n)
    input_data['ngpu']=np.array(ngpu)
    input_data['niter']=np.array(niter)
    input_data['npt_']=np.array(npt_)
    input_data['nphc']=np.array(nphc)
    input_data['ntri']=np.array(ntri)
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
        input_data['v_tor']=np.append(input_data['v_tor'],file_buffer[:,4,0])
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
  
    print("finished reading final particle list from LOCUST")

    return input_data


def dump_final_particle_list_LOCUST(output_data,filepath): 
    """
    writes final particle list to LOCUST format
    
    notes:

    """

    pass

################################################################## Final_Particle_List class

class Final_Particle_List(base_output.LOCUST_output):
    """
    class describing final particle list output for LOCUST
    
    inZ_1Derited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'final particle list'
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
        my_final_particle_list['status_flags'] contains a guide to the values a particle's status flag may contain 
    """

    LOCUST_output_type='final particle list'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,properties=None):
        """
        read final particle list from file 

        notes:
        """

        if utils.none_check(self.ID,self.LOCUST_output_type,"cannot read_data - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not utils.none_check(self.ID,self.LOCUST_output_type,"cannot read_data from LOCUST - filename required\n",filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties=properties
                self.data=read_final_particle_list_LOCUST(self.filepath) #read the file
        else:
            print("cannot read_data - please specify a compatible data_format (LOCUST)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write final particle list to file

        notes: 
        """
        if utils.none_check(self.ID,self.LOCUST_output_type,"cannot dump_data - self.data and data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not utils.none_check(self.ID,self.LOCUST_output_type,"cannot dump_data to LOCUST - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_final_particle_list_LOCUST(self.data,filepath)
        else:
            print("cannot dump_data - please specify a compatible data_format (LOCUST)\n")


#################################

##################################################################

###################################################################################################
