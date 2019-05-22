#process_output.py

'''
Samuel Ward
25/1/2018
----
processing routines for LOCUST output data
---
notes:
    https://stackoverflow.com/questions/9455111/python-define-method-outside-of-class-definition
---
'''

##################################################################
#Preamble

import sys

try:
    import scipy.integrate
    import numpy as np
    import pathlib
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
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


##################################################################
#Main Code

def dfn_transform(some_dfn,axes=['R','Z']):
    """
    transforms and integrates the distribution function according to pre-defined configurations 
    
    args:
        axes - the dimensions over which to transform the DFN to
    notes:
        remember dimensions of unedited dfn are my_dfn['dfn'][P,V/E,V_pitch,R,Z]
        assumes unedited dfn
        assumes the bin widths for a given dimension are constant
        assumes toroidal symmetry (no toroidal dimension in dfn)
        if an array of indices is given, then slice the dfn accordingly and return without any integration
            note for an infinite slice, axes will need to contain slice() objects e.g. axes=[0,0,0,slice(None),slice(None)] for all R,Z values
        dimension P is meaningless in EBASE mode 

    axes options:
        R,Z - integrate over pitch, gyrophase and velocity [m]^-3
        E,V_pitch - integrate over space and transform to [eV]^-1[dpitch]^-1 
        E - [eV]^-1 
        R - [m]^-3
        N - total # 
        list of indices and slices
    """

    dfn=copy.deepcopy(some_dfn) #make deep copy here since functions designed to repeatedly take fresh DFNs would otherwise permanently change it

    #begin list of specific options

    if dfn.properties['EBASE'] is True: #if LOCUST dfn is against energy

        if axes==['R','Z']:
            dfn['dfn']*=dfn['dE']*dfn['dV_pitch'] #integrate
            for counter in range(3): #sum
                dfn['dfn']=np.sum(dfn['dfn'],axis=0)

        elif axes==['E','V_pitch']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.*constants.pi*dfn['dR']*dfn['dZ']
            #then need to integrate over the unwanted coordinates
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over Z
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over R
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over P

        elif axes==['E']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.*constants.pi*dfn['dR']*dfn['dZ']
            dfn['dfn']*=dfn['dV_pitch'] #integrate over pitch
            
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #sum over Z
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #sum over R
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #sum over V_pitch
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over P

        elif axes==['R']:
            dfn['dfn']*=dfn['dE']*dfn['dV_pitch'] #integrate
            for counter in range(3): #sum over gyrophase, pitch and energy
                dfn['dfn']=np.sum(dfn['dfn'],axis=0)
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #sum over Z

        elif axes==['V_pitch']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.*constants.pi*dfn['dR']*dfn['dZ']*dfn['dE']
            #then need to integrate over the unwanted coordinates
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over Z
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over R
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over P
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over E

        elif axes==['N']:
            #applying full Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.*constants.pi*dfn['dR']*dfn['dZ']*dfn['dV_pitch']*dfn['dE']
            for all_axes in range(dfn['dfn'].ndim): #sum over all dimensions
                dfn['dfn']=np.sum(dfn['dfn'],axis=0) 

        #general option
        elif len(axes)==dfn['dfn'].ndim: #if user supplies all axes then slice WITHOUT integrating
            dfn['dfn']=dfn['dfn'][tuple(axes)]
            #XXX need to then reset dfn['nV'],dfn['R'] etc data here?
        else:
            print("ERROR: dfn_transform given invalid axes arguement: "+str(axes))

    else: #if LOCUST dfn is against velocity

        if axes==['R','Z']:
            #apply velocity space Jacobian
            for v in range(len(dfn['V'])):
                dfn['dfn'][:,v,:,:,:]*=dfn['V'][v]**2
            dfn['dfn']*=dfn['dV']*dfn['dV_pitch']*dfn['dP']

            #then need to integrate over the first 3 dimensions which we do not need
            for counter in range(3):
                dfn['dfn']=np.sum(dfn['dfn'],axis=0) #sum over gyrophase then V then V_pitch

        elif axes==['E','V_pitch']:
            #applying velocity space and gyrophase Jacobian
            for v in range(len(dfn['V'])):
                dfn['dfn'][:,v,:,:,:]*=dfn['V'][v]
            dfn['dfn']*=dfn['dP']*constants.species_charge/(constants.species_mass)

            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.0*constants.pi*dfn['dR']*dfn['dZ']

            #then need to integrate over the unwanted coordinates
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over gyrophase
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over Z
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over R

        elif axes==['E']:
            #applying velocity space and gyrophase Jacobian
            for v in range(len(dfn['V'])):
                dfn['dfn'][:,v,:,:,:]*=dfn['V'][v]
            dfn['dfn']*=dfn['dP']*dfn['dV_pitch']*constants.species_charge/(constants.species_mass)

            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.0*constants.pi*dfn['dR']*dfn['dZ']

            #then need to integrate over the unwanted coordinates
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over gyrophase
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over Z
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over R
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over V_pitch

        elif axes==['R']:
            #apply velocity space Jacobian
            for v in range(len(dfn['V'])):
                dfn['dfn'][:,v,:,:,:]*=dfn['V'][v]**2
            dfn['dfn']*=dfn['dV']*dfn['dV_pitch']*dfn['dP']

            #then need to integrate over the first 3 dimensions which we do not need
            for counter in range(3):
                dfn['dfn']=np.sum(dfn['dfn'],axis=0) #sum over gyrophase then V then V_pitch
            dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #sum over Z

        elif axes==['N']:
            #apply velocity space Jacobian
            for v in range(len(dfn['V'])):
                dfn['dfn'][:,v,:,:,:]*=dfn['V'][v]**2
            dfn['dfn']*=dfn['dV']*dfn['dV_pitch']*dfn['dP']

            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn['R'])):
                dfn['dfn'][:,:,:,r,:]*=dfn['R'][r]*2.0*constants.pi*dfn['dR']*dfn['dZ']

            #sum over all axes
            for all_axes in range(dfn['dfn'].ndim):
                dfn['dfn']=np.sum(dfn['dfn'],axis=0)

        #general option
        elif len(axes)==dfn['dfn'].ndim: #if user supplies all axes then slice WITHOUT integrating
            dfn['dfn']=dfn['dfn'][tuple(axes)]
            #XXX need to then reset dfn['nV'],dfn['R'] etc data here?
        else:
            print("ERROR: dfn_transform given invalid axes arguement: "+str(axes))


    return dfn

def dfn_crop(some_dfn,**kwargs):
    """
    notes:
        warning! to work, dfn_index and 1D dfn axes must accurately reflect dfn which is still stored e.g. dfn['dfn'][r,z] must contain dfn['R'],dfn['Z'] and dfn['dfn_index']=['R','Z']
        always maintains shape of dfn
    args:
        kwargs - axes and their limits 
    usage:
        new_dfn=dfn_crop(R=[1]) generates dfn at point closest to R=1
        new_dfn=dfn_crop(R=[0,1]) crops dfn between 0<R<1
    """

    dfn=copy.deepcopy(some_dfn)

    keys=list(kwargs.keys())
    values=list(kwargs.values())

    for key,value in zip(keys,values):
        if key not in dfn['dfn_index']:
            print("ERROR: dfn_crop supplied invalid axis name ({}) - see ['dfn_index'] for possible axes".format(key))    
        else:

            dimension_to_edit=dfn['dfn_index'].tolist().index(key) #figure out which dimension we are cropping over
            
            if len(value)==2: #user has supplied range
                i=np.where((value[0]<dfn[key])&(dfn[key]<value[1])) #get new indices which satisfy range
                i=i[0] #get first element of returned tuple
            elif len(value)==1: #user has supplied single value for nearest neighbour
                i=[np.abs(dfn[key]-value[0]).argmin()]

            dfn[key]=dfn[key][i] #crop 1D arrays accordingly
            nkey='n{}'.format(key) #reset associated nkey values too e.g. reset nR if cropping R
            dfn[nkey]=np.array(len(dfn[key]))

            dfn['dfn']=np.moveaxis(dfn['dfn'],dimension_to_edit,0) #move desired axis of dfn array to front to crop
            dfn['dfn']=dfn['dfn'][i] #crop dfn
            dfn['dfn']=np.moveaxis(dfn['dfn'],0,dimension_to_edit) #move desired axis of dfn array back to original position             

    return dfn

def particle_list_compression(filepath,coordinates=['R','phi','Z','time','status_flag'],dump=False):
    """
    opens huge particle lists in memory-efficient way

    args:
        coordinates - the particle coordinates to read in
        dump - toggle to re-dump to ASCII afterwards (NOTE: NOT YET IMPLEMENTED)
    notes:
        this code will break if the file line length > number of entries for each coordinate
        currently only reads the first phc values i.e. phc index = 0  
    """

    print("compressing final particle list file: "+str(filepath))

    indices_coordinate={} #and need a dictionary to refer to locations of quantities in file
    indices_coordinate['0']='R'
    indices_coordinate['1']='phi'
    indices_coordinate['2']='Z'
    indices_coordinate['3']='V_R'
    indices_coordinate['4']='V_tor'
    indices_coordinate['5']='V_Z'
    indices_coordinate['6']='time'
    indices_coordinate['7']='status_flag'
    indices_coordinate['8']='additional_flag1'
    indices_coordinate['9']='additional_flag2'
    indices_coordinate['10']='additional_flag3'
    indices_coordinate['11']='additional_flag4'
    indices_coordinate['12']='additional_flag5'
    indices_coordinate['13']='additional_flag6'
    indices_coordinate['14']='additional_flag7'
    indices_coordinate['15']='additional_flag8'
    indices_coordinate['16']='additional_flag9'

    input_data={} #need to initialise dictionary to hold data we read from file

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

    for coordinate in coordinates: #set up arrays to hold the data
        input_data[coordinate]=np.array([])
        
    with open(filepath) as file:
        file_buffer=[] #to hold 1D chunks of the file 
        dimension_counter=0 #since the file is laid out as R_1 R_2 ... R_n phi_1 phi_2 ... phi_n where n-> particle number, can keep track of which dimension we're currently reading 
                
        line=file.readline().split() #read header line first
        n=int(line[0]) 
        ngpu=int(line[1]) #n*ngpu=number of particles
        niter=int(line[2]) #time iterations
        npt_=int(line[3]) #info slots
        nphc=int(line[4]) #levels in -DSPLIT split cache (always use the first)
        ntri=int(line[5]) #triangle grid "dimension" 

        for line_number,line in enumerate(file): #read line-by-line for memory efficiency
            line=line.split()
            
            file_buffer.extend([float(number) for number in line]) #read a line and add to the file buffer 
            if len(file_buffer) > n*ngpu: #if we have read in all values for this coordinate
                
                coordinate_completed=indices_coordinate[str(dimension_counter)] #check which coordinate we just finished reading in
                if coordinate_completed in coordinates: #if we want this coordinate then keep
                    input_data[coordinate_completed]=np.append(input_data[coordinate_completed],file_buffer[0:n*ngpu])

                del(file_buffer[0:n*ngpu]) #clear the bit of the buffer which we have read in (might still contain some of next coordinate)
                dimension_counter+=1
                if dimension_counter==17: #stop after we read the first 16*number_particles values i.e. only phc=0

                    input_data['n']=np.array(n)
                    input_data['ngpu']=np.array(ngpu)
                    input_data['niter']=np.array(niter)
                    input_data['npt_']=np.array(npt_)
                    input_data['nphc']=np.array(nphc)
                    input_data['ntri']=np.array(ntri)
                    input_data['number_particles']=np.array(n*ngpu)

                    print("finished compressing final particle list file: "+str(filepath))

                    return input_data

'''
def extract_DFN_particle_list(some_particle_list,some_equilibrium,some_bins=None):
    """
    generates a new distribution function object from data stored in a final particle list

    notes:
        produces unweighted binned distribution function
    args:
        some_particle_list - 
        some_equilibrium - 
        some_bins - 
    """

    V_pitch=processing.utils.pitch_calc_2D(some_particle_list,some_equilibrium)
    dummy_gyrophase=np.zeros(len(some_particle_list['R'])) #if uniform in gyrophase
    if 'V' not in some_particle_list.data:
        some_particle_list['V']=np.array(np.sqrt(some_particle_list['V_R']**2+some_particle_list['V_tor']**2+some_particle_list['V_Z']**2))

    if some_bins:
        dfn=np.histogramdd((dummy_gyrophase,some_particle_list['V'],some_particle_list['V_pitch'],some_particle_list['R'],some_particle_list['Z']))        
        edges=
        centres=
        d=
    else:
        dfn=np.histogramdd((dummy_gyrophase,some_particle_list['V'],some_particle_list['V_pitch'],some_particle_list['R'],some_particle_list['Z']))
        edges=
        centres=
        d=

    dfn/=d  ##apply real space jacobian
    dfn/=d  ##apply velocity space jacobian

    P,V,V_pitch,R,Z
'''
#################################

##################################################################

###################################################################################################