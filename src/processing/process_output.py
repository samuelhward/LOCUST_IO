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
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
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
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
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
            print("ERROR: dfn_transform given invalid axes argument: "+str(axes))

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
            print("ERROR: dfn_transform given invalid axes argument: "+str(axes))


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

'''

XXX needs to be adapted to read in different nesting order I think - I think file is all particles' R value, THEN all particles' phi value and so on...
XXX 99% certain this is true based on the reshaping in particle_list and the observed quantities for the compression toggle on and off (first value of each matches)
XXX also nphc>1 usually...

def particle_list_compression(filepath,coordinates=['R','phi','Z','time','dt','V_R','V_phi','V_Z'],dump=False):
    """
    reads selected data from particle lists in memory-efficient way

    args:
        filepath - full path of file
        coordinates - the particle coordinates to read in
        dump - toggle to re-dump to ASCII afterwards (NOTE: NOT YET IMPLEMENTED)
    notes:
        currently only reads the first phc values i.e. phc index = 0  
    """

    print("compressing final particle list file: "+str(filepath))

    coordinates,dump=run_scripts.utils.literal_eval(coordinates,dump) #in case taking command line input and coordinates are string

    indices_coordinate={} #and need a dictionary to refer to locations of quantities in file
    indices_coordinate['R']=0
    indices_coordinate['phi']=1
    indices_coordinate['Z']=2
    indices_coordinate['V_R']=3
    indices_coordinate['V_phi']=4
    indices_coordinate['V_Z']=5
    indices_coordinate['time']=6
    indices_coordinate['dt']=7
    indices_coordinate['FG']=8
    indices_coordinate['tet']=9
    indices_coordinate['psi']=10
    indices_coordinate['V_R_next']=11
    indices_coordinate['V_phi_next']=12
    indices_coordinate['V_Z_next']=13
    indices_coordinate['R_next']=14
    indices_coordinate['phi_next']=15
    indices_coordinate['Z_next']=16

    input_data={} #need to initialise dictionary to hold data we read from file

    input_data['status_flags']={} #nested dictionary to hold possible status flags for the particle list
    input_data['status_flags'][0.0]='undefined'
    input_data['status_flags'][-1.0]='left_space_grid'
    input_data['status_flags'][-1000.0]='not_poss_on_1st_call '
    input_data['status_flags'][-2000.0]='track_failure'
    input_data['status_flags'][-3.0]='unresolved_hit'
    input_data['status_flags'][-3000.0]='left_mesh'
    input_data['status_flags'][-4000.0]='track_problem'
    input_data['status_flags'][-4.0]='PFC_intercept_2D'
    input_data['status_flags'][-5000.0]='ptcl_disconnect'
    input_data['status_flags'][-5.0]='PFC_intercept_3D'
    input_data['status_flags'][-6.0]='left_field_grid'
    input_data['status_flags'][-7000.0]='goose_fail'
    input_data['status_flags'][-8.0]='left_plasma'
    input_data['status_flags'][-9.0]='thermalised'
    input_data['status_flags'][-10000.0]='coll_op_fail'
    input_data['status_flags'][-10.0]='GC_calc_fail '
    input_data['status_flags'][-11.0]='CX_loss'
    input_data['status_flags'][-11000.0]='gc_init_fail'
    input_data['status_flags'][-12.0]='bin_fail_soft'
    input_data['status_flags'][-13000.0]='bin_fail_hard_1'
    input_data['status_flags'][-14.0]='time_limit_reached'
    input_data['status_flags'][-15.0]='cross_open_face'
    input_data['status_flags'][-16000.0]='bin_fail_hard_2'
    input_data['status_flags'][-99999.]='generic_fail_hard'

    def status_flags_dispatch(value):
        try:
            return input_data['status_flags'][value]
        except:
            return 'ok'
    
    with open(filepath) as file:
        
        file_buffer=[] #to hold 1D chunks of the file 
        particle_number=0 #since the file is laid out as R_1 Phi_1 Z_1...R_n Phi_n Z_n... where n-> particle number, can keep track of which particke we are currently reading 
                
        line=file.readline().split() #read header line first
        n=int(line[0]) 
        ngpu=int(line[1]) #n*ngpu=number of particles
        niter=int(line[2]) #time iterations
        npt_=int(line[3]) #info slots
        nphc=int(line[4]) #levels in -DSPLIT split cache (always use the first)
        ntri=int(line[5]) #triangle grid "dimension" 

        for coordinate in coordinates: #set up arrays to hold the data
            input_data[coordinate]=np.zeros(n*ngpu*niter)
        input_data['status_flag']=[]

        line_number=0
        line=True
        while line:
        
            while len(file_buffer)<npt_: #always make sure buffer has at least one particle's worth of data in
                line=file.readline() #read line-by-line for memory efficiency
                line_number+=1        
                file_buffer.extend([float(number) for number in line.split()]) #read a line and add to the file buffer 
            
            #we have read in all values for this particle and know the first N entries in file_buffer correspond to the particle coordinates
            for coordinate in coordinates: #read the requested data from the buffer
                input_data[coordinate][particle_number]=file_buffer[indices_coordinate[coordinate]]
                if coordinate is 'dt': input_data['status_flag'].append(status_flags_dispatch(input_data['dt'][particle_number])) #determine status_flag of marker

            del(file_buffer[0:npt_]) #delete this particle's data from the buffer
            
            particle_number+=1

            if particle_number==n*ngpu*niter: #stop after we read the all the particles

                input_data['n']=np.array(n)
                input_data['ngpu']=np.array(ngpu)
                input_data['niter']=np.array(niter)
                input_data['npt_']=np.array(npt_)
                input_data['nphc']=np.array(nphc)
                input_data['ntri']=np.array(ntri)
                input_data['number_particles']=np.array(npt_)
                input_data['status_flag']=np.array(input_data['status_flag'])

                try:
                    input_data['weight']=np.full(len(input_data['R']),1.)
                except:
                    pass

                if all([quant in input_data for quant in ['V_R','V_phi','V_Z']]):
                    input_data['E']=.5*constants.species_mass*(input_data['V_R']**2+input_data['V_phi']**2+input_data['V_Z']**2)/constants.species_charge

                print("finished compressing final particle list file: "+str(filepath))

                return input_data
'''

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
        some_particle_list['V']=np.array(np.sqrt(some_particle_list['V_R']**2+some_particle_list['V_phi']**2+some_particle_list['V_Z']**2))

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