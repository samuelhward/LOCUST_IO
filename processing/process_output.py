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

try:
    import scipy.integrate
    import numpy as np
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
    
pi=np.pi
e_charge=1.602e-19 #define electron charge
mass_neutron=1.661e-27 #define mass of neutron

##################################################################
#Main Code

def dfn_transform(some_dfn,axes=['R','Z']):
    """
    transforms and integrates the distribution function according to pre-defined configurations 
    
    args:
        axes - the dimensions over which to transform the DFN to
    notes:
        remember dimensions of dfn are my_dfn['dfn'][P,V,V_pitch,R,Z]
        always assumes the bin widths for a given dimension are constant
        assumes deuterium in energy conversion to energy space
        always assumes toroidal symmetry (no toroidal dimension in dfn)
        if an array of indices is given, then slice the dfn accordingly and return without any manipulation
            note for an infinite slice, axes will need to contain slice() objects e.g. axes=[0,0,0,slice(None),slice(None)] for all R,Z values

    axes options:
        R,Z - integrate over pitch, gyrophase, velocity and toroidal angle [m]^-2
        E,V_pitch - [eV]^-1[pitch bin]^-1  
    """

    dfn=copy.deepcopy(some_dfn) #make deep copy here since functions designed to repeatedly take fresh DFNs would otherwise permanently change it

    #begin list of specific options

    if axes==['R','Z']:
        #integrating over the velocity space
        for v in range(int(some_dfn['nV'])):
            dfn['dfn'][:,v,:,:,:]*=some_dfn['V'][v]**2
        dfn['dfn']*=some_dfn['dV']*some_dfn['dV_pitch']*some_dfn['dP']

        #then need to collapse over the first 3 dimensions which we do not need
        for counter in range(3):
            dfn['dfn']=np.sum(dfn['dfn'],axis=0) #sum over gyrophase then V then V_pitch


    elif axes==['E','V_pitch']:
        #integrating over gyrophase and applying velocity space Jacobian
        for v in range(int(some_dfn['nV'])):
            dfn['dfn'][:,v,:,:,:]*=some_dfn['V'][v]
        dfn['dfn']*=some_dfn['dP']*e_charge/(2.*mass_neutron)

        #integrating over all real space
        for r in range(int(some_dfn['nR'])):
            dfn['dfn'][:,:,:,r,:]*=some_dfn['R'][r]*2.0*pi*some_dfn['dR']*some_dfn['dZ']

        #then need to collapse over the unwanted coordinates
        dfn['dfn']=np.sum(dfn['dfn'],axis=0) #over gyrophase
        dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over Z
        dfn['dfn']=np.sum(dfn['dfn'],axis=-1) #over R

    #general option
    
    elif len(axes)==dfn['dfn'].ndim: #if user supplies all axes then slice
        dfn['dfn']=dfn['dfn'][tuple(axes)]
        #XXX need to then reset dfn['nV'],dfn['R'] etc data here?
    else:
        print("ERROR: dfn_transform given invalid axes arguement: "+str(axes))

    return dfn

def dfn_crop(some_dfn,**kwargs):
    """
    notes:
        assumes full 3D dfn
    args:
        kwargs - axes and their limits e.g R=[0,1] crops dfn between 0<R<1
    """

    dfn=copy.deepcopy(some_dfn)

    keys=list(kwargs.keys())
    values=list(kwargs.values())

    for key,value in zip(keys,values):
        if key not in dfn['dfn_index']:
            print("ERROR: dfn_crop supplied invalid axis name - see ['dfn_index'] for possible axes")    
        else:

            dimension_to_edit=dfn['dfn_index'].tolist().index(key) #figure out which dimension we are cropping over

            i=(np.where((value[0]<dfn[key])&(dfn[key]<value[1]))) #get new indices which satisfy range

            dfn[key]=dfn[key][i] #crop 1D arrays accordingly

            dfn['dfn']=np.moveaxis(dfn['dfn'],dimension_to_edit,0) #move desired axis of dfn array to front to crop
            dfn['dfn']=dfn['dfn'][i,:,:,:,:] #crop dfn
            dfn['dfn']=np.moveaxis(dfn['dfn'],0,dimension_to_edit) #move desired axis of dfn array back to original position             

    return dfn

def particle_list_compression(filepath,coordinates=['R','phi','Z','V_R','V_tor','V_Z','status_flag'],dump=False):
    """
    opens huge particle lists in memory-efficient way

    args:
        coordinates - the particle coordinates to read in
        dump - toggle to re-dump to ASCII afterwards (NOTE: NOT YET IMPLEMENTED)
    notes:
        this code will break if the file line length > number of entries for each coordinate
        currently only reads the first phc values i.e. phc index = 0  
    """

    print("compressing final particle list file: "+filepath)

    indices_coordinate={} #and need a dictionary to refer to locations of quantities in file
    indices_coordinate['0']='R'
    indices_coordinate['1']='phi'
    indices_coordinate['2']='Z'
    indices_coordinate['3']='V_R'
    indices_coordinate['4']='V_tor'
    indices_coordinate['5']='V_Z'
    indices_coordinate['6']='t'
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

    for coordinate in coordinates: #set up arrays to hold the data
        input_data[coordinate]=np.array([])
        
    with open(filepath) as file:
        file_buffer=[] #to hold 1D chunks of the file 
        dimension_counter=0 #since the file is laid out as R_1 R_2 ... R_n phi_1 phi_2 ... phi_n where n-> particle number, can keep track of which dimension we're currently reading 
        
        for line_number,line in enumerate(file): #read line-by-line for memory efficiency
            line=line.split()
            
            if line_number==0: #we are still on the header line
                n=int(line[0]) 
                ngpu=int(line[1]) #n*ngpu=number of particles
                niter=int(line[2]) #time iterations
                npt_=int(line[3]) #info slots
                nphc=int(line[4]) #levels in -DSPLIT split cache (always use the first)
                ntri=int(line[5]) #triangle grid "dimension" 

            else:
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
                        input_data['number_particles']=n*ngpu

                        print("finished compressing final particle list file: "+filepath)

                        return input_data

#################################

##################################################################

###################################################################################################