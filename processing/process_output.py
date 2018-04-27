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

##################################################################
#Main Code

def dfn_integrate(some_dfn):
    """
    integrate a distribution function to get particles/cell and real velocity

    notes:
        returns the Dfn array
    """

    dfn=copy.deepcopy(some_dfn['dfn'])

    for p in range(int(some_dfn['nP'])): #denormalise the velocity dimension
        for v in range(int(some_dfn['nV'])):
            dfn[p,v,:,:,:]*=some_dfn['V'][v]**2

    for r in range(int(some_dfn['nR'])): #convert from markers per s^3/m^6 to markers per cell
        if some_dfn['nP']>1:
            dfn[:,:,:,:,r]*=2.0*pi*some_dfn['R'][r]*(some_dfn['R'][1]-some_dfn['R'][0])*(some_dfn['Z'][1]-some_dfn['Z'][0])*(some_dfn['V'][1]-some_dfn['V'][0])*(some_dfn['L'][1]-some_dfn['L'][0])*(some_dfn['P'][1]-some_dfn['P'][0])
        else:
            dfn[:,:,:,:,r]*=2.0*pi*some_dfn['R'][r]*(some_dfn['R'][1]-some_dfn['R'][0])*(some_dfn['Z'][1]-some_dfn['Z'][0])*(some_dfn['V'][1]-some_dfn['V'][0])*(some_dfn['L'][1]-some_dfn['L'][0])*2.*pi

    return dfn

def dfn_collapse(some_dfn,coordinates=['R','Z']):
    """
    integrate and collapse the dfn to the supplied to coordinates

    notes:
        returns the Dfn array
    """

    dfn=copy.deepcopy(some_dfn['dfn'])

    coordinate_indices=[] #collapse dfn along specified axes
    if 'Z' not in coordinates: #if Z is not in the coordinates we want to keep
        coordinate_indices.extend([4]) #then mark it as a dimension to integrate over
    if 'R' not in coordinates:
        coordinate_indices.extend([3])
    if 'L' not in coordinates: #pitch
        coordinate_indices.extend([2]) 
    if 'V' not in coordinates: #velocity
        coordinate_indices.extend([1])
    if 'P' not in coordinates: #special s
        coordinate_indices.extend([0])

    for coordinate in coordinate_indices: #axis denotes which coordinate will be collapsed, so go in descending to get array shape correct
        dfn=np.sum(dfn,axis=coordinate)
    
    return dfn

def particle_list_compression(filepath,coordinates=['R','phi','Z','V_R','V_tor','V_Z','status_flag'],dump=False):
    """
    opens huge particle lists in memory-efficient way

    notes:
        coordinates - the particle coordinates to read in
        dump - toggle to re-dump to ASCII afterwards (NOTE: NOT YET IMPLEMENTED)

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