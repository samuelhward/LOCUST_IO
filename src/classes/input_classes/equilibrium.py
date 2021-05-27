#equilibrium.py
 
'''
Samuel Ward
02/11/2017
----
class to handle LOCUST equilibrium input data
---
usage:
    see README.md for usage
 
notes:         
---
'''


###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import re
    import ast
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)   
try:
    import processing.process_input 
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/process_input.py could not be imported!\nreturning\n")
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
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## Equilibrium read functions

def read_equilibrium_GEQDSK(filepath,**properties): 
    """ 
    generic function for reading a G-EQDSK-formatted equilibrium file
 
    notes:
        originally written by Ben Dudson and edited by Nick Walkden
    """

    input_data = {}

    def file_numbers(ingf):#instance of generators are objects hence can use .next() method on them
        """
        generator to read numbers in a file
     
        notes:
            originally written by Ben Dudson
     
            calling get_next() on a generator produced with this function will permanently advance the iterator in this generator
            when generators reach the end, that's permanently it!
        """
        toklist = []
        while True: #runs forever until break statement (until what is read is not a line) below will cause a return - which permanently stops a generator
            line = ingf.readline() #read in line from INGF
            if not line: break #break if readline() doesnt return a line i.e. end of file
            line = line.replace("NaN","-0.00000e0") #replaces NaNs with 0s
            pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?' #regular expression to find numbers
            toklist = re.findall(pattern,line) #toklist now holds all the numbers in a file line (is essentially a single file line)
            for tok in toklist: #yield every number in that one line individually
                yield tok #so in terms of iteration, using get_next() will cause another line to be read

    print("reading equilibrium from GEQDSK")
    
    with open(filepath,'r') as file: #open file
 
        line = file.readline() #first line should be case, id number and dimensions
        if not line:
            raise IOError("ERROR: read_equilibrium_GEQDSK() cannot read from "+str(filepath))
         
        #extract case, id number and dimensions  
        conts = line.split()    #split by white space (no argument in .split())
        input_data['nZ_1D'] = np.asarray(int(conts[-1])) #same as nyefit or height dimension
        input_data['nR_1D'] = np.asarray(int(conts[-2])) #same as nxefit or width dimension
        input_data['idum'] = np.asarray(int(conts[-3]))
        #input_data['time'] = np.asarray(int(consts[-4])) #time slice time (in ms)
        #input_data['EFITshot'] = np.asarray(int(consts[-5])) #shot number

        flags = {} #read the case in (basically just ID data)
        flags['case'] = conts[0:-4] 
     
        #now use generator to read numbers from remaining lines in file
        token = file_numbers(file) #token now holds all lines in the file, containing all values. each processing.utils.get_next() call will grab the next number in a line (since processing.utils.get_next() returns the next value in the yield loop? check this plz)
         
        float_keys = [
        'rdim','zdim','rcentr','rleft','zmid',
        'rmaxis','zmaxis','simag','sibry','bcentr',
        'current','simag','xdum','rmaxis','xdum',
        'zmaxis','xdum','sibry','xdum','xdum']
     
        #read in all 0D floats
        for key in float_keys:                              
            input_data[key] = np.asarray(float(processing.utils.get_next(token))) #processing.utils.get_next(token) always yields just a single value, convert this to numpy array
     
        def read_1d(n): 
            """
            """
            input_data = np.zeros(n) #initialise blank lists
            for i in np.arange(n): #instead of using linspace or something makes a temporary numpy array of dimension n to iterate through
                input_data[i] = float(processing.utils.get_next(token))
            return input_data
     
        def read_2d(nx,ny):
            """
            notes:
                edited this function as it was not consistent with GEQDSK storage format
            """
            input_data = np.zeros((nx,ny))
            for j in np.arange(ny): #same technique here, iterate through one dimension and call read_1d along the other dimension
                input_data[:,j] = read_1d(nx)
            return input_data

        #read in the arrays
        input_data['fpol'] = read_1d(input_data['nR_1D']) #remember input_data['nR_1D'] holds width, input_data['nZ_1D'] holds height
        input_data['pres'] = read_1d(input_data['nR_1D'])
        input_data['ffprime'] = read_1d(input_data['nR_1D'])
        input_data['pprime'] = read_1d(input_data['nR_1D'])
        input_data['psirz'] = read_2d(input_data['nR_1D'],input_data['nZ_1D'])
        input_data['qpsi'] = read_1d(input_data['nR_1D'])
     
        #now deal with boundaries
        input_data['lcfs_n'] = np.array(int(processing.utils.get_next(token)))
        input_data['limitr'] = np.array(int(processing.utils.get_next(token)))
     
        def read_bndy(nb,nl): #number of boundaries and limiters
            if nb > 0: #read in boundaries
                rb = np.zeros(nb)
                zb = np.zeros(nb)
                for i in np.arange(nb):
                    rb[i] = float(processing.utils.get_next(token)) #read in R,Z pairs
                    zb[i] = float(processing.utils.get_next(token))
            else:
                rb = np.array(0)
                zb = np.array(0)
         
            if nl > 0: #read in limiters
                rl = np.zeros(nl)
                zl = np.zeros(nl)
                for i in np.arange(nl):
                    rl[i] = float(processing.utils.get_next(token))
                    zl[i] = float(processing.utils.get_next(token))
            else:
                rl = np.array(0)
                zl = np.array(0)
     
            return rb,zb,rl,zl
     
        input_data['lcfs_r'],input_data['lcfs_z'],input_data['rlim'],input_data['zlim'] = read_bndy(input_data['lcfs_n'],input_data['limitr'])

        GEQDSKFIX_factor=1.
        if 'GEQDSKFIX1' in properties and properties['GEQDSKFIX1'] is True: #apply LOCUST flag transformations
            input_data['psirz']*=-1.
            input_data['sibry']*=-1.
            input_data['simag']*=-1.
            input_data['current']*=-1. #not done in LOCUST
            GEQDSKFIX_factor*=-1.
        if 'GEQDSKFIX2' in properties and properties['GEQDSKFIX2'] is True: 
            input_data['fpol']*=-1.
            input_data['bcentr']*=-1. #not done in LOCUST
        PSI_sclh=1./(input_data['sibry']-input_data['simag'])
        IPDIRh=-1. if PSI_sclh > 0. else 1. 
        if 'BPFLIP' in properties and properties['BPFLIP'] is True:
            IPDIRh*=-1.
        ITDIRh=-1. if input_data['fpol'][0]<0. else 1.
        if 'BTFLIP' in properties and properties['BTFLIP'] is True:
            ITDIRh*=-1.

        #additional data
        input_data['R_1D']=np.linspace(input_data['rleft'],input_data['rleft']+input_data['rdim'],num=input_data['nR_1D'])     
        input_data['Z_1D']=np.linspace(input_data['zmid']-0.5*input_data['zdim'],input_data['zmid']+0.5*input_data['zdim'],num=input_data['nZ_1D']) 
        input_data['flux_pol']=np.linspace(input_data['simag'],input_data['sibry'],input_data['ffprime'].size) #all 1D profiles are defined against a flux grid, so use any 1D profile's length
        input_data['flux_tor']=processing.process_input.QTP_calc(Q=input_data['qpsi'],P=input_data['flux_pol'])*GEQDSKFIX_factor #if we flipped poloidal flux this will also flip toroidal flux since we assume Q sign always +ve by convention

        print("finished reading equilibrium from GEQDSK")

    return input_data

def read_equilibrium_IDS(shot,run,**properties): 
    """
    reads relevant LOCUST equilibrium data from an equilibrium IDS and returns as a dictionary
 
    notes:
        idum not read
        assumes dim1, dim2 of IDS are R,Z respectively
        reads limiter and LCFS from same source
    """
 
    print("reading equilibrium from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_equilibrium_IDS could not import IMAS module!\nreturning\n")
        return

    input_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    if 'username' not in properties: properties['username']=settings.username
    if 'imasdb' not in properties: properties['imasdb']=settings.imasdb
    if 'imas_version' not in properties: properties['imas_version']=settings.imas_version
    input_IDS.open_env(properties['username'],properties['imasdb'],properties['imas_version'])
    input_IDS.equilibrium.get() #open the file and get all the data from it
 
    input_data = {} #initialise blank dictionary to hold the data
 
    #easy bits
    #0D data
    input_data['rcentr']=np.asarray(input_IDS.equilibrium.vacuum_toroidal_field.r0).reshape([]) #need reshape to ensure shape consistency
    input_data['rmaxis']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.r).reshape([])
    input_data['zmaxis']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.z).reshape([])
    input_data['simag']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.psi_axis/(2.0*constants.pi)).reshape([]) #convert to Wb/rad
    input_data['sibry']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.psi_boundary/(2.0*constants.pi)).reshape([]) #convert to Wb/rad
    input_data['current']=np.abs(np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.ip).reshape([])) #force +ve as Yueqiang does (LOCUST does not care)
 
    #1D data
    input_data['bcentr']=np.abs(np.asarray(input_IDS.equilibrium.vacuum_toroidal_field.b0)) #force +ve as Yueqiang does (LOCUST does not care)
    input_data['fpol']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.f) #flux grid data
    input_data['pres']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.pressure)
    input_data['ffprime']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.f_df_dpsi)
    input_data['pprime']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.dpressure_dpsi)
    input_data['qpsi']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.q)
    input_data['rlim']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.lcfs.r) #boundaries, lcfs is obsolete in latest IMAS
    input_data['zlim']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.lcfs.z)
    input_data['lcfs_r']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.outline.r) 
    input_data['lcfs_z']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.outline.z)
    input_data['flux_pol']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.psi/(2.0*constants.pi)) #convert to Wb/rad
    input_data['flux_tor']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.phi/(2.0*constants.pi)) #convert to Wb/rad
    input_data['flux_tor_coord']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.rho_tor)
    R_1D=input_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim1 #dim1=R values,dim2=Z values
    Z_1D=input_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim2
    input_data['R_1D']=np.asarray(R_1D)
    input_data['Z_1D']=np.asarray(Z_1D)

    #2D data    
    input_data['psirz']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_2d[0].psi)/(2.0*constants.pi) #convert to Wb/rad
    input_data['phirz']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_2d[0].phi)/(2.0*constants.pi) #convert to Wb/rad
 
    #harder bits (values derived from grids and profiles)
    input_data['lcfs_n']=np.asarray(len(input_IDS.equilibrium.time_slice[0].boundary.outline.z)).reshape([])
    input_data['limitr']=np.asarray(len(input_IDS.equilibrium.time_slice[0].boundary.outline.z)).reshape([])
    input_data['nbbbs']=np.asarray(len(input_IDS.equilibrium.time_slice[0].boundary.lcfs.z)).reshape([])
    input_data['nR_1D']=np.asarray(len(R_1D)).reshape([])
    input_data['nZ_1D']=np.asarray(len(Z_1D)).reshape([])
    input_data['rleft']=np.asarray(min(R_1D)).reshape([])
    input_data['rdim']=np.asarray(abs(max(R_1D)-min(R_1D))).reshape([])
    input_data['zdim']=np.asarray(abs(max(Z_1D)-min(Z_1D))).reshape([])
    input_data['zmid']=np.asarray(0.5*(max(Z_1D)+min(Z_1D))).reshape([])
 
    input_IDS.close()
    print("finished reading equilibrium from IDS")

    return input_data

def read_equilibrium_UDA(shot,time,**properties):
    """
    reads equilibrium from analysed EFIT signals from MAST UDA database

    notes:

    args:
        shot - target shot number to read
        time - fetch equilibrium at this time
    """

    print("reading equilibrium from UDA")

    input_data={}

    try:
        import pyuda
        udaClient=pyuda.Client()
        getdata=udaClient.get
    except:
        raise ImportError("ERROR: read_equilibrium_UDA could not import pyuda!\nreturning\n")

    #time dependent data first
    input_data['psirz']=getdata('efm_psi(r,z)',shot)
    input_data['simag']=getdata('efm_psi_axis',shot)
    input_data['sibry']=getdata('efm_psi_boundary',shot)
    input_data['current']=getdata('efm_plasma_curr(C)',shot)
    input_data['bcentr']=getdata('efm_bvac_val',shot)
    input_data['rcentr']=getdata('efm_bvac_r',shot)
    input_data['pres']=getdata('efm_p(psi)_(c)',shot)
    input_data['fpol']=getdata('efm_f(psi)_(c)',shot)
    input_data['rmaxis']=getdata('efm_magnetic_axis_r',shot)
    input_data['zmaxis']=getdata('efm_magnetic_axis_z',shot)
    input_data['ffprime']=getdata('efm_ffprime',shot)
    input_data['pprime']=getdata('efm_pprime',shot)
    input_data['qpsi']=getdata('efm_q(psi)_(c)',shot)
    input_data['lcfs_n']=getdata('efm_lcfs(n)_(c)',shot)
    input_data['lcfs_r']=getdata('efm_lcfs(r)_(c)',shot)
    input_data['lcfs_z']=getdata('efm_lcfs(z)_(c)',shot)
    input_data['limitr']=getdata('efm_limiter(n)',shot) #limiter data not on time base but contains single time coordinate, so read in here
    input_data['rlim']=getdata('efm_limiter(r)',shot)
    input_data['zlim']=getdata('efm_limiter(z)',shot)

    #redefine read time based on available converged EFIT equilibria constructions 
    time_index=np.abs(input_data['psirz'].time.data-time).argmin()
    time_new=input_data['psirz'].time.data[time_index]

    #go through time dependent data and pick closest times
    for key in input_data:
        #print(key) #to check times
        time_index=np.abs(input_data[key].time.data-time_new).argmin()
        #print(input_data[key].time.data[time_index]) #to check times
        input_data[key]=np.array(input_data[key].data[time_index])

    input_data['lcfs_r']=input_data['lcfs_r'][0:input_data['lcfs_n']] #LCFS data seems to have junk on the end, so crop accordingly
    input_data['lcfs_z']=input_data['lcfs_z'][0:input_data['lcfs_n']]

    input_data['psirz']=input_data['psirz'].T

    input_data['R_1D']=np.asarray(getdata('efm_grid(r)',shot).data[0,:])
    input_data['Z_1D']=np.asarray(getdata('efm_grid(z)',shot).data[0,:])
    input_data['nR_1D']=np.asarray(getdata('efm_nw',shot).data)
    input_data['nZ_1D']=np.asarray(getdata('efm_nh',shot).data)
    input_data['rdim']=input_data['R_1D'][-1]-input_data['R_1D'][0]
    input_data['rleft']=input_data['R_1D'][0]
    input_data['zdim']=input_data['Z_1D'][-1]-input_data['Z_1D'][0]
    input_data['zmid']=np.asarray((input_data['Z_1D'][-1]+input_data['Z_1D'][0])*.5)
    input_data['flux_pol']=np.linspace(input_data['simag'],input_data['sibry'],100)
    
    '''
    input_data['rlim'] = np.array([
    0.195, 0.195, 0.280, 0.280, 0.280, 0.175,
    0.175, 0.190, 0.190, 0.330, 0.330, 0.535,
    0.535, 0.755, 0.755, 0.755, 1.110, 1.655,
    1.655, 1.655, 2.000, 2.000, 2.000, 2.000,
    2.000, 1.655, 1.655, 1.655, 1.110, 0.755,
    0.755, 0.755, 0.535, 0.535, 0.330, 0.330,
    0.190, 0.190, 0.175, 0.175, 0.280, 0.280,
    0.280, 0.195, 0.195])

    input_data['zlim'] = np.array([
    0.000, 1.080, 1.220, 1.450, 1.670, 1.670,
    1.720, 1.820, 1.905, 2.000, 2.000, 2.000,
    2.000, 2.000, 1.975, 1.826, 1.826, 1.826,
    1.975, 2.000, 2.000, 0.300, 0.000,-0.300,
    -2.000,-2.000,-1.975,-1.826,-1.826,-1.826,
    -1.975,-2.000,-2.000,-2.000,-2.000,-2.000,
    -1.905,-1.820,-1.720,-1.670,-1.670,-1.450,
    -1.220,-1.080, 0.000])'''

    print("finished reading equilibrium from UDA")

    return input_data

################################################################## Equilibrium write functions

def dump_equilibrium_GEQDSK(output_data,filepath,**properties):
    """
    generic function for writing GEQDSK-formatted data to file
 
    notes:
        originally written by Ben Dudson and edited by Nick Walkden
    """ 

    try:
        import time
        import itertools
    except:
        raise ImportError("ERROR: dump_equilibrium_GEQDSK could not import necessary modules!\nreturning\n")
        return

    print("writing equilibrium to GEQDSK")

    def write_number(file,number,counter):
        
        if number==0:
            separator = " "
            number = np.abs(number)
        elif number < 0:
            separator = "-"
            number = np.abs(number)
        else:
            separator = " "
        if processing.utils.get_next(counter) == 4:
            last = "\n"
        else:
            last = ""
         
        string = '{:.9e}'.format(float(number))
        #mant,exp = string.split('E')
        file.write(separator+string+last)
 
    def write_1d(file,array,counter):
        for num in array:
            write_number(file,num,counter)
 
    def write_2d(file,array,counter):
        ny = array.shape[1]
        for j in np.arange(ny):
            write_1d(file,array[:,j],counter)
     
    def write_bndry(file,R,Z,counter):
        for i in np.arange(len(list(R))):
            write_number(file,R[i],counter)
            write_number(file,Z[i],counter)
    
    with open(filepath,'w') as file:
        
        cnt = itertools.cycle([0,1,2,3,4]) #initialise counter

        EFIT_shot=19113 #just 'make up' a shot number and time (in ms) for now
        EFIT_time=23
        output_data['xdum']=np.array(0)
        line = "LOCUSTIO   "+time.strftime("%d/%m/%y")+"      #"+processing.utils.fortran_string(EFIT_shot,6)+processing.utils.fortran_string(EFIT_time,6)+processing.utils.fortran_string(int(output_data['xdum']),14)+processing.utils.fortran_string(int(output_data['nR_1D']),4)+processing.utils.fortran_string(int(output_data['nZ_1D']),4)+"\n"
        file.write(line)
 
        float_keys = [
        'rdim','zdim','rcentr','rleft','zmid',
        'rmaxis','zmaxis','simag','sibry','bcentr',
        'current','simag','xdum','rmaxis','xdum',
        'zmaxis','xdum','sibry','xdum','xdum']
        for key in float_keys:
            write_number(file,output_data[key],cnt)
    
        #need to interpolate onto equally spaced grid according to GEQDSK convention
        flux_pol_norm_old=(output_data['flux_pol']-output_data['flux_pol'][0])/(output_data['flux_pol'][-1]-output_data['flux_pol'][0])
        flux_pol_norm_new=np.linspace(0,1,output_data['nR_1D'])
        interpolated_variables={}
        for var in ['fpol','pres','ffprime','pprime','qpsi']:
            interpolator=processing.utils.interpolate_1D(flux_pol_norm_old,output_data[var])
            interpolated_variables[var]=interpolator(flux_pol_norm_new)

        write_1d(file,interpolated_variables['fpol'],cnt)
        if processing.utils.get_next(cnt) != 0: #if last character was not newline, write newline
            file.write('\n') 
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_1d(file,interpolated_variables['pres'],cnt)
        if processing.utils.get_next(cnt) != 0: #if last character was not newline, write newline
            file.write('\n') 
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_1d(file,-interpolated_variables['ffprime'],cnt) #Yueqiang reverses this (possibly ad hoc?) but LOCUST does not use this data
        if processing.utils.get_next(cnt) != 0: #if last character was not newline, write newline
            file.write('\n') 
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_1d(file,-interpolated_variables['pprime'],cnt) #Yueqiang reverses this (possibly ad hoc?) but LOCUST does not use this data
        if processing.utils.get_next(cnt) != 0: #if last character was not newline, write newline
            file.write('\n') 
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_2d(file,output_data['psirz'],cnt)    
        if processing.utils.get_next(cnt) != 0: #if last character was not newline, write newline
            file.write('\n') 
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_1d(file,interpolated_variables['qpsi'],cnt) 
        if processing.utils.get_next(cnt) == 0: #check if \n was written last
            file.write(processing.utils.fortran_string(len(output_data['lcfs_r']),5)+processing.utils.fortran_string(len(output_data['rlim']),5)+'\n') #write out number of limiter/plasma boundary points
        else:
            file.write('\n'+processing.utils.fortran_string(len(output_data['lcfs_r']),5)+processing.utils.fortran_string(len(output_data['rlim']),5)+'\n') #write out number of limiter/plasma boundary points
        cnt = itertools.cycle([0,1,2,3,4]) #reinitialise counter
        write_bndry(file,output_data['lcfs_r'],output_data['lcfs_z'],cnt)
        write_bndry(file,output_data['rlim'],output_data['zlim'],cnt)
        file.write("\n") #file has a newline here

        print("finished writing equilibrium to GEQDSK")
 
def dump_equilibrium_IDS(ID,output_data,shot,run,**properties):
    """
    writes relevant LOCUST equilibrium data to an equilibrium IDS
 
    notes:
        current format in line with NEMO input
        currently only for rectangular equilibria 
        currently overwrites pre-existing IDSs
        idum not dumped
    """
 
    print("writing equilibrium to IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_equilibrium_IDS could not import IMAS module!\nreturning\n")
        return

    output_IDS=imas.ids(int(shot),int(run)) 
    if 'username' not in properties: properties['username']=settings.username
    if 'imasdb' not in properties: properties['imasdb']=settings.imasdb
    if 'imas_version' not in properties: properties['imas_version']=settings.imas_version
    output_IDS.open_env(settings.username,settings.imasdb,'3') #open the IDS
    output_IDS.equilibrium.get()
 
    #write out code properties
    output_IDS.equilibrium.ids_properties.comment=ID #write out identification
    output_IDS.equilibrium.code.name="LOCUST_IO"
    if settings.commit_hash_default_LOCUST_IO: output_IDS.equilibrium.code.commit=str(settings.commit_hash_default_LOCUST_IO)
    output_IDS.equilibrium.code.version=support.LOCUST_IO_version
    output_IDS.equilibrium.ids_properties.homogeneous_time=1   #must set homogeneous_time variable
    output_IDS.equilibrium.time=np.array([0.0])

    #add a time_slice and set the time of this slice
    if len(output_IDS.equilibrium.time_slice)==0: output_IDS.equilibrium.time_slice.resize(1) #just add one time_slice i.e. static equilibrium
    output_IDS.equilibrium.time_slice[0].time=0.0
    output_IDS.equilibrium.time=np.array(0.0,ndmin=1) #set the global time (required by vacuum_toroidal_field.b0)
 
    #write out the easy stuff - global quantities, some 1D profiles and the boundaries
    output_IDS.equilibrium.vacuum_toroidal_field.r0=output_data['rcentr'].item(0) 
    output_IDS.equilibrium.vacuum_toroidal_field.b0=np.array(output_data['bcentr'],ndmin=1) #this needs to be 1D and match the dimensions of output_IDS.equilibrium.time (above)
    output_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.r=output_data['rmaxis'].item(0)   
    output_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.z=output_data['zmaxis'].item(0)    
    output_IDS.equilibrium.time_slice[0].global_quantities.psi_axis=output_data['simag'].item(0)  
    output_IDS.equilibrium.time_slice[0].global_quantities.psi_boundary=output_data['sibry'].item(0)    
    output_IDS.equilibrium.time_slice[0].global_quantities.ip=output_data['current'].item(0)
 
    output_IDS.equilibrium.time_slice[0].boundary.type=0 #boundary type (0 for limiter, 1 for diverted)
    output_IDS.equilibrium.time_slice[0].boundary.outline.r=output_data['lcfs_r'] 
    output_IDS.equilibrium.time_slice[0].boundary.outline.z=output_data['lcfs_z']
    output_IDS.equilibrium.time_slice[0].boundary.lcfs.r=output_data['rlim'] #NOTE this is apparently obsolete - need to figure out where to write to 
    output_IDS.equilibrium.time_slice[0].boundary.lcfs.z=output_data['zlim']
     
    #write out the uniform flux grid output_data
    output_IDS.equilibrium.time_slice[0].profiles_1d.f=output_data['fpol'] 
    output_IDS.equilibrium.time_slice[0].profiles_1d.pressure=output_data['pres'] 
    output_IDS.equilibrium.time_slice[0].profiles_1d.f_df_dpsi=output_data['ffprime'] 
    output_IDS.equilibrium.time_slice[0].profiles_1d.dpressure_dpsi=output_data['pprime'] 
    output_IDS.equilibrium.time_slice[0].profiles_1d.q=output_data['qpsi'] 
    output_IDS.equilibrium.time_slice[0].profiles_1d.psi=output_data['flux_pol']
    output_IDS.equilibrium.time_slice[0].profiles_1d.phi=output_data['flux_tor']
 
    #now define the R,Z simulation grid
    output_IDS.equilibrium.time_slice[0].profiles_2d.resize(1) #add an element onto the profiles_2d struct_array to define this grid
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid_type.name='rectangular grid' #add some identifiers for this particular grid
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid_type.description=''
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid_type.index=1 #1 for rectangular (R,Z)
  
    #write out R,Z grid coordinate arrays
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim1=output_data['R_1D'] #dim1=R values/dim2=Z values
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim2=output_data['Z_1D']
    Z_2D,R_2D=np.meshgrid(output_data['Z_1D'],output_data['R_1D']) #generate 2D arrays of R,Z values
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].r=R_2D
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].z=Z_2D
     
    #write out 2D profiles
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].psi=output_data['psirz'] 

    #write out magnetic field if present
    #XXX this combination is deprecated but is compatible with NEMO
    if 'B_field_R' in output_data: output_IDS.equilibrium.time_slice[0].profiles_2d[0].b_r=output_data['B_field_R']
    if 'B_field_tor' in output_data: output_IDS.equilibrium.time_slice[0].profiles_2d[0].b_field_tor=output_data['B_field_tor']
    if 'B_field_Z' in output_data: output_IDS.equilibrium.time_slice[0].profiles_2d[0].b_z=output_data['B_field_Z']
    
    #calculate and write out 2D toroidal flux profile 
    if 'flux_pol' in output_data and 'flux_tor' in output_data:
        output_IDS.equilibrium.time_slice[0].profiles_2d[0].phi=2.*np.pi*processing.utils.flux_func_to_RZ(output_data['flux_pol'],output_data['flux_tor'],output_data) 

    #'put' all the output_data into the file and close
    output_IDS.equilibrium.put()
    output_IDS.close()

    print("finished writing equilibrium to IDS")

def dump_equilibrium_ASCOT(output_data,filepath,**properties):
    """
    writes equilibrium data to ASCOT input files

    notes:
    """

    dump_equilibrium_GEQDSK(output_data,filepath,**properties) #ASCOT also requires GEQDSK so write out just in case

    #XXX under construction

    '''eqd = read_eqdsk_riplos1('g157418.03000');
                eqd.psirz = -eqd.psirz;
            
            
                [bkg,Xguess] = eqdsk2magn_bkg(eqd, [], [], [], [], [1.3 -1.1], true);
            
            
            
                head = eqdsk2ascbkg(eqd, [], [], [], [], [], Xguess);
            
            
            
                write_magn_bkg(bkg, 'input.magn_bkg')
            
            
            
                write_magn_header(head, 'input.magn_header')'''

################################################################## Equilibrium class
 
class Equilibrium(classes.base_input.LOCUST_input):
    """
    class describing the equilibrium input for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'equilibrium'
    class data
        self.data_format            data format of original data e.g. GEQDSK
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
        MACHINE - GEQDSKFIX
            MAST    0
            DIII-D  1
            ITER    
            AUG     0
    """
 
    LOCUST_input_type='equilibrium'
  
    def read_data(self,data_format=None,filename=None,shot=None,run=None,time=None,**properties): #always supply all possible arguments for reading in data, irrespective of read in type
        """
        read equilibrium from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='GEQDSK': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from GEQDSK - filename required\n".format(self.ID),filename): #check we have all info for reading GEQDSKs
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_equilibrium_GEQDSK(self.filepath,**properties)
            
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from equilibrium IDS - shot and run required\n".format(self.ID),shot,run):
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_equilibrium_IDS(self.shot,self.run,**properties)

        elif data_format=='UDA':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from UDA - shot and time required\n".format(self.ID),shot,time):
                self.data_format=data_format
                self.shot=shot
                self.time=time
                self.properties={**properties}
                self.data=read_equilibrium_UDA(self.shot,self.time,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (GEQDSK/IDS/UDA)\n".format(self.ID))
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write equilibrium to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID = {})".format(self.ID)) 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"dump_data requires self.data and data_format\n",self.data,data_format):
            pass
         
        elif data_format=='GEQDSK':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to GEQDSK - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_equilibrium_GEQDSK(self.data,filepath,**properties)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to equilibrium IDS - shot and run required\n".format(self.ID),shot,run):
                dump_equilibrium_IDS(self.ID,self.data,shot,run,**properties)
 
        elif data_format=='ASCOT':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to ASCOT - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_equilibrium_ASCOT(self.data,filepath,**properties)

        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (GEQDSK/IDS/ASCOT)\n".format(self.ID))

    def plot(self,key='psirz',LCFS=True,limiters=False,number_bins=20,fill=True,vminmax=None,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,label='',ax=False,fig=False):
        """
        plots equilibrium
        
        notes:
            
        args:
            key - selects which data in equilibrium to plot
            LCFS - toggles plasma boundary on/off in 2D plots
            limiters - toggles limiters on/off in 2D plots
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            vminmax - set mesh Vmin/Vmax values
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            line_style - set 1D line style
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        number_bins,colmap_val,vminmax=run_scripts.utils.literal_eval(number_bins,colmap_val,vminmax)

        import scipy
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

        #0D data
        if self[key].ndim==0:
            print(self[key])
            return
        
        #>0D data is plottable
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)
        ax.set_title(self.ID)

        #1D data
        if self[key].ndim==1:
            ax.plot(self[key],color=colmap(colmap_val),linestyle=line_style,label=label)
            ax.set_ylabel(key)

            
                

        #2D data
        elif self[key].ndim==2:

            X=self['R_1D'] #make a mesh
            Y=self['Z_1D']
            dx,dy=X[1]-X[0],Y[1]-Y[0]

            Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres
            Z=self[key] #2D array (nR_1D,nZ_1D) of poloidal flux

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.amin(Z)
                vmax=np.amax(Z)
            
            #2D plot
            if fill is True:
                mesh=ax.contourf(X,Y,Z,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax)
                #mesh=ax.pcolormesh(X,Y,Z,edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax,cmap=colmap)
                for c in mesh.collections: #for use in contourf
                    c.set_edgecolor("face")
            else:
                mesh=ax.contour(X,Y,Z,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax)
                if settings.plot_contour_labels:
                    ax.clabel(mesh,inline=1,fontsize=10)
                
            #mesh=ax.pcolormesh(X,Y,Z,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))

            #3D plot
            #ax=ax.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
            
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')
            ax.set_aspect('equal')
            ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
            ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D']))
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

            if LCFS:
                ax.plot(self['lcfs_r'],self['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
            if limiters: #add boundaries if desired
                ax.plot(self['rlim'],self['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall') 

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()


    def plot_field_line(self,axes=['X','Y','Z'],LCFS=False,limiters=False,number_field_lines=1,angle=2.0*constants.pi,plot_full=False,start_mark=True,start_coord=None,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,label='',ax=False,fig=False):
        """
        plots random field lines for 'angle' radians around the tokamak

        notes:
            essentially uses the Euler method of integration
        args:
            axes - list of strings specifying which axes should be plotted
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            number_field_lines - the number of field lines to plot
            angle - plot field line for this many radians around central column
            plot_full - choose whether each field line will trace the whole plasma topology (see below also)
            start_mark - add marker for start of the field line
            start_coord - optional choose [R,phi,Z] starting position 
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots
            line_style - set 1D line style
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        import scipy
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

        if not np.all([component in self.data.keys() for component in ['B_field_R','B_field_tor','B_field_Z']]):
            print("plot_field_line - found no B_field in equilibrium - calculating!")
            self.B_calc()

        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)
        ax.set_title(self.ID)
        ax.set_aspect('equal')
        
        dr=np.abs(self['R_1D'][1]-self['R_1D'][0])
        dz=np.abs(self['Z_1D'][1]-self['Z_1D'][0])
        
        if plot_full is True: #set integration path step size
            #if this option, then dl chosen to cause slight numerical drift that such that field line strays outward onto
            #different flux surfaces - thus tracing the whole field topology
            dl=3.0*np.sqrt(dr**2+dz**2) 
            angle=constants.pi*200.0
            number_field_lines=1
            ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
            ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D'])) 
        else:
            #if this option chosen, then dl is reduced to stop numerical drift - thus a single flux surface is plotted
            dl=0.005*np.sqrt(dr**2+dz**2)
            if ax_flag is False and len(axes)==3:
                ax = fig.gca(projection='3d')
            
        print('plot_B_field_line - generating B field interpolators')
        B_field_R_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_R'])
        B_field_tor_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_tor'])
        B_field_Z_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_Z'])
        print('plot_B_field_line - finished generating B field interpolators')

        for line in range(number_field_lines): 

            R_points=np.array([]) #reset the arrays that will hold marker coordinates along a single field line
            Z_points=np.array([])
            tor_points=np.array([])

            if not start_coord:
                R_point=float(self['rmaxis']) #pick some starting points                                                       
                tor_point=np.random.uniform(0.0,2.0*constants.pi*R_point) 
                if plot_full is True: 
                    Z_point=float(1.05*self['zmaxis'])
                else:
                    Z_point=np.random.uniform(1.05*self['zmaxis'],0.9*np.max(self['lcfs_z']))      
            else:
                R_point,tor_point,Z_point=[coord for coord in start_coord]

            R_points=np.append(R_points,R_point) #add this position to our array of points along trajectory
            Z_points=np.append(Z_points,Z_point)
            tor_points=np.append(tor_points,tor_point)     
            tor_point_start=tor_point #remember where we started toroidally

            while np.abs(tor_point-tor_point_start)<self['rmaxis']*angle: #keep going until we rotate 'angle' radians around the tokamak
                
                B_field_R=B_field_R_interpolator(R_point,Z_point)
                B_field_tor=B_field_tor_interpolator(R_point,Z_point)
                B_field_Z=B_field_Z_interpolator(R_point,Z_point)

                #could save computing power by just dividing by largest B component*some_constant - since we only need to divide the B vector by something we know is going to be larger than than the largest component
                #and remember (in 3D) the magnitude cannot ever be more than sqrt(3) times larger than the largest component of the field, so if we divide by at least sqrt(3)*np.max(field) then we always know our normalised values will be < 1.0
                B_field_mag=np.sqrt(B_field_R**2+B_field_Z**2+B_field_tor**2)
                
                B_field_R/=B_field_mag #normalise the vector magnetic field
                B_field_Z/=B_field_mag
                B_field_tor/=B_field_mag

                Z_point+=B_field_Z*dl
                tor_point+=B_field_tor*dl
                R_point+=B_field_R*dl 

                R_points=np.append(R_points,R_point)
                Z_points=np.append(Z_points,Z_point)
                tor_points=np.append(tor_points,tor_point)    

            R_points=R_points[1::2] #take every other value to help the visuals
            Z_points=Z_points[1::2]
            tor_points=tor_points[1::2]
            X_points=R_points*np.cos(tor_points/R_points) #transform to cartesian
            Y_points=R_points*np.sin(tor_points/R_points) 

            if plot_full is True: #if wanting to trace the flux surfaces, then plot in r,z plane
                ax.plot(R_points,Z_points,color=colmap(colmap_val),label=label)
                if start_mark: 
                    ax.scatter(R_points[0],Z_points[0],color=settings.colour_start_mark,s=10)
            else:
                
                if axes==['R','Z']: #poloidal plot
                    ax.plot(R_points,Z_points,color=colmap(colmap_val),label=label)
                    if start_mark: 
                        ax.scatter(R_points[0],Z_points[0],color=settings.colour_start_mark,s=10)
                    if LCFS: #plot plasma boundary
                        ax.plot(self['lcfs_r'],self['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
                    if limiters: #add boundaries if desired
                        ax.plot(self['rlim'],self['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')       
            
                    ax.set_xlabel('R [m]')
                    ax.set_ylabel('Z [m]')
                    ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
                    ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D']))

                elif axes==['X','Y']: #top-down plot
                    ax.plot(X_points,Y_points,color=colmap(colmap_val),label=label)
                    if start_mark: 
                        ax.scatter(X_points[0],Y_points[0],color=settings.colour_start_mark,s=10)
                    if LCFS: #plot plasma boundary
                        plasma_max_R=np.max(self['lcfs_r'])
                        plasma_min_R=np.min(self['lcfs_r'])
                        ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                        ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
                    if limiters: #add boundaries if desired
                        ax.set_xlim(-1.0*np.max(self['rlim']),np.max(self['rlim']))
                        ax.set_ylim(-1.0*np.max(self['rlim']),np.max(self['rlim']))
                        limiters_max_R=np.max(self['rlim'])
                        limiters_min_R=np.min(self['rlim'])
                        ax.plot(limiters_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')
                        ax.plot(limiters_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')   
                        
                    ax.set_xlabel('X [m]')
                    ax.set_ylabel('Y [m]')
                    ax.set_xlim(-1.0*np.max(self['R_1D']),np.max(self['R_1D']))
                    ax.set_ylim(-1.0*np.max(self['R_1D']),np.max(self['R_1D']))
                
                else: #3D plot
                    if start_mark: 
                        ax.scatter(X_points[0],Y_points[0],Z_points[0],color=settings.colour_start_mark,s=10)
                    if LCFS: #plot periodic poloidal cross-sections in 3D
                        for angle in np.linspace(0.0,2.0*constants.pi,4,endpoint=False):
                            x_points=self['lcfs_r']*np.cos(angle)
                            y_points=self['lcfs_r']*np.sin(angle)
                            z_points=self['lcfs_z']
                            ax.plot(x_points,y_points,zs=z_points,color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                    if limiters: #plot periodic poloidal cross-sections in 3D
                        for angle in np.linspace(0.0,2.0*constants.pi,4,endpoint=False):
                            x_points=self['rlim']*np.cos(angle)
                            y_points=self['rlim']*np.sin(angle)
                            z_points=self['zlim']
                            ax.plot(x_points,y_points,zs=z_points,color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')

                    ax.set_xlim(-1.0*np.max(self['R_1D']),np.max(self['R_1D']))
                    ax.set_ylim(-1.0*np.max(self['R_1D']),np.max(self['R_1D']))
                    ax.set_zlim(np.min(self['Z_1D']),np.max(self['Z_1D'])) 
                    ax.plot(X_points,Y_points,zs=Z_points,color=colmap(colmap_val),label=label)

        if ax_flag is False and fig_flag is False:
            plt.show()


    def plot_field_stream(self,LCFS=True,limiters=True,colmap=settings.cmap_default,line_style=settings.plot_line_style,label='',ax=False,fig=False):
        """
        stream plot of magnetic field in R,Z plane

        args:
            LCFS - toggles plasma boundary on/off in 2D plots
            limiters - toggles limiters on/off in 2D plots
            colmap - set the colour map (use get_cmap names)
            line_style - set 1D line style
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        notes:
            take transpose due to streamplot index convention
        """

        import scipy
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
            ax.set_title(self.ID)
        ax.set_aspect('equal')

        if not np.all([component in self.data.keys() for component in ['B_field_R','B_field_tor','B_field_Z']]): #calculate B field if missing
            print("plot_field_stream - found no B_field in equilibrium - calculating!")
            self.B_calc()

        B_mag=np.sqrt(self['B_field_R']**2+self['B_field_Z']**2) #calculate poloidal field magnitude
        strm = ax.streamplot(self['R_1D'],self['Z_1D'],self['B_field_R'].T,self['B_field_Z'].T, color=B_mag.T, linewidth=1, cmap=colmap)

        if LCFS:
            ax.plot(self['lcfs_r'],self['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
        if limiters: #add boundaries if desired
            ax.plot(self['rlim'],self['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall') 

        if fig_flag is False:    
            fig.colorbar(strm.lines,ax=ax,orientation='horizontal')
        ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
        ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D']))

        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')

        if ax_flag is False and fig_flag is False:
            plt.show()

    def fpolrz_calc(self):
        """
        interpolates 1D poloidal flux function fpol onto an r,z grid given  

        notes:
            fpol is defined over the poloidal flux points flux_pol
            assumes flux_pol and fpol are defined from magnetic axis to plasma boundary
            if the grid extends outside of the plasma radius, where poloidal flux function is not defined, then use a given vacuum toroidal field and work
            backwards to get the value of the flux function (since f is constant in vacuum)
        """

        print("fpolrz_calc - calculating 2D flux function")

        fpolrz=processing.utils.flux_func_to_RZ(self['flux_pol'],self['fpol'],self)
        Z,R=np.meshgrid(self['Z_1D'],self['R_1D'])
        within_LCFS=processing.utils.within_LCFS(R.flatten(),Z.flatten(),self).reshape(self['nR_1D'],self['nZ_1D'])

        for w in np.arange(self['nR_1D']): #loop over 2D grid and any points outside LCFS set to vacuum field
            for h in np.arange(self['nZ_1D']):
                if not within_LCFS[w,h]:
                    fpolrz[w,h]=self['bcentr']*self['rcentr']

        print("fpolrz_calc - finished calculating 2D flux function")

        self.set(fpolrz=fpolrz)

    def B_calc(self): #ordering OK - for Z[X,Y].shape=(2,3), gradient(Z,X,Y) --> gradient[0] along X, gradient[1] along Y
        """
        calculates r, phi, z axisymmetric magnetic field at coordinates R_1D, Z_1D

        notes:
            see Wesson, based on RH coordinate system [R,Phi,Z]
            calculates B_field_i(r,z) (3 2D arrays holding component i of B field at grid indices r, z)

            remember - A[i,j,k] is an array nested such as i*[j*[k*[some_number]]]
                     - [row,column] in numpy
        """
        
        if 'fpolrz' not in self.data: 
            print("WARNING: B_calc - fpolrz missing in equilibrium object - calculating!")
            self.fpolrz_calc() #if fpolrz is missing calculate it

        print("B_calc - calculating 2D magnetic field")

        gradient=np.gradient(self['psirz'],self['R_1D'],self['Z_1D']) #calculate gradient along both axes (gradient[i] is 2D)

        one_R=1.0/self['R_1D']
        one_R=one_R[:,np.newaxis]

        B_field_Z=(gradient[0])*one_R
        B_field_R=(gradient[1])*one_R*(-1.0)
        B_field_tor=self['fpolrz']*one_R

        self.set(B_field_R=B_field_R)
        self.set(B_field_tor=B_field_tor)
        self.set(B_field_Z=B_field_Z)

        print("B_calc - finished calculating 2D magnetic field")

    def mag_axis_calc(self,threshold=1.e-7):
        """
        args:
            threshold - final result is accurate to threshold*grid_spacing
        usage:
            eq['rmaxis'],eq['zmaxis'],eq['simag']=eq.mag_axis_calc()
        notes:
        """

        interpolator=processing.utils.interpolate_2D(        
                                        self['R_1D'],
                                        self['Z_1D'],
                                        self['psirz'])

        psirz_norm=np.abs((self['psirz']-self['simag'])/(self['sibry']-self['simag']))
        
        rmaxis,zmaxis,simag_norm=processing.utils.local_minimum_2D(
                                                    self['rmaxis'],
                                                    self['zmaxis'],
                                                    self['R_1D'],
                                                    self['Z_1D'],
                                                    psirz_norm,
                                                    threshold=threshold)
        simag=interpolator(self['rmaxis'],self['zmaxis'])[0][0]

        return rmaxis,zmaxis,simag

    def evaluate(self,R,Z):
        """
        returns the three components of magnetic field at a point in the plasma 
        
        args:
            R - array of R coordinates to calculate magnetic field components at 
            Z - array of Z coordinates to calculate magnetic field components at
        notes:

        usage:
            B_R,B_tor,B_Z=my_equilibrium.evaluate(R=[1,2,3],Z=[1,2,3])
        """
        
        if not np.all([component in self.data.keys() for component in ['B_field_R','B_field_tor','B_field_Z']]): #calculate B field if missing
            print("evaluate found no B_field in equilibrium - calculating!")
            self.B_calc()

        print("evaluate generating B_field interpolators")
        B_field_R_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_R']) #construct interpolators here
        B_field_tor_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_tor'])
        B_field_Z_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['B_field_Z'])
        print("evaluate finished generating B_field interpolators")

        B_R=[]
        B_tor=[]
        B_Z=[]  

        for R_point,Z_point in zip(R,Z):
            B_R.append(float(B_field_R_interpolator(R_point,Z_point)))
            B_tor.append(float(B_field_tor_interpolator(R_point,Z_point)))
            B_Z.append(float(B_field_Z_interpolator(R_point,Z_point)))
        
        B_R=np.asarray(B_R)
        B_tor=np.asarray(B_tor)
        B_Z=np.asarray(B_Z)

        return B_R,B_tor,B_Z

    def calc_r_LCFS(self):
        """
        calculate the minor radius at the outboard LCFS
        
        notes:

        """

        #first thing to do is accurately determine the magnetic axis location
        rmaxis_equilibrium,zmaxis_equilibrium,simag=self.mag_axis_calc(threshold=1.e-7)

        #create data along outboard midplane for interpolation 
        midplane_radius_major=np.linspace(rmaxis_equilibrium,self['R_1D'][-1],100)
        midplane_poloidal_flux=processing.utils.value_at_RZ(
        R=midplane_radius_major,
        Z=np.full(len(midplane_radius_major),
            zmaxis_equilibrium),
        quantity=self['psirz'],
        grid=self)
    
        interpolator_r_min_of_psi=processing.utils.interpolate_1D(midplane_poloidal_flux,midplane_radius_major-rmaxis_equilibrium,function='linear',type='interp1d') #create r_minor(psi) interpolator
        
        r_min_sibry=interpolator_r_min_of_psi(self['sibry'])

        return r_min_sibry

#################################
 
##################################################################
 
###################################################################################################
