#equilibrium.py
 
"""
Samuel Ward
02/11/2017
----
class to handle LOCUST equilibrium input data
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
    print("WARNING: IMAS module could not be imported!\n")
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/ could not be imported!\nreturning\n")
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
try:
    from processing import process_input 
except:
    raise ImportError("ERROR: process_input.py could not be imported!\nreturning\n")
    sys.exit(1)

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
 
pi=np.pi



################################################################## Equilibrium functions

'''
def dump_2Dwall_LOCUST(output_data,filepath):
    """
    dumps 2D wall profile to LOCUST format

    notes:
        output_data must be GEQDSK-style object
        2D wall points must be monotonic in poloidal angle
    """
    
    with open(filepath,'r') as file:
        #write the number of points (must be 3600 always)
        #write the major radius
        #write out the radii of limiter points as we rotate poloidally 3600 times anti-clockwise from outboard side 
'''

def read_equilibrium_GEQDSK(filepath): 
    """ 
    generic function for reading a G-EQDSK-formatted equilibrium file
 
    notes:
        originally written by Ben Dudson and edited by Nick Walkden
    """
 
    input_data = {}

    print("reading equilibrium from GEQDSK")
    
    with open(filepath,'r') as file: #open file
 
        line = file.readline() #first line should be case, id number and dimensions
        if not line:
            raise IOError("ERROR: read_equilibrium_GEQDSK() cannot read from "+filepath)
         
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
        token = processing.utils.file_numbers(file) #token now holds all lines in the file, containing all values. each processing.utils.get_next() call will grab the next number in a line (since processing.utils.get_next() returns the next value in the yield loop? check this plz)
         
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
        
        #additional data
        input_data['R_1D']=np.linspace(input_data['rleft'],input_data['rleft']+input_data['rdim'],num=input_data['nR_1D'])     
        input_data['Z_1D']=np.linspace(input_data['zmid']-0.5*input_data['zdim'],input_data['zmid']+0.5*input_data['zdim'],num=input_data['nZ_1D']) 
        input_data['flux_pol']=np.linspace(input_data['simag'],input_data['sibry'],input_data['ffprime'].size) #all 1D profiles are defined against a flux grid, so use any 1D profile's length
        input_data['flux_tor']=process_input.QTP_calc(Q=input_data['qpsi'],P=input_data['flux_pol'])

        print("finished reading equilibrium from GEQDSK")

    return input_data
     
def dump_equilibrium_GEQDSK(output_data,filepath):
    """
    generic function for writing G-EQDSK-formatted data to file
 
    notes:
        originally written by Ben Dudson and edited by Nick Walkden
    """ 
    def write_number(file,number,counter):
        
        if number==0:
            separator = "  "
            number = np.abs(number)
        elif number < 0:
            separator = " -"
            number = np.abs(number)
        else:
            separator = "  "
        if processing.utils.get_next(counter) == 4:
            last = "\n"
        else:
            last = ""
         
        string = '{:.8e}'.format(float(number))
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
     

    print("writing equilibrium to GEQDSK")
    
    with open(filepath,'w') as file:
        
        cnt = itertools.cycle([0,1,2,3,4]) #initialise counter

        EFIT_shot=19113 #just 'make up' a shot number and time (in ms) for now
        EFIT_time=23
        line = "LOCUSTIO   "+time.strftime("%d/%m/%y")+"      #"+processing.utils.fortran_string(EFIT_shot,6)+processing.utils.fortran_string(EFIT_time,6)+processing.utils.fortran_string(output_data['xdum'],14)+processing.utils.fortran_string(output_data['nR_1D'],4)+processing.utils.fortran_string(output_data['nZ_1D'],4)+"\n"
        file.write(line)
 
        output_data['xdum']=0.
        float_keys = [
        'rdim','zdim','rcentr','rleft','zmid',
        'rmaxis','zmaxis','simag','sibry','bcentr',
        'current','simag','xdum','rmaxis','xdum',
        'zmaxis','xdum','sibry','xdum','xdum']
        for key in float_keys:
            write_number(file,output_data[key],cnt)
 
        cnt = itertools.cycle([0,1,2,3,4]) #reset the counter (otherwise our newlines will be out of sync)
        write_1d(file,output_data['fpol'],cnt)
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_1d(file,output_data['pres'],cnt)
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again        
        write_1d(file,output_data['ffprime'],cnt)
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_1d(file,output_data['pprime'],cnt)
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_2d(file,output_data['psirz'],cnt)    
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_1d(file,output_data['qpsi'],cnt) 
        
        file.write("\n"+processing.utils.fortran_string(len(output_data['lcfs_r']),5)+processing.utils.fortran_string(len(output_data['rlim']),5)) #write out number of limiter/plasma boundary points
        
        file.write("\n")
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_bndry(file,output_data['lcfs_r'],output_data['lcfs_z'],cnt)
        
        file.write("\n") #file has a newline here
        cnt = itertools.cycle([0,1,2,3,4]) #reset again
        write_bndry(file,output_data['rlim'],output_data['zlim'],cnt)

        print("finished writing equilibrium to GEQDSK")
 
def read_equilibrium_IDS(shot,run): 
    """
    reads relevant LOCUST equilibrium data from an equilibrium IDS and returns as a dictionary
 
    notes:
        idum not read
        assumes dim1, dim2 of IDS are R,Z respectively
    """
 
    print("reading equilibrium from IDS")

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.open()
    input_IDS.equilibrium.get() #open the file and get all the data from it
 
    input_data = {} #initialise blank dictionary to hold the data
 
    #easy bits
    #0D data
    input_data['rcentr']=np.asarray(input_IDS.equilibrium.vacuum_toroidal_field.r0).reshape([]) #need reshape to ensure shape consistency
    input_data['bcentr']=np.asarray(input_IDS.equilibrium.vacuum_toroidal_field.b0).reshape([])
    input_data['rmaxis']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.r).reshape([])
    input_data['zmaxis']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.z).reshape([])
    input_data['simag']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.psi_axis/(2.0*pi)).reshape([]) #convert to Wb/rad
    input_data['sibry']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.psi_boundary/(2.0*pi)).reshape([]) #convert to Wb/rad
    input_data['current']=np.asarray(input_IDS.equilibrium.time_slice[0].global_quantities.ip).reshape([])
 
    #1D data
    input_data['fpol']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.f) #flux grid data
    input_data['pres']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.pressure)
    input_data['ffprime']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.f_df_dpsi)
    input_data['pprime']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.dpressure_dpsi)
    input_data['qpsi']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.q)
    input_data['rlim']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.lcfs.r) #boundaries, lcfs is obsolete in latest IMAS
    input_data['zlim']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.lcfs.z)
    input_data['lcfs_r']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.outline.r) 
    input_data['lcfs_z']=np.asarray(input_IDS.equilibrium.time_slice[0].boundary.outline.z)
    input_data['flux_pol']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.psi/(2.0*pi)) #convert to Wb/rad
    input_data['flux_tor']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_1d.phi/(2.0*pi)) #convert to Wb/rad
    R_1D=input_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim1 #dim1=R values,dim2=Z values
    Z_1D=input_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim2
    input_data['R_1D']=np.asarray(R_1D)
    input_data['Z_1D']=np.asarray(Z_1D)

    #2D data    
    input_data['psirz']=np.asarray(input_IDS.equilibrium.time_slice[0].profiles_2d[0].psi)/(2.0*pi) #convert to Wb/rad
 
    #harder bits (values derived from grids and profiles)
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
 
def dump_equilibrium_IDS(ID,output_data,shot,run):
    """
    writes relevant LOCUST equilibrium data to an equilibrium IDS
 
    notes:
        currently only for rectangular equilibria 
        currently overwrites pre-existing IDSs
        idum not dumped
 
    """
 
    print("writing equilibrium to IDS")

    output_IDS=imas.ids(shot,run) 
    output_IDS.create() #this will overwrite any existing IDS for this shot/run
 
    #write out code properties
    output_IDS.equilibrium.ids_properties.comment=ID #write out identification
    output_IDS.equilibrium.code.name="LOCUST_IO"
    output_IDS.equilibrium.code.version=support.LOCUST_IO_version
    output_IDS.equilibrium.ids_properties.homoegeneous_time=0   #must set homogeneous_time variable
 
    #add a time_slice and set the time of this slice
    output_IDS.equilibrium.time_slice.resize(1) #just add one time_slice i.e. static equilibrium
    output_IDS.equilibrium.time_slice[0].time=0.0
    output_IDS.equilibrium.time=np.array(0.0,ndmin=1) #set the global time (required by vacuum_toroidal_field.b0)
 
    #write out the easy stuff - global quantities, some 1D profiles and the boundaries
    output_IDS.equilibrium.vacuum_toroidal_field.r0=output_data['rcentr'] 
    output_IDS.equilibrium.vacuum_toroidal_field.b0=np.array(output_data['bcentr'],ndmin=1) #this needs to be 1D and match the dimensions of output_IDS.equilibrium.time (above)
    output_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.r=output_data['rmaxis']   
    output_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.z=output_data['zmaxis']    
    output_IDS.equilibrium.time_slice[0].global_quantities.psi_axis=output_data['simag']  
    output_IDS.equilibrium.time_slice[0].global_quantities.psi_boundary=output_data['sibry']    
    output_IDS.equilibrium.time_slice[0].global_quantities.ip=output_data['current']
 
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
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid_type.index=1 #1 for rectangular (R,Z), 0 for inverse (psi,theta)
  
    #write out R,Z grid coordinate arrays
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim1=output_data['R_1D'] #dim1=R values/dim2=Z values
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].grid.dim2=output_data['Z_1D']
    R_2D,Z_2D=np.meshgrid(output_data['R_1D'],output_data['Z_1D']) #generate 2D arrays of R,Z values
    R_2D,Z_2D=R_2D.T,Z_2D.T #since things are defined r,z need to take transpose here
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].r=R_2D
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].z=Z_2D
     
    #write out 2D profiles
    output_IDS.equilibrium.time_slice[0].profiles_2d[0].psi=output_data['psirz'] 
     
    #'put' all the output_data into the file and close
    output_IDS.equilibrium.put()
    output_IDS.close()

    print("finished writing equilibrium to IDS")

def dump_equilibrium_ASCOT(output_data,filepath):
    """
    writes equilibrium data to ASCOT input files

    notes:
    """

    dump_equilibrium_GEQDSK(output_data,filepath) #ASCOT also requires GEQDSK so write out just in case

    #XXX under construction

    '''eqd = read_eqdsk_riplos1('g157418.03000');
                eqd.psirz = -eqd.psirz;
            
            
                [bkg,Xguess] = eqdsk2magn_bkg(eqd, [], [], [], [], [1.3 -1.1], true);
            
            
            
                head = eqdsk2ascbkg(eqd, [], [], [], [], [], Xguess);
            
            
            
                write_magn_bkg(bkg, 'input.magn_bkg')
            
            
            
                write_magn_header(head, 'input.magn_header')'''


def read_equilibrium_UDA(shot,time):
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
        import numpy as np
    except:
        raise ImportError("ERROR: read_equilibrium_UDA could not import pyuda")

    #time dependent data first
    input_data['psirz']=getdata('efm_psi(r,z)',       shot)
    input_data['simag']=getdata('efm_psi_axis',       shot)
    input_data['sibry']=getdata('efm_psi_boundary',   shot)
    input_data['current']=getdata('efm_plasma_curr(C)', shot)
    input_data['bcentr']=getdata('efm_bvac_val',       shot)
    input_data['rcentr']=getdata('efm_bvac_r',       shot)
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

    #redefine read time based on available converged EFIT equilibria constructions 
    time_index=np.abs(input_data['psirz'].time.data-time).argmin()
    time_new=input_data['psirz'].time.data[time_index]

    #go through time dependent data and pick closest times
    for key in input_data:
        #print(key) #to check times
        time_index=np.abs(input_data[key].time.data-time_new).argmin()
        #print(input_data[key].time.data[time_index]) #to check times
        input_data[key]=np.asarray(input_data[key].data[time_index])

    input_data['R_1D']=np.asarray(getdata('efm_grid(r)',shot).data[0,:])
    input_data['Z_1D']=np.asarray(getdata('efm_grid(z)',shot).data[0,:])
    input_data['nR_1D']=np.asarray(getdata('efm_nw',shot))
    input_data['nZ_1D']=np.asarray(getdata('efm_nh',shot))
    input_data['rdim']=input_data['R_1D'][-1]-input_data['R_1D'][0]
    input_data['rleft']=input_data['R_1D'][0]
    input_data['zdim']=input_data['Z_1D'][-1]-input_data['Z_1D'][0]
    input_data['zmid']=(input_data['Z_1D'][-1]+input_data['Z_1D'][0])*.5

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
    -1.220,-1.080, 0.000])

    print("finished reading equilibrium from UDA")

    return input_data

################################################################## Equilibrium class
 
class Equilibrium(base_input.LOCUST_input):
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
 
    """
 
    LOCUST_input_type='equilibrium'
  
    def read_data(self,data_format=None,filename=None,shot=None,run=None,time=None,**properties): #always supply all possible arguments for reading in data, irrespective of read in type
        """
        read equilibrium from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='GEQDSK': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from GEQDSK - filename required\n",filename): #check we have all info for reading GEQDSKs
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_equilibrium_GEQDSK(self.filepath) #read the file
            
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from equilibrium IDS - shot and run required\n",shot,run):
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_equilibrium_IDS(self.shot,self.run)

        elif data_format=='UDA':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from UDA - shot and time required\n",shot,time):
                self.data_format=data_format
                self.shot=shot
                self.time=time
                self.properties={**properties}
                self.data=read_equilibrium_UDA(self.shot,self.time)

        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (GEQDSK/IDS/UDA)\n")
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write equilibrium to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"dump_data requires self.data and data_format\n",self.data,data_format):
            pass
         
        elif data_format=='GEQDSK':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to GEQDSK - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_equilibrium_GEQDSK(self.data,filepath)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to equilibrium IDS - shot and run required\n",shot,run):
                dump_equilibrium_IDS(self.ID,self.data,shot,run)
 
        elif data_format=='ASCOT':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_equilibrium_ASCOT(self.data,filepath)

        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (GEQDSK/IDS/ASCOT)\n")

 
#################################
 
##################################################################
 
###################################################################################################