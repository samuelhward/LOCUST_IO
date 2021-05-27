#beam_deposition.py
 
'''
Samuel Ward
02/11/2017
----
class to handle LOCUST beam deposition input data
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
    import copy
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


################################################################## Beam_Deposition read functions
 
def read_beam_depo_LOCUST_full_orbit(filepath,**properties):
    """
    reads birth profile stored in LOCUST format - R phi Z V_R V_phi V_Z

    notes:
        calculates energy in eV
    """
 
    print("reading full orbit beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+str(filepath))
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['absorption_fraction']=np.array(float(lines[0]))
        input_data['absorption_scaling']=np.array(float(lines[1]))
        del(lines[0])
        del(lines[0])
     
        input_data['R']=[] #initialise the arrays 
        input_data['phi']=[]
        input_data['Z']=[]
        input_data['V_R']=[]
        input_data['V_phi']=[]
        input_data['V_Z']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['R'].append(float(split_line[0]))
            input_data['phi'].append(float(split_line[1]))
            input_data['Z'].append(float(split_line[2]))
            input_data['V_R'].append(float(split_line[3]))
            input_data['V_phi'].append(float(split_line[4]))
            input_data['V_Z'].append(float(split_line[5]))
     
        input_data['R']=np.asarray(input_data['R']) #convert to arrays
        input_data['phi']=np.asarray(input_data['phi'])
        input_data['Z']=np.asarray(input_data['Z'])
        input_data['V_R']=np.asarray(input_data['V_R'])
        input_data['V_phi']=np.asarray(input_data['V_phi'])
        input_data['V_Z']=np.asarray(input_data['V_Z'])
        input_data['E']=np.asarray((input_data['V_R']**2+input_data['V_phi']**2+input_data['V_Z']**2)*.5*constants.species_mass/constants.species_charge)
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading full orbit beam deposition from LOCUST")
 
    return input_data

def read_beam_depo_LOCUST_full_orbit_weighted(filepath,**properties):
    """
    reads birth profile stored in LOCUST format - R phi Z V_R V_phi V_Z weight

    notes:
        calculates energy in eV
    """
 
    print("reading weighted full orbit beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+str(filepath))
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['absorption_fraction']=np.array(float(lines[0]))
        input_data['absorption_scaling']=np.array(float(lines[1]))
        del(lines[0])
        del(lines[0])
     
        input_data['R']=[] #initialise the arrays 
        input_data['phi']=[]
        input_data['Z']=[]
        input_data['V_R']=[]
        input_data['V_phi']=[]
        input_data['V_Z']=[]
        input_data['weight']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['R'].append(float(split_line[0]))
            input_data['phi'].append(float(split_line[1]))
            input_data['Z'].append(float(split_line[2]))
            input_data['V_R'].append(float(split_line[3]))
            input_data['V_phi'].append(float(split_line[4]))
            input_data['V_Z'].append(float(split_line[5]))
            input_data['weight'].append(float(split_line[6]))
     
        input_data['R']=np.asarray(input_data['R']) #convert to arrays
        input_data['phi']=np.asarray(input_data['phi'])
        input_data['Z']=np.asarray(input_data['Z'])
        input_data['V_R']=np.asarray(input_data['V_R'])
        input_data['V_phi']=np.asarray(input_data['V_phi'])
        input_data['V_Z']=np.asarray(input_data['V_Z'])
        input_data['weight']=np.asarray(input_data['weight'])
        input_data['E']=np.asarray((input_data['V_R']**2+input_data['V_phi']**2+input_data['V_Z']**2)*.5*constants.species_mass/constants.species_charge)
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading weighted full orbit beam deposition from LOCUST")
 
    return input_data

#def read_beam_depo_LOCUST_guiding_centre(): missing this option in LOCUST

def read_beam_depo_LOCUST_guiding_centre_weighted(filepath,**properties):
    """
    reads birth profile stored in LOCUST format - R phi Z V_parallel V weight

    notes:
        calculates energy in eV
    """

    print("reading weighted guiding centre beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+str(filepath))
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['absorption_fraction']=np.array(float(lines[0]))
        input_data['absorption_scaling']=np.array(float(lines[1]))
        del(lines[0])
        del(lines[0])
     
        input_data['R']=[] #initialise the arrays 
        input_data['phi']=[]
        input_data['Z']=[]
        input_data['V_pitch']=[]
        input_data['E']=[]
        input_data['weight']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['R'].append(float(split_line[0]))
            input_data['phi'].append(float(split_line[1]))
            input_data['Z'].append(float(split_line[2]))
            input_data['V_pitch'].append(float(split_line[3])/float(split_line[4]))
            input_data['E'].append(float(split_line[4])**2)
            input_data['weight'].append(float(split_line[5]))
     
        input_data['R']=np.asarray(input_data['R']) #convert to arrays
        input_data['phi']=np.asarray(input_data['phi'])
        input_data['Z']=np.asarray(input_data['Z'])
        input_data['V_pitch']=np.asarray(input_data['V_pitch'])
        input_data['E']=np.asarray(input_data['E'])*.5*constants.species_mass/constants.species_charge
        input_data['weight']=np.asarray(input_data['weight'])
        input_data['number_particles']=np.asarray(len(input_data['R']))
     
    print("finished reading weighted guiding centre beam deposition from LOCUST")
 
    return input_data

def read_beam_depo_IDS(shot,run,**properties):
    """
    reads birth profile from a distribution_sources IDS and returns as a dictionary
 
    notes:
        retains coordinate names from IMAS entries - must overwrite these manually
        reads in an arbitrary number of coordinates and injectors for each source        
        assumes that all sources have the same coordinate structure
        assumes markers hold only one time slice, at ...markers[0]
    """
    print("reading beam deposition from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_beam_depo_IDS could not import IMAS module!\nreturning\n")
        return

    input_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    if 'username' not in properties: properties['username']=settings.username
    if 'imasdb' not in properties: properties['imasdb']=settings.imasdb
    if 'imas_version' not in properties: properties['imas_version']=settings.imas_version
    input_IDS.open_env(properties['username'],properties['imasdb'],properties['imas_version'])
    input_IDS.distribution_sources.get() #open the file and get all the data from it

    input_data = {} #initialise blank dictionary to hold the data
    input_data['weight']=[]
    for identifier in input_IDS.distribution_sources.source[0].markers[0].coordinate_identifier: #generate keys for input_data by looking at the coordinates of the particle markers
        input_data[identifier.name.replace('\x00','').strip()]=[] #need to remove the unicode bits

    for source in input_IDS.distribution_sources.source: #cycle through all possible sources    
        if len(source.markers)>0:
            if len(source.markers[0].positions)>0:

                for coordinate_index in range(len(source.markers[0].positions[0,:])): #loop over the possible coordinate types e.g. r, phi, z
                    coordinate_name=source.markers[0].coordinate_identifier[coordinate_index].name.replace('\x00','').strip()

                    for marker in source.markers[0].positions[:,coordinate_index]: #this range should/must be the same for all values of coordinate_index

                        input_data[coordinate_name].extend([marker])    

            if len(source.markers[0].weights)>0: #if markers have defined weights
                input_data['weight'].extend(source.markers[0].weights)

    for key in input_data: #convert to numpy arrays
        input_data[key]=np.asarray(input_data[key])

    #check for common field names to convert to LOCUST_IO variable names

    #possible field names for different codes
    nemo_names=['Energy','Rhotor','V_PHI','Pitch angle'] 
    bbnbi_names=['z','vx','vy','vz','energy','pitch']

    #matching LOCUST_IO fields that we want to have instead
    locust_io_names=[
    ['E','rho_tor','V_phi','V_pitch'],
    ['Z','V_X','V_Y','V_Z','E','V_pitch']
    ]

    for code_counter,code_names in enumerate([nemo_names,bbnbi_names]):
        for code_name,locust_io_name in zip(code_names,locust_io_names[code_counter]):
            if code_name in input_data.keys():
                input_data[locust_io_name]=input_data.pop(code_name)
     
    input_IDS.close()

    print("finished reading beam deposition from IDS")
 
    return input_data

def read_beam_depo_TRANSP_fbm(filepath,**properties):
    """
    reads birth profile from TRANSP ASCII file at particle position
 
    notes:
        assumes all distances in cm in TRANSP file
    """

    print("reading beam deposition from TRANSP FBM format")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_TRANSP_fbm() cannot read from "+str(filepath))

        for counter,line in enumerate(lines): #look for the start of the data, marked by a certain string
            if str(line.split()[0])=='<start-of-data>':
                del(lines[0:counter+1])
                break

        input_data = {} #initialise the dictionary to hold the data
        input_data['X']=[] #initialise the arrays
        input_data['Y']=[]
        input_data['Z']=[]
        input_data['R']=[]  
        input_data['phi']=[]
        input_data['V_X']=[]
        input_data['V_Y']=[]
        input_data['V_Z']=[]
        input_data['V_R']=[]
        input_data['V_phi']=[]

        for line in lines:
            split_line=line.split()
            split_line[0]=float(split_line[0])
            split_line[1]=float(split_line[1])
            split_line[2]=float(split_line[2])
            split_line[3]=float(split_line[3])
            split_line[4]=float(split_line[4])
            split_line[5]=float(split_line[5])

            input_data['X'].append(split_line[0]) #only read in x,y,z with append for speed
            input_data['Y'].append(split_line[1])
            input_data['Z'].append(split_line[2])
            input_data['V_X'].append(split_line[3])
            input_data['V_Y'].append(split_line[4])
            input_data['V_Z'].append(split_line[5])

        input_data['X']=0.01*np.asarray(input_data['X']) #convert to arrays and from cm to m
        input_data['Y']=0.01*np.asarray(input_data['Y'])
        input_data['Z']=0.01*np.asarray(input_data['Z'])
        input_data['V_X']=0.01*np.asarray(input_data['V_X'])
        input_data['V_Y']=0.01*np.asarray(input_data['V_Y'])
        input_data['V_Z']=0.01*np.asarray(input_data['V_Z'])
        
        input_data['R']=np.asarray(np.sqrt(input_data['X']**2+input_data['Y']**2)) #need to convert from x,y,z 
        input_data['phi']=np.asarray(np.arctan2(input_data['Y'],input_data['X']))
        input_data['V_R']=np.asarray(input_data['V_X']*np.cos(input_data['phi'])+input_data['V_Y']*np.sin(input_data['phi']))
        input_data['V_phi']=np.asarray(-input_data['V_X']*np.sin(input_data['phi'])+input_data['V_Y']*np.cos(input_data['phi']))
        input_data['E']=np.asarray((input_data['V_X']**2+input_data['V_Y']**2+input_data['V_Z']**2)*.5*constants.species_mass)
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading beam deposition from TRANSP FBM format")

    return input_data

def read_beam_depo_TRANSP_fbm_guiding_centre(filepath,**properties):
    """
    reads birth profile from TRANSP ASCII file at guiding centre
 
    notes:
    """

    print("reading beam deposition from TRANSP FBM guiding centre format")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_TRANSP_fbm_guiding_centre() cannot read from "+str(filepath))

        for counter,line in enumerate(lines): #look for the start of the data, marked by a certain string
            if str(line.split()[0])=='<start-of-data>':
                del(lines[0:counter+1])
                break

        input_data = {} #initialise the dictionary to hold the data
        input_data['R']=[] #initialise the arrays
        input_data['Z']=[]
        input_data['V_pitch']=[]  
        input_data['E']=[]
        input_data['phi']=[]

        for line in lines:
            split_line=line.split()
            split_line[0]=float(split_line[0])
            split_line[1]=float(split_line[1])
            split_line[2]=float(split_line[2])
            split_line[3]=float(split_line[3])
            split_line[4]=float(split_line[4])

            input_data['R'].append(split_line[0]) #only read in x,y,z with append for speed
            input_data['Z'].append(split_line[1])
            input_data['V_pitch'].append(split_line[2])
            input_data['E'].append(split_line[3])
            input_data['phi'].append(split_line[4])

        input_data['R']=0.01*np.asarray(input_data['R']) #convert to arrays and from cm to m
        input_data['Z']=0.01*np.asarray(input_data['Z'])
        input_data['V_pitch']=np.asarray(input_data['V_pitch'])
        input_data['E']=np.asarray(input_data['E'])
        input_data['phi']=2.0*constants.pi*np.asarray(input_data['phi'])/360.0 #convert to radians
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading beam deposition from TRANSP FBM guiding centre format")

    return input_data

def read_beam_depo_TRANSP_birth(filepath,**properties):
    """
    reads birth profile from TRANSP birth CDF file at particle position

    notes:
        currently missing velocity vector - waiting for Princeton
    """


    print("reading beam deposition from TRANSP birth CDF format")

    try:
        from scipy.io import netcdf as ncdf
    except:
        raise ImportError("ERROR: scipy.io.netcdf could not be imported!\nreturning\n")
        return

    file=ncdf.netcdf_file(filepath,'r')
    input_data={}

    input_data['R']=file.variables['bs_r_D_MCBEAM'].data*.01 #convert to metres
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*constants.pi/360 #convert to radians
    input_data['Z']=file.variables['bs_z_D_MCBEAM'].data*.01 
    input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
    input_data['E']=file.variables['bs_einj_D_MCBEAM'].data
    input_data['V_pitch']=-1.*file.variables['bs_xksid_D_MCBEAM'].data
    input_data['weight']=file.variables['bs_wght_D_MCBEAM'].data
    input_data['number_particles']=np.asarray(len(input_data['R']))

    file.close()
    file.close()

    print("finished reading beam deposition from TRANSP birth CDF format")
    
    return input_data

def read_beam_depo_TRANSP_birth_guiding_centre(filepath,**properties):
    """
    reads birth profile from TRANSP birth CDF file at guiding centre

    notes:
    """

    print("reading beam deposition from TRANSP birth CDF guiding centre format")
    
    try:
        from scipy.io import netcdf as ncdf
    except:
        raise ImportError("ERROR: scipy.io.netcdf could not be imported!\nreturning\n")
        return

    file=ncdf.netcdf_file(filepath,'r')

    input_data={}

    input_data['R']=file.variables['bs_rgc_D_MCBEAM'].data*.01
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*constants.pi/360 #convert to radians
    input_data['Z']=file.variables['bs_zgc_D_MCBEAM'].data*.01
    input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
    input_data['E']=file.variables['bs_einj_D_MCBEAM'].data
    input_data['V_pitch']=-1.*file.variables['bs_xksid_D_MCBEAM'].data
    input_data['weight']=file.variables['bs_wght_D_MCBEAM'].data
    input_data['number_particles']=np.asarray(len(input_data['R']))

    file.close()
    file.close()

    print("finished reading beam deposition from TRANSP birth CDF guiding centre format")

    return input_data

def read_beam_depo_ASCOT_full_orbit(filepath,**properties):
    """
    reads birth profile from full orbit ASCOT birth ASCII file

    notes:
        typically named input.particles
    """

    with open(filepath,'r') as file: #open file

        print("reading beam deposition from ASCOT full orbit input.particles format")

        for line in file:
            if 'Number of particles' in line:
                number_particles=int(line.split()[0])
            if 'Number of different fields' in line:
                number_fields=int(line.split()[0])
                break

        fields=[] #this will hold the names of the quantities stored in the file - in order
        counter=0
        for line in file:
            fields.append(line.split()[0])
            counter+=1
            if counter==number_fields:
                break

        blank_line=file.readline() #remove read blank line XXX check this

        raw_data={} #create data dictionaries for holding all data and final returned data
        input_data={}
        for field in fields:
            raw_data[field]=[]
        
        for line_number in range(number_particles):
            line=file.readline() #XXX check this
            for number,field in zip(line.split(),fields):
                raw_data[field].extend([float(number)])

        #rename variables to LOCUST_IO conventions
        ascot_names=['energy','rho','phiprt','Rprt','zprt','vphi','vR','vz','weight'] #possible ASCOT fields
        locust_io_names=['E','rho','phi','R','Z','V_phi','V_R','V_Z','weight'] #corresponding LOCUST_IO fields that we want to retain
        for ascot_name,locust_io_name in zip(ascot_names,locust_io_names):
            if ascot_name in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])

        input_data['phi']*=2.*constants.pi/360.
        input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading beam deposition from ASCOT full orbit input.particles format")

    return input_data

def read_beam_depo_ASCOT_guiding_centre(filepath,**properties):
    """
    reads birth profile from guiding centre ASCOT birth ASCII file 

    notes:
        typically named input.particles
    """

    with open(filepath,'r') as file: #open file

        print("reading beam deposition from ASCOT guiding centre input.particles format")

        for line in file:
            if 'Number of particles' in line:
                number_particles=int(line.split()[0])
            if 'Number of different fields' in line:
                number_fields=int(line.split()[0])
                break

        fields=[] #this will hold the names of the quantities stored in the file - in order
        counter=0
        for line in file:
            fields.append(line.split()[0])
            counter+=1
            if counter==number_fields:
                break

        blank_line=file.readline() #remove read blank line XXX check this

        raw_data={} #create data dictionaries for holding all data and final returned data
        input_data={}
        for field in fields:
            raw_data[field]=[]
        
        for line_number in range(number_particles):
            line=file.readline() #XXX check this
            for number,field in zip(line.split(),fields):
                raw_data[field].extend([float(number)])

        #rename variables to LOCUST_IO conventions
        ascot_names=['energy','pitch' ,'rho','phi','R','z','vphi','vR','vz','weight'] #possible ASCOT fields
        locust_io_names=['E','V_pitch','rho','phi','R','Z','V_phi','V_R','V_Z','weight'] #corresponding LOCUST_IO fields that we want to retain
        for ascot_name,locust_io_name in zip(ascot_names,locust_io_names):
            if ascot_name in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])
    
        input_data['phi']*=2.*constants.pi/360.
        input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading beam deposition from ASCOT guiding centre input.particles format")

    return input_data

def read_beam_depo_ASCOT_full_orbit_hdf5_ini(filepath,**properties):
    """
    reads full orbit birth profile from ASCOT inistate stored in output hdf5 file 

    notes:
    """

    print("reading beam deposition from ASCOT full orbit hdf5 inistate")

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_beam_depo_ASCOT_full_orbit_hdf5_ini could not import h5py module!\n") 
        return

    with h5py.File(filepath,'r') as file:
        
        ascot_names=['energy','rho','phiprt','Rprt','zprt','vphi','vR','vz','weight'] #possible ASCOT fields
        locust_io_names=['E','rho','phi','R','Z','V_phi','V_R','V_Z','weight'] #corresponding LOCUST_IO fields that we want to retain
        input_data={}

        for ascot_name,locust_io_name in zip(ascot_names,locust_io_names):    
            input_data[locust_io_name]=file['inistate/'+ascot_name].value

        input_data['phi']*=2.*constants.pi/360.

    print("finished reading beam deposition from ASCOT full orbit hdf5 inistate")

    return input_data

def read_beam_depo_ASCOT_guiding_centre_hdf5_ini(filepath,**properties):
    """
    reads guiding centre birth profile from ASCOT inistate stored in output hdf5 file

    notes:
    """

    print("reading beam deposition from ASCOT guiding centre hdf5 inistate")

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_beam_depo_ASCOT_guiding_centre_hdf5_ini could not import h5py module!\n") 
        return

    with h5py.File(filepath,'r') as file:

        ascot_names=['energy','pitch' ,'rho','phi','R','z','vphi','vR','vz','weight'] #possible ASCOT fields
        locust_io_names=['E','V_pitch','rho','phi','R','Z','V_phi','V_R','V_Z','weight'] #corresponding LOCUST_IO fields that we want to retain
        input_data={}

        for ascot_name,locust_io_name in zip(ascot_names,locust_io_names):    
            input_data[locust_io_name]=file['inistate/'+ascot_name].value

        input_data['phi']*=2.*constants.pi/360.

    print("finished reading beam deposition from ASCOT guiding centre hdf5 inistate")
    
    return input_data

def read_beam_depo_SPIRAL_FO(filepath,**properties):
    """
    reads birth profile from full orbit SPIRAL birth ASCII file 

    notes:
    """

    with open(filepath,'r') as file:

        print("reading beam deposition from SPIRAL full orbit format")

        input_data={}
        column_numbers={} #keep track of which quantity is stored in which column

        split_line=file.readline().split()
        input_data['A']=np.asarray(split_line[1])
        input_data['Z']=np.asarray(split_line[3])
        split_line=file.readline() #read number of particles line
        split_line=file.readline().split() #read dat field headers

        for counter,quantity in enumerate(split_line): #add field headers as new quantities in input_data
            input_data[quantity]=[]
            column_numbers[str(counter)]=quantity

        split_line=file.readline().strip()#assume this is the <start-of-data> line
        if split_line!='<start-of-data>':
            print("ERROR: read_beam_depo_SPIRAL_FO encountered unknown file format!\nreturning!\n")
            return

        lines=file.readlines()
        for line in lines:
            split_line=line.split()
            for counter,number in enumerate(split_line):
                input_data[column_numbers[str(counter)]].extend([float(number)])

        #rename and recalculate variables to LOCUST_IO conventions
        spiral_names=['x[cm]','y[cm]','z[cm]','v_x[cm/s]','v_y[cm/s]','v_z[cm/s]'] #possible SPIRAL fields
        locust_io_names=['X','Y','Z','V_X','V_Y','V_Z'] #corresponding LOCUST_IO fields that we want to retain
        
        for spiral_name,locust_io_name in zip(spiral_names,locust_io_names):
            if spiral_name in input_data.keys():
                input_data[spiral_name]=0.01*np.asarray(input_data[spiral_name]) #convert everything from cm to m and to np arrays
                input_data[locust_io_name]=input_data.pop(spiral_name)

        input_data['R']=np.asarray(np.sqrt(input_data['X']**2+input_data['Y']**2)) #need to convert from x,y,z 
        input_data['phi']=np.asarray(np.arctan2(input_data['Y'],input_data['X']))
        input_data['V_R']=np.asarray(input_data['V_X']*np.cos(input_data['phi'])+input_data['V_Y']*np.sin(input_data['phi']))
        input_data['V_phi']=np.asarray(-input_data['V_X']*np.sin(input_data['phi'])+input_data['V_Y']*np.cos(input_data['phi']))
        input_data['E']=np.asarray((input_data['V_X']**2+input_data['V_Y']**2+input_data['V_Z']**2)*.5*constants.species_mass)/constants.species_charge
        input_data['number_particles']=np.asarray(len(input_data['R']))

    print("finished reading beam deposition from SPIRAL full orbit format")

    return input_data

################################################################## Beam_Deposition write functions 

def dump_beam_depo_LOCUST_full_orbit(output_data,filepath,**properties):
    """
    writes birth profile to LOCUST format - R phi Z V_R V_phi V_Z 
    
    args:
        output_data - dict holding data to dump
        filepath - full path of target file
        shuffle - if True then shuffle particle list 
    notes:
        if absorption_fraction and absorption_scaling missing in output_data, assumes 1.0
        if absorption_fraction or absorption_scaling have length>1 it writes first value
    """
 
    print("writing full orbit beam deposition to LOCUST")

    if 'shuffle' not in properties: properties['shuffle']=True
    if properties['shuffle']: 
        R,phi,Z,V_R,V_phi,V_Z=processing.utils.knuth_shuffle(output_data['R'],output_data['phi'],output_data['Z'],output_data['V_R'],output_data['V_phi'],output_data['V_Z'])
    else:
        R,phi,Z,V_R,V_phi,V_Z=output_data['R'],output_data['phi'],output_data['Z'],output_data['V_R'],output_data['V_phi'],output_data['V_Z']

    with open(filepath,'w') as file: #open file
 
        for quantity in ['absorption_fraction','absorption_scaling']:
            if quantity in output_data:
                file.write("{}\n".format(processing.utils.fortran_string(output_data[quantity].item(0),13,decimals=1))) #use .item() since disregards shape of array - could be length 0 or greater
            else:
                file.write("{}\n".format(processing.utils.fortran_string(1.0,13)))

        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays
            if phi[this_particle]<0: phi[this_particle]+=2.*np.pi 
            file.write("{r}{phi}{z}{v_r}{V_phi}{v_z}\n".format(r=processing.utils.fortran_string(R[this_particle],14,6),phi=processing.utils.fortran_string(phi[this_particle],14,6),z=processing.utils.fortran_string(Z[this_particle],14,6),v_r=processing.utils.fortran_string(V_R[this_particle],14,6),V_phi=processing.utils.fortran_string(V_phi[this_particle],14,6),v_z=processing.utils.fortran_string(V_Z[this_particle],14,6)))
    
    print("finished writing full orbit beam deposition to LOCUST") 

def dump_beam_depo_LOCUST_full_orbit_weighted(output_data,filepath,**properties):
    """
    writes birth profile to LOCUST format - R phi Z V_R V_phi V_Z weight

    args:
        output_data - dict holding data to dump
        filepath - full path of target file
        shuffle - if True then shuffle particle list 
    notes:
        -DWLIST -DWREAL are corresponding LOCUST flags
        if absorption_fraction and absorption_scaling missing in output_data, assumes 1.0
    """
 
    print("writing weighted full orbit beam deposition to LOCUST")

    if 'shuffle' not in properties: properties['shuffle']=True
    if properties['shuffle']: 
        R,phi,Z,V_R,V_phi,V_Z,weight=processing.utils.knuth_shuffle(output_data['R'],output_data['phi'],output_data['Z'],output_data['V_R'],output_data['V_phi'],output_data['V_Z'],output_data['weight'])
    else:
        R,phi,Z,V_R,V_phi,V_Z,weight=output_data['R'],output_data['phi'],output_data['Z'],output_data['V_R'],output_data['V_phi'],output_data['V_Z'],output_data['weight']

    with open(filepath,'w') as file: #open file
 
        for quantity in ['absorption_fraction','absorption_scaling']:
            if quantity in output_data:
                file.write("{}\n".format(processing.utils.fortran_string(output_data[quantity].item(0),13,decimals=1)))
            else:
                file.write("{}\n".format(processing.utils.fortran_string(1.0,13)))
 
        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays
            if phi[this_particle]<0: phi[this_particle]+=2.*np.pi 
            file.write("{r}{phi}{z}{v_r}{V_phi}{v_z}{weight}\n".format(r=processing.utils.fortran_string(R[this_particle],14,6),phi=processing.utils.fortran_string(phi[this_particle],14,6),z=processing.utils.fortran_string(Z[this_particle],14,6),v_r=processing.utils.fortran_string(V_R[this_particle],14,6),V_phi=processing.utils.fortran_string(V_phi[this_particle],14,6),v_z=processing.utils.fortran_string(V_Z[this_particle],14,6),weight=processing.utils.fortran_string(weight[this_particle],14,6)))
    
    print("finished writing weighted full orbit beam deposition to LOCUST") 

#def dump_beam_depo_LOCUST_guiding_centre(): missing this option in LOCUST

def dump_beam_depo_LOCUST_guiding_centre_weighted(output_data,filepath,equilibrium,**properties):
    """
    writes weighted birth profile to LOCUST format - R phi Z V_parallel V weight
     
    args:
        output_data - dict holding data to dump
        filepath - full path of target file
        shuffle - if True then shuffle particle list 
    notes:
        assumes R,Z,V_parallel are at the guiding centre
        if absorption_fraction and absorption_scaling missing in output_data, assumes 1.0
    """
 
    print("writing weighted guiding centre beam deposition to LOCUST")

    if 'V_pitch' not in output_data:
        print("dump_beam_depo_LOCUST_weighted found no V_pitch in output_data - calculating!")
        output_data['V_pitch']=processing.utils.pitch_calc(particle_list=output_data,equilibria=[equilibrium])

    with open(filepath,'w') as file: #open file
 
        for quantity in ['absorption_fraction','absorption_scaling']:
            if quantity in output_data:
                file.write("{}\n".format(processing.utils.fortran_string(output_data[quantity].item(0),13,decimals=1)))
            else:
                file.write("{}\n".format(processing.utils.fortran_string(1.0,13)))

        V=np.sqrt(constants.species_charge*output_data['E']*2./constants.species_mass)
        V_parallel=V*output_data['V_pitch']
 
    if 'shuffle' not in properties: properties['shuffle']=True
    if properties['shuffle']: 
        R,phi,Z,V,V_parallel,weight=processing.utils.knuth_shuffle(output_data['R'],output_data['phi'],output_data['Z'],V,V_parallel,output_data['weight'])
    else:
        R,phi,Z,V,V_parallel,weight=output_data['R'],output_data['phi'],output_data['Z'],V,V_parallel,output_data['weight']

        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays
            if phi[this_particle]<0: phi[this_particle]+=2.*np.pi
            file.write("{r}{phi}{z}{V_parallel}{V}{weight}\n".format(r=processing.utils.fortran_string(R[this_particle],14,6),phi=processing.utils.fortran_string(phi[this_particle],14,6),z=processing.utils.fortran_string(Z[this_particle],14,6),V_parallel=processing.utils.fortran_string(V_parallel[this_particle],14,6),V=processing.utils.fortran_string(V[this_particle],14,6),weight=processing.utils.fortran_string(weight[this_particle],14,6)))
    
    print("finished writing weighted guiding centre beam deposition to LOCUST") 
    
def dump_beam_depo_IDS(ID,output_data,shot,run,**properties):
    """
    writes birth profile to a distribution_sources IDS
 
    notes:
    """
    
    print("writing beam deposition to IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_beam_depo_IDS could not import IMAS module!\nreturning\n")
        return

    output_IDS=imas.ids(int(shot),int(run)) 
    if 'username' not in properties: properties['username']=settings.username
    if 'imasdb' not in properties: properties['imasdb']=settings.imasdb
    if 'imas_version' not in properties: properties['imas_version']=settings.imas_version
    output_IDS.open_env(settings.username,settings.imasdb,'3') #open the IDS
    output_IDS.distribution_sources.get()
 
    #write out code properties
    output_IDS.distribution_sources.ids_properties.comment=ID #write out identification
    output_IDS.distribution_sources.code.name="LOCUST_IO"
    if settings.commit_hash_default_LOCUST_IO: output_IDS.distribution_sources.code.commit=str(settings.commit_hash_default_LOCUST_IO)
    output_IDS.distribution_sources.code.version=support.LOCUST_IO_version
    output_IDS.distribution_sources.ids_properties.homogeneous_time=1   #must set homogeneous_time variable
    output_IDS.distribution_sources.time=np.array([0.0])
     
    #add a type of source and add a time_slice for this source
    output_IDS.distribution_sources.source.resize(1) #adds a type of source here
    output_IDS.distribution_sources.source[0].markers.resize(1) #adds a time_slice here    
    output_IDS.distribution_sources.source[0].markers[0].time=0.0 #set the time of this time_slice
 
    #add definition of our coordinate basis - r,z,phi,v_r,v_z,V_phi in this case
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier.resize(9)
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].name="R" #name of coordinate
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].index=0 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].description="major radius coordinate [m]]" #description of coordinate

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].name="phi" 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].index=1 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].description="toroidal angle coordinate [rad]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="Z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=2 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="vertical coordinate [m]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[3].name="E"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[3].index=3 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[3].description="Energy [eV]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[4].name="pitch"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[4].index=4 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[4].description="V_par/V"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[5].name="rhotor"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[5].index=5 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[5].description="toroidal rho"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[6].name="V_R"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[6].index=6
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[6].description="radial velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[7].name="V_Z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[7].index=7
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[7].description="vertical velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[8].name="V_phi"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[8].index=8
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[8].description="toroidal velocity [m/s]"

    #start storing particle data
    output_IDS.distribution_sources.source[0].markers[0].weights=output_data['weight'] 
    rho_tor=np.ones(output_data['R'].size)*-1.
    output_IDS.distribution_sources.source[0].markers[0].positions=np.array([output_data['R'],
        output_data['phi'],output_data['Z'],output_data['E'],output_data['V_pitch'],rho_tor,
        output_data['V_R'],output_data['V_Z'],output_data['V_phi']]).T
 
    #'put' all the output_data into the file and close
    output_IDS.distribution_sources.put()
    output_IDS.close()

    print("finished writing beam deposition to IDS")

def dump_beam_depo_ASCOT_full_orbit(output_data,filepath,**properties):
    """
    dumps birth profile to ASCOT ACII format at particle position

    notes:
    """
    print("writing beam deposition to ASCOT format")

    with open(filepath,'w') as file:

        #write headers
        file.write(" PARTICLE INITIAL DATA FOR ASCOT\n")
        file.write(" 4 VERSION =====================\n")
        file.write("\n")
        
        file.write(" 2  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)\n")
        file.write("This is an example file for ASCOT4 particle input.\n")
        file.write("$HeadURL: https://solps-mdsplus.aug.ipp.mpg.de/repos/ASCOT/trunk/ascot4/input.particles $")
        file.write("$LastChangedRevision: 9738 $\n")
        file.write("\n")

        file.write(" {number_particles} # Number of particles\n".format(number_particles=output_data['R'].size))
        file.write("\n")

        file.write(" 14 # Number of different fields for each particle [10 first letters are significant]\n")
        file.write("Anum      - mass number of particle        (integer)\n")
        file.write("mass      - mass of the particle           (amu)\n")
        file.write("Znum      - charge number of particle      (integer)\n")
        file.write("charge    - charge of the particle         (elemental charge)\n")
        #file.write("pitch    - pitch angle cosine of particle (vpar/vtot)\n")
        file.write("phiprt    - toroidal angle of particle     (deg)\n")
        file.write("Rprt      - R of particle                  (m)\n")
        file.write("zprt      - z of particle                  (m)\n")
        file.write("vphi      - toroidal velocity of particle  (m/s)\n")
        file.write("vR        - radial velocity of particle    (m/s)\n")
        file.write("vz        - vertical velocity of particle  (m/s)\n")
        file.write("origin    - origin of the particle         ()\n")
        file.write("weight    - weight factor of particle      (particle/second)\n")
        file.write("id        - unique identifier of particle  (integer)\n")
        file.write("Tmax      - maximum time to follow the prt (s)\n")
        file.write("\n")

        if 'weight' in output_data:
            weight=output_data['weight']
        else: #if weight not supplied in data, hard code a pseudo-weight based on 1W deposited power
            beam_power=1.0
            energies_sum=np.sum(output_data['E'])
            weight=np.zeros(len(output_data['E']))+beam_power/energies_sum

        i=0 #counter for particle identifier
        for phi,R,Z,V_phi,V_R,V_Z,W in zip(output_data['phi'],output_data['R'],output_data['Z'],output_data['V_phi'],output_data['V_R'],output_data['V_Z'],weight): 
            
            line=''
            line+=processing.utils.fortran_string(2,6,0,False) #mass and charge
            line+=processing.utils.fortran_string(constants.species_mass_amu,14,5)
            line+=processing.utils.fortran_string(1,6,0,False)
            line+=processing.utils.fortran_string(1.0,14,5)
            
            #line+=processing.utils.fortran_string(pitch,13,5,False)
            
            line+=processing.utils.fortran_string(360.0*(phi/(2.*constants.pi)),18,9) #position
            line+=processing.utils.fortran_string(R,18,9)
            line+=processing.utils.fortran_string(Z,18,9)

            line+=processing.utils.fortran_string(V_phi,18,9) #velocity
            line+=processing.utils.fortran_string(V_R,18,9)
            line+=processing.utils.fortran_string(V_Z,18,9)

            line+=processing.utils.fortran_string(1.0,9,0,False) #origin
            line+=processing.utils.fortran_string(W,18,9) #weight
            line+=processing.utils.fortran_string(i,10,0,False) #ID
            line+=processing.utils.fortran_string(999.0,18,9) #Tmax

            line+="\n"
            
            file.write(line)
            i+=1

        file.write("#EOF\n")

    print("finished writing beam deposition to ASCOT format")

def dump_beam_depo_ASCOT_guiding_centre(output_data,filepath,equilibrium,**properties):
    """
    dumps birth profile to ASCOT ACII format at guiding centre

    notes:
    """

    print("writing beam deposition to ASCOT guiding centre format")

    with open(filepath,'w') as file:

        #write headers
        file.write(" PARTICLE INITIAL DATA FOR ASCOT\n")
        file.write(" 4 VERSION =====================\n")
        file.write("\n")
        
        file.write(" 2  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)\n")
        file.write("This is an example file for ASCOT4 particle input.\n")
        file.write("$HeadURL: https://solps-mdsplus.aug.ipp.mpg.de/repos/ASCOT/trunk/ascot4/input.particles $")
        file.write("$LastChangedRevision: 9738 $\n")
        file.write("\n")

        file.write(" {number_particles} # Number of particles\n".format(number_particles=output_data['R'].size))
        file.write("\n")

        file.write(" 16 # Number of different fields for each particle [10 first letters are significant]\n")
        file.write("Anum      - mass number of particle        (integer)\n")
        file.write("mass      - mass of the particle           (amu)\n")
        file.write("Znum      - charge number of particle      (integer)\n")
        file.write("charge    - charge of the particle         (elemental charge)\n")
        file.write("energy    - kinetic energy of particle     (eV)\n")
        file.write("pitch     - pitch angle cosine of particle (vpar/vtot)\n")
        file.write("phi       - toroidal angle of particle     (deg)\n")
        file.write("R         - R of particle                  (m)\n")
        file.write("z         - z of particle                  (m)\n")
        #file.write("vphi      - toroidal velocity of particle  (m/s)\n")
        #file.write("vR        - radial velocity of particle    (m/s)\n")
        #file.write("vz        - vertical velocity of particle  (m/s)\n")
        file.write("origin    - origin of the particle         ()\n")
        file.write("weight    - weight factor of particle      (particle/second)\n")
        file.write("id        - unique identifier of particle  (integer)\n")
        file.write("Tmax      - maximum time to follow the prt (s)\n")
        file.write("Bphi      - toroidal magnetic field        (T)\n")
        file.write("BR        - radial magnetic field          (T)\n")
        file.write("Bz        - vertical magnetic field        (T)\n")
        file.write("\n")

        #calculate particle energies if missing
        if 'E' not in output_data.keys():
            output_data['E']=0.5*constants.species_mass*(output_data['V_R']**2+output_data['V_phi']**2+output_data['V_Z']**2)/constants.species_charge

        #interpolate B field to particle locations with supplied equilibrium
        if not np.all([component in equilibrium.data.keys() for component in ['B_field_R','B_field_tor','B_field_Z']]):
            print("dump_beam_depo_ASCOT_gc found no B_field in equilibrium - calculating!")
            equilibrium.B_calc()

        print("dump_beam_depo_ASCOT_gc generating B_field interpolators")
        B_field_R_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field_R']) #construct interpolators here
        B_field_tor_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field_tor'])
        B_field_Z_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field_Z'])
        print("dump_beam_depo_ASCOT_gc finished generating B_field interpolators")

        if 'V_pitch' not in output_data:
            print("dump_beam_depo_ASCOT_gc found no V_pitch in output_data - calculating!")
            output_data['V_pitch']=processing.utils.pitch_calc(particle_list=output_data,equilibria=[equilibrium])

        if 'weight' in output_data:
            weight=output_data['weight']
        else: #if weight not supplied in data, hard code a pseudo-weight based on 1W deposited power
            beam_power=1.0
            energies_sum=np.sum(output_data['E'])
            weight=np.zeros(len(output_data['E']))+beam_power/energies_sum

        i=0 #counter for particle identifier
        for E,V_pitch,phi,R,Z,W in zip(output_data['E'],output_data['V_pitch'],output_data['phi'],output_data['R'],output_data['Z'],weight): 
            
            line=''
            line+=processing.utils.fortran_string(2,6,0,False) #mass and charge
            line+=processing.utils.fortran_string(constants.species_mass_amu,14,5)
            line+=processing.utils.fortran_string(1,6,0,False)
            line+=processing.utils.fortran_string(1.0,14,5)
            
            line+=processing.utils.fortran_string(E,18,9) #energy
            line+=processing.utils.fortran_string(V_pitch,13,5,False)
            
            line+=processing.utils.fortran_string(360.0*(phi/(2.*constants.pi)),18,9) #position
            line+=processing.utils.fortran_string(R,18,9)
            line+=processing.utils.fortran_string(Z,18,9)

            #line+=processing.utils.fortran_string(V_phi,18,9) #velocity
            #line+=processing.utils.fortran_string(V_R,18,9)
            #line+=processing.utils.fortran_string(V_Z,18,9)

            line+=processing.utils.fortran_string(1.0,9,0,False) #origin
            line+=processing.utils.fortran_string(W,18,9) #weight
            line+=processing.utils.fortran_string(i,10,0,False) #ID
            line+=processing.utils.fortran_string(999.0,18,9) #Tmax
            
            line+=processing.utils.fortran_string(float(B_field_tor_interpolator(R,Z)),18,9) #B field
            line+=processing.utils.fortran_string(float(B_field_R_interpolator(R,Z)),18,9) #must interpolate ad hoc to save memory
            line+=processing.utils.fortran_string(float(B_field_Z_interpolator(R,Z)),18,9)
            line+="\n"
            
            file.write(line)
            i+=1

        file.write("#EOF\n")

    print("finished writing beam deposition to ASCOT guiding centre format")

################################################################## Beam_Deposition class
 
class Beam_Deposition(classes.base_input.LOCUST_input):
    """
    class describing neutral beam deposition profile input for LOCUST
 
    inherited from LOCUST_input:
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
        data is stored such that the coordinate 'R' for all particles is stored in my_beam_deposition['R']
        therefore the phase space position of particle p is:
            (my_beam_deposition['R'][p], my_beam_deposition['phi'][p], my_beam_deposition['Z'][p], my_beam_deposition['V_R'][p], my_beam_deposition['V_phi'][p], my_beam_deposition['V_Z'][p])
    """
 
    LOCUST_input_type='beam_deposition'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties): 
        """
        read beam_deposition from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST_FO': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_FO - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST_full_orbit(self.filepath,**properties)

        elif data_format=='LOCUST_FO_weighted': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_FO_weighted - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST_full_orbit_weighted(self.filepath,**properties)

        elif data_format=='LOCUST_GC_weighted': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_GC_weighted - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST_guiding_centre_weighted(self.filepath,**properties)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from distribution_sources IDS - shot and run required\n".format(self.ID),shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_beam_depo_IDS(self.shot,self.run,**properties)

        elif data_format=='TRANSP_fbm':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from TRANSP_fbm - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm(self.filepath,**properties)

        elif data_format=='TRANSP_fbm_gc':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from TRANSP_fbm_gc - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm_guiding_centre(self.filepath,**properties)

        elif data_format=='TRANSP_birth':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from TRANSP_birth - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_birth(self.filepath,**properties) 

        elif data_format=='TRANSP_birth_gc':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from TRANSP_birth_gc - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_birth_guiding_centre(self.filepath,**properties) 

        elif data_format=='ASCOT_FO':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_FO - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_ASCOT_full_orbit(self.filepath,**properties) 

        elif data_format=='ASCOT_GC':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_GC - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_ASCOT_guiding_centre(self.filepath,**properties)

        elif data_format=='ASCOT_FO_h5_ini':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_FO_h5_ini - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_ASCOT_full_orbit_hdf5_ini(self.filepath,**properties) 

        elif data_format=='ASCOT_GC_h5_ini':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_GC_h5_ini - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_ASCOT_guiding_centre_hdf5_ini(self.filepath,**properties)

        elif data_format=='SPIRAL_FO':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from SPIRAL_FO - filename required\n".format(self.ID),filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.properties={**properties}
                self.data=read_beam_depo_SPIRAL_FO(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST_FO/LOCUST_FO_weighted/LOCUST_GC_weighted/IDS/TRANSP_fbm/TRANSP_fbm_gc/TRANSP_birth/TRANSP_birth_gc/ASCOT_FO/ASCOT_GC/ASCOT_FO_h5_ini/ASCOT_GC_h5_ini/SPIRAL_FO)\n".format(self.ID))
 
        #populate any missing data if possible
        if any([quantity not in self.data for quantity in ['X','Y']]) and all([quantity in self.data for quantity in ['R','phi']]):
            self.data['X'],self.data['Y']=processing.utils.RphiZ_to_XYZ(self.data['R'],self.data['phi'],RH=True)

        if any([quantity not in self.data for quantity in ['R','phi']]) and all([quantity in self.data for quantity in ['X','Y']]):
            self.data['R'],self.data['phi']=processing.utils.XYZ_to_RphiZ(self.data['X'],self.data['Y'])

        if any([quantity not in self.data for quantity in ['V_R','V_Z','V_phi']]) and all([quantity in self.data for quantity in ['X','Y','V_X','V_Y']]):
            self.data['V_R'],self.data['V_phi']=processing.utils.V_XYZ_to_V_RphiZ(self.data['X'],self.data['Y'],self.data['V_X'],self.data['V_Y'])

        if any([quantity not in self.data for quantity in ['V_X','V_Y']]) and all([quantity in self.data for quantity in ['phi','V_R','V_phi']]):
            self.data['V_X'],self.data['V_Y']=processing.utils.V_RphiZ_to_V_XYZ(self.data['phi'],self.data['V_R'],self.data['V_phi'])

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,equilibrium=None,**properties):
        """
        write beam_deposition to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID = {})".format(self.ID)) 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
         
        elif data_format=='LOCUST_FO':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST_FO - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_beam_depo_LOCUST_full_orbit(self.data,filepath,**properties)

        elif data_format=='LOCUST_FO_weighted':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST_FO_weighted - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_beam_depo_LOCUST_full_orbit_weighted(self.data,filepath,**properties)

        elif data_format=='LOCUST_GC_weighted':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST_GC_weighted - filename and equilibrium required\n".format(self.ID),filename,equilibrium):
                filepath=support.dir_input_files / filename
                dump_beam_depo_LOCUST_guiding_centre_weighted(self.data,filepath,equilibrium,**properties) 
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to distribution_sources IDS - shot and run required\n".format(self.ID),shot,run):
                dump_beam_depo_IDS(self.ID,self.data,shot,run,**properties)

        elif data_format=='ASCOT_FO':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to ASCOT_FO - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_beam_depo_ASCOT_full_orbit(self.data,filepath,**properties)

        elif data_format=='ASCOT_GC':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to ASCOT_GC - filename and equilibrium required\n".format(self.ID),filename,equilibrium):
                filepath=support.dir_input_files / filename
                dump_beam_depo_ASCOT_guiding_centre(self.data,filepath,equilibrium,**properties)
 
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST_FO/LOCUST_FO_weighted/LOCUST_GC_weighted/IDS/ASCOT_FO/ASCOT_GC)\n".format(self.ID))

    def plot(self,grid=False,style='histogram',weight=True,number_bins=20,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,fill=True,quivers=False,label='',ax=False,fig=False):
        """
        plots beam deposition

        notes:
            grid - grid-like object containing same 'axes' to bin against e.g. distribution_function object with ['R'] and ['Z'] data
            style - choose from scatter or histogram
            weight - toggle whether to include marker weights in histograms
            number_bins - set number of bins or levels
            axes - list of strings specifying which axes should be plotted
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - sets r,z scale to real tokamak cross section
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            line_style - set 1D line style
            fill - toggle contour fill on 2D plots
            quivers - toggle whether to add velocity vectors
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,number_bins,colmap_val=run_scripts.utils.literal_eval(axes,number_bins,colmap_val)

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
            polar=True if axes==['phi','R'] else False
            ax = fig.add_subplot(111,polar=polar)
            
        quiver_velocity_axes=['V_X','V_Y'] if ax.name is 'polar' else [f'V_{axes[0]}',f'V_{axes[0]}']
        ax.set_title(self.ID)

        if ('X' in axes and 'X' not in self.data) or ('Y' in axes and 'Y' not in self.data):
            try:
                self['X'],self['Y']=processing.utils.RphiZ_to_XYZ(self['R'],self['phi'])
            except:
                pass
        if ((('X' in axes) or (('V_X' in axes) or ('V_X' in quiver_velocity_axes))) and 'V_X' not in self.data) or (('Y' in axes) or (('V_Y' in axes) or ('V_Y' in quiver_velocity_axes)) and 'V_Y' not in self.data):
            try:
                self['V_X'],self['V_Y']=processing.utils.V_RphiZ_to_V_XYZ(self['phi'],self['V_R'],self['V_phi'])
            except:
                pass

        ndim=len(axes) #infer how many dimensions user wants to plot
        if ndim==1: #plot 1D histograms
            if weight:
                try:
                    self_binned,self_binned_edges=np.histogram(self[axes[0]],bins=number_bins,weights=self['weight'])
                except:
                    print("ERROR: beam_deposition.plot could not find weight\n")
            else:
                self_binned,self_binned_edges=np.histogram(self[axes[0]],bins=number_bins)
            self_binned_centres=(self_binned_edges[:-1]+self_binned_edges[1:])*0.5
            ax.plot(self_binned_centres,self_binned,color=colmap(colmap_val),linestyle=line_style,label=label)
            ax.set_xlabel(axes[0])

        elif ndim==2: #plot 2D histograms

            if real_scale is True: #set x and y plot limits to real scales
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

            if style=='histogram':
                if grid is not False: #bin according to pre-defined grid
                    if weight:
                        self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]],self[axes[1]],bins=[grid[axes[0]],grid[axes[1]]],weights=self['weight'])
                    else:
                        self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]],self[axes[1]],bins=[grid[axes[0]],grid[axes[1]]])
                else:
                    if weight:
                        self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]],self[axes[1]],bins=number_bins,weights=self['weight'])
                    else:
                        self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]],self[axes[1]],bins=number_bins)
                #self_binned_x and self_binned_x are first edges then converted to centres
                self_binned_x=(self_binned_x[:-1]+self_binned_x[1:])*0.5
                self_binned_y=(self_binned_y[:-1]+self_binned_y[1:])*0.5
                dx,dy=self_binned_x[1]-self_binned_x[0],self_binned_y[1]-self_binned_y[0]            
                    
                self_binned_y,self_binned_x=np.meshgrid(self_binned_y-dy/2.,self_binned_x-dx/2.) #offset ticks onto bin centres
                
                if fill:
                    ax.set_facecolor(colmap(np.amin(self_binned)))
                    mesh=ax.pcolormesh(self_binned_x,self_binned_y,self_binned,cmap=colmap,vmin=np.amin(self_binned),vmax=np.amax(self_binned))
                else:
                    mesh=ax.contour(self_binned_x,self_binned_y,self_binned,levels=np.linspace(np.amin(self_binned),np.amax(self_binned),num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=np.amin(self_binned),vmax=np.amax(self_binned))
                    if settings.plot_contour_labels:
                        ax.clabel(mesh,inline=1,fontsize=10)

                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')

            elif style=='scatter':
                if quivers:
                    V_mag=np.sqrt(self['V_R']**2+self['V_phi']**2+self['V_Z']**2)
                    #ax.quiver(self[axes[0]],self[axes[1]],-1.*self[f'V_{quiver_velocity_axes[0]}']/V_mag,-1.*self[f'V_{quiver_velocity_axes[1]}']/V_mag,color=colmap(colmap_val),scale=0.000001,headwidth=0.000001,headlength=0.000001) #tail
                    ax.quiver(self[axes[0]],self[axes[1]],self[quiver_velocity_axes[0]]/V_mag,self[quiver_velocity_axes[1]]/V_mag,color=colmap(colmap_val))
                ax.scatter(self[axes[0]],self[axes[1]],color='red',marker='x',s=1,label=self.ID)

            if axes==['R','Z']:
                if real_scale is True:                    
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                if LCFS: #plot plasma boundary
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')

            elif axes==['X','Y']:
                if real_scale is True: #set x and y plot limits to real scales
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                if LCFS: #plot plasma boundary
                    plasma_max_R=np.max(LCFS['lcfs_r'])
                    plasma_min_R=np.min(LCFS['lcfs_r'])
                    ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                    ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')          
                if limiters: #add boundaries if desired
                    ax.set_xlim(-1.0*np.max(limiters['rlim']),np.max(limiters['rlim']))
                    ax.set_ylim(-1.0*np.max(limiters['rlim']),np.max(limiters['rlim']))
                    limiters_max_R=np.max(limiters['rlim'])
                    limiters_min_R=np.min(limiters['rlim'])
                    ax.plot(limiters_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')
                    ax.plot(limiters_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')           
            
            if (ax_flag is True or fig_flag is True) and style=='histogram': #return the plot object
                return mesh

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
           
        elif ndim==3: #plot 3D scatter - assume X,Y,Z

            if style!='scatter':
                print("ERROR: beam_deposition.plot() can only plot scatter style in 3D!")
                return

            if ax_flag is False and len(axes)==3:
                ax = fig.gca(projection='3d')
            
            if LCFS: #plot periodic poloidal cross-sections in 3D
                for angle in np.linspace(0.0,2.0*constants.pi,4,endpoint=False):
                    x_points=LCFS['lcfs_r']*np.cos(angle)
                    y_points=LCFS['lcfs_r']*np.sin(angle)
                    z_points=LCFS['lcfs_z']
                    ax.plot(x_points,y_points,zs=z_points,color=settings.plot_colour_LCFS,label='LCFS')

            if limiters: #plot periodic poloidal cross-sections in 3D
                for angle in np.linspace(0.0,2.0*constants.pi,4,endpoint=False):
                    x_points=limiters['rlim']*np.cos(angle)
                    y_points=limiters['rlim']*np.sin(angle)
                    z_points=limiters['zlim']
                    ax.plot(x_points,y_points,zs=z_points,color=settings.plot_colour_limiters,label='wall')

            if real_scale is True:
                ax.set_aspect('equal')
 
            mesh=ax.scatter(self[axes[0]],self[axes[1]],self[axes[2]],color=colmap(colmap_val),s=0.1,label=self.ID)
        
        if ax_flag is False and fig_flag is False:
            plt.show() 

    def combine(self,*targets):
        """
        combine multiple beam_deposition particle lists into single object

        args:
            targets - target objects to steal data from

        usage:
            my_beam_depo.combine(another_beam_depo,yet_another_beam_depo) #take everything

        notes:
            takes ALL data from target and generates new data fields in self if they do not exist already
            able to generate empty object and then populate with data from other objects this way
        """

        keys_self=self.data.keys()
        for target in targets: #best to loop this way around since inner loop can be random order
            keys_target=target.data.keys() #keys_target is dict_keys
            for key in keys_target:
                if key not in keys_self:
                    self[key]=np.asarray([])
                    keys_self=self.data.keys() #reset keys_self after new key is created
                self[key]=np.append(self[key],target[key])


#################################
 
##################################################################
 
###################################################################################################