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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import processing.process_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_input 
except:
    raise ImportError("ERROR: LOCUST_IO/classes/base_input.py could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from constants import *
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


################################################################## Beam_Deposition read functions
 
def read_beam_depo_LOCUST(filepath):
    """
    reads birth profile stored in LOCUST format - R Z phi V_R V_Z V_tor

    notes:
    """
 
 
    print("reading beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+filepath)
     
        input_data = {} #initialise the dictionary to hold the data
        input_data['absorption_fraction']=np.array(float(lines[0]))
        input_data['absorption_scaling']=np.array(float(lines[1]))
        del(lines[0])
        del(lines[0])
     
        input_data['R']=[] #initialise the arrays 
        input_data['phi']=[]
        input_data['Z']=[]
        input_data['V_R']=[]
        input_data['V_tor']=[]
        input_data['V_Z']=[]
     
        for line in lines:
     
            split_line=line.split()
            input_data['R'].append(float(split_line[0]))
            input_data['phi'].append(float(split_line[1]))
            input_data['Z'].append(float(split_line[2]))
            input_data['V_R'].append(float(split_line[3]))
            input_data['V_tor'].append(float(split_line[4]))
            input_data['V_Z'].append(float(split_line[5]))
     
        input_data['R']=np.asarray(input_data['R']) #convert to arrays
        input_data['phi']=np.asarray(input_data['phi'])
        input_data['Z']=np.asarray(input_data['Z'])
        input_data['V_R']=np.asarray(input_data['V_R'])
        input_data['V_tor']=np.asarray(input_data['V_tor'])
        input_data['V_Z']=np.asarray(input_data['V_Z'])
        input_data['E']=np.asarray(np.sqrt(input_data['V_R']**2+input_data['V_tor']**2+input_data['V_Z']**2)*.5*mass_deuterium)

    print("finished reading beam deposition from LOCUST")
 
    return input_data

def read_beam_depo_LOCUST_weighted(filepath):
    """
    reads birth profile stored in LOCUST format - R Z phi V_parallel V weight

    notes:
        calculates energy in eV
    """

    print("reading weighted beam deposition from LOCUST")
    
    with open(filepath,'r') as file:
     
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+filepath)
     
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
        input_data['E']=np.asarray(input_data['E'])*.5*mass_deuterium/e_charge
        input_data['weight']=np.asarray(input_data['weight'])
        
    print("finished reading weighted beam deposition from LOCUST")
 
    return input_data

def read_beam_depo_IDS(shot,run):
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

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.open()
    input_IDS.distribution_sources.get() #open the file and get all the data from it

    input_data = {} #initialise blank dictionary to hold the data
    input_data['weight']=[]
    for identifier in input_IDS.distribution_sources.source[0].markers[0].coordinate_identifier: #generate keys for input_data by looking at the coordinates of the particle markers
        input_data[identifier.name.replace('\x00','').strip()]=[] #need to remove the unicode bits

    for source in input_IDS.distribution_sources.source: #cycle through all possible sources
        if len(source.markers[0].positions)>0:

            for coordinate_index in range(len(source.markers[0].positions[0,:])): #loop over the possible coordinate types e.g. r, phi, z
                coordinate_name=source.markers[0].coordinate_identifier[coordinate_index].name.replace('\x00','').strip()

                for marker in source.markers[0].positions[:,coordinate_index]: #this range should/must be the same for all values of coordinate_index

                    input_data[coordinate_name].extend([marker])    

        if len(source.markers[0].weights)>0: #if markers have defined weights
            input_data['weight'].extend(source.markers[0].weights)

    for key in input_data: #convert to numpy arrays
        input_data[key]=np.asarray(input_data[key])
 
    input_IDS.close()

    print("finished reading beam deposition from IDS")
 
    return input_data

def read_beam_depo_TRANSP_fbm(filepath):
    """
    reads birth profile from TRANSP ASCII file at particle position
 
    notes:
        assumes all distances in cm in TRANSP file
    """

    print("reading beam deposition from TRANSP FBM format")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_TRANSP_fbm() cannot read from "+filepath)

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
        input_data['V_tor']=[]

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
        input_data['V_tor']=np.asarray(-input_data['V_X']*np.sin(input_data['phi'])+input_data['V_Y']*np.cos(input_data['phi']))
        input_data['E']=np.asarray(np.sqrt(input_data['V_X']**2+input_data['V_Y']**2+input_data['V_Z']**2)*.5*mass_deuterium)

    print("finished reading beam deposition from TRANSP FBM format")

    return input_data

def read_beam_depo_TRANSP_fbm_gc(filepath):
    """
    reads birth profile from TRANSP ASCII file at guiding centre
 
    notes:
    """

    print("reading beam deposition from TRANSP FBM guiding centre format")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_beam_depo_TRANSP_fbm_gc() cannot read from "+filepath)

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
        input_data['phi']=2.0*pi*np.asarray(input_data['phi'])/360.0

    print("finished reading beam deposition from TRANSP FBM guiding centre format")

    return input_data

def read_beam_depo_TRANSP_birth(filepath):
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

    input_data['R']=file.variables['bs_r_D_MCBEAM'].data*.01
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*pi/360
    input_data['Z']=file.variables['bs_z_D_MCBEAM'].data*.01
    input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
    input_data['E']=file.variables['bs_einj_D_MCBEAM'].data
    input_data['V_pitch']=-1.*file.variables['bs_xksid_D_MCBEAM'].data
    input_data['weight']=file.variables['bs_wght_D_MCBEAM'].data

    file.close()
    file.close()

    print("finished reading beam deposition from TRANSP birth CDF format")
    
    return input_data

def read_beam_depo_TRANSP_birth_gc(filepath):
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
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*pi/360
    input_data['Z']=file.variables['bs_zgc_D_MCBEAM'].data*.01
    input_data['X'],input_data['Y']=processing.utils.RphiZ_to_XYZ(input_data['R'],input_data['phi'])
    input_data['E']=file.variables['bs_einj_D_MCBEAM'].data
    input_data['V_pitch']=-1.*file.variables['bs_xksid_D_MCBEAM'].data
    input_data['weight']=file.variables['bs_wght_D_MCBEAM'].data

    file.close()
    file.close()

    print("finished reading beam deposition from TRANSP birth CDF guiding centre format")

    return input_data

################################################################## Beam_Deposition write functions 

def dump_beam_depo_LOCUST(output_data,filepath):
    """
    writes birth profile to LOCUST format - R Z phi V_R V_Z V_tor
     
    notes:

    """
 
    print("writing beam deposition to LOCUST")

    with open(filepath,'w') as file: #open file
 
        file.write("{}\n".format(processing.utils.fortran_string(1.0,13))) #re-insert absorption fraction lines
        file.write("{}\n".format(processing.utils.fortran_string(1.0,13)))
 
        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays

            file.write("{r}{phi}{z}{v_r}{v_tor}{v_z}\n".format(r=processing.utils.fortran_string(output_data['R'][this_particle],14,6),phi=processing.utils.fortran_string(output_data['phi'][this_particle],14,6),z=processing.utils.fortran_string(output_data['Z'][this_particle],14,6),v_r=processing.utils.fortran_string(output_data['V_R'][this_particle],14,6),v_tor=processing.utils.fortran_string(output_data['V_tor'][this_particle],14,6),v_z=processing.utils.fortran_string(output_data['V_Z'][this_particle],14,6)))
    
    print("finished writing beam deposition to LOCUST") 

def dump_beam_depo_LOCUST_weighted(output_data,filepath):
    """
    writes weighted birth profile to LOCUST format - R Z phi V_parallel V weight
     
    notes:
        assumes quantities are at the guiding centre
    """
 
    print("writing weighted beam deposition to LOCUST")

    if 'V_pitch' not in output_data:
        print("dump_beam_depo_LOCUST_weighted found no V_pitch in output_data - calculating!")
        output_data['V_pitch']=processing.utils.pitch_calc_2D(output_data=output_data,some_equilibrium=equilibrium)

    with open(filepath,'w') as file: #open file
 
        file.write("{}\n".format(processing.utils.fortran_string(1.0,13))) #re-insert absorption fraction and scaling factor lines
        file.write("{}\n".format(processing.utils.fortran_string(1.0,13)))

        V=np.sqrt(e_charge*output_data['E']*2./mass_deuterium)
        V_parallel=V*output_data['V_pitch']
 
        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays

            file.write("{r}{phi}{z}{V_parallel}{V}{weight}\n".format(r=processing.utils.fortran_string(output_data['R'][this_particle],14,6),phi=processing.utils.fortran_string(output_data['phi'][this_particle],14,6),z=processing.utils.fortran_string(output_data['Z'][this_particle],14,6),V_parallel=processing.utils.fortran_string(V_parallel[this_particle],14,6),V=processing.utils.fortran_string(V[this_particle],14,6),weight=processing.utils.fortran_string(output_data['weight'][this_particle],14,6)))
    
    print("finished writing weighted beam deposition to LOCUST") 
  
def dump_beam_depo_IDS(ID,output_data,shot,run):
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
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].name="R" #name of coordinate
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].index=0 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[0].description="major radius coordinate [m]]" #description of coordinate

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].name="phi" 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].index=1 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[1].description="toroidal angle coordinate [rad]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="Z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=2 
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="vertical coordinate [m]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="V_R"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=3
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="radial velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="V_tor"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=4
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="toroidal velocity [m/s]"

    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].name="V_Z"
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].index=5
    output_IDS.distribution_sources.source[0].markers[0].coordinate_identifier[2].description="vertical velocity [m/s]"

    #start storing particle data
    output_IDS.distribution_sources.source[0].markers[0].weights=np.ones(output_data['R'].size) #define the weights, i.e. number of particles per marker 
    positions=np.array([output_data['R'],output_data['phi'],output_data['Z'],output_data['V_R'],output_data['V_tor'],output_data['V_Z']]) #create 2D array of positions
    output_IDS.distribution_sources.source[0].markers[0].positions=np.transpose(positions) #swap the indices due to data dictionary convention
 
    #'put' all the output_data into the file and close
    output_IDS.distribution_sources.put()
    output_IDS.close()

    print("finished writing beam deposition to IDS")

def dump_beam_depo_ASCOT(output_data,filepath):
    """
    dumps birth profile to ASCOT ACII format at particle position

    notes:
        assumes output_data stores energy in eV, phi in rad
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

        #hard-code weight calculation for now
        beam_power=1.0
        energies_sum=np.sum(output_data['E'])
        weight=beam_power/energies_sum

        print("writing particle list to file")
        i=0 #counter for particle identifier
        for phi,R,Z,V_tor,V_R,V_Z in zip(output_data['phi'],output_data['R'],output_data['Z'],output_data['V_tor'],output_data['V_R'],output_data['V_Z']): 
            
            line=''
            line+=processing.utils.fortran_string(2,6,0,False) #mass and charge
            line+=processing.utils.fortran_string(mass_deuterium_amu,14,5)
            line+=processing.utils.fortran_string(1,6,0,False)
            line+=processing.utils.fortran_string(1.0,14,5)
            
            #line+=processing.utils.fortran_string(pitch,13,5,False)
            
            line+=processing.utils.fortran_string(360.0*(phi/(2.*pi)),18,9) #position
            line+=processing.utils.fortran_string(R,18,9)
            line+=processing.utils.fortran_string(Z,18,9)

            line+=processing.utils.fortran_string(V_tor,18,9) #velocity
            line+=processing.utils.fortran_string(V_R,18,9)
            line+=processing.utils.fortran_string(V_Z,18,9)

            line+=processing.utils.fortran_string(1.0,9,0,False) #origin
            line+=processing.utils.fortran_string(weight,18,9) #weight
            line+=processing.utils.fortran_string(i,10,0,False) #ID
            line+=processing.utils.fortran_string(999.0,18,9) #Tmax

            line+="\n"
            
            file.write(line)
            i+=1

        file.write("#EOF\n")

    print("finished writing beam deposition to ASCOT format")

def dump_beam_depo_ASCOT_gc(output_data,filepath,equilibrium):
    """
    dumps birth profile to ASCOT ACII format at guiding centre

    notes:
        assumes output_data stores energy in eV, phi in rad
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
            output_data['E']=0.5*mass_deuterium*(output_data['V_R']**2+output_data['V_tor']**2+output_data['V_Z']**2)/e_charge

        #interpolate B field to particle locations with supplied equilibrium
        if 'B_field' not in equilibrium.data.keys(): #calculate B field if missing
            print("dump_beam_depo_ASCOT_gc found no B_field in equilibrium - calculating!")
            if 'fpolrz' not in equilibrium.data.keys(): #calculate flux if missing
                print("dump_beam_depo_ASCOT_gc found no fpolrz in equilibrium - calculating!")
                equilibrium.set(fpolrz=processing.process_input.fpolrz_calc(equilibrium))
            equilibrium.set(B_field=processing.process_input.B_calc(equilibrium))

        print("dump_beam_depo_ASCOT_gc generating B_field interpolators")
        B_field_R_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,0]) #construct interpolators here
        B_field_tor_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,1])
        B_field_Z_interpolator=processing.utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,2])
        print("dump_beam_depo_ASCOT_gc finished generating B_field interpolators")

        if 'V_pitch' not in output_data:
            print("dump_beam_depo_ASCOT_gc found no V_pitch in output_data - calculating!")
            output_data['V_pitch']=processing.utils.pitch_calc_2D(output_data=output_data,some_equilibrium=equilibrium)

        if 'weight' in output_data:
            weight=output_data['weight']
        else: #if weight not supplied in data, hard code a pseudo-weight based on 1W deposited power
            beam_power=1.0
            energies_sum=np.sum(output_data['E'])
            weight=np.zeros(len(output_data['E']))+beam_power/energies_sum

        print("writing particle list to file")
        i=0 #counter for particle identifier
        for E,V_pitch,phi,R,Z,W in zip(output_data['E'],output_data['V_pitch'],output_data['phi'],output_data['R'],output_data['Z'],weight): 
            
            line=''
            line+=processing.utils.fortran_string(2,6,0,False) #mass and charge
            line+=processing.utils.fortran_string(mass_deuterium_amu,14,5)
            line+=processing.utils.fortran_string(1,6,0,False)
            line+=processing.utils.fortran_string(1.0,14,5)
            
            line+=processing.utils.fortran_string(E,18,9) #energy
            line+=processing.utils.fortran_string(V_pitch,13,5,False)
            
            line+=processing.utils.fortran_string(360.0*(phi/(2.*pi)),18,9) #position
            line+=processing.utils.fortran_string(R,18,9)
            line+=processing.utils.fortran_string(Z,18,9)

            #line+=processing.utils.fortran_string(V_tor,18,9) #velocity
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
            (my_beam_deposition['R'][p], my_beam_deposition['phi'][p], my_beam_deposition['Z'][p], my_beam_deposition['V_R'][p], my_beam_deposition['V_tor'][p], my_beam_deposition['V_Z'][p])
    """
 
    LOCUST_input_type='beam_deposition'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties): 
        """
        read beam_deposition from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST(self.filepath) #read the file

        elif data_format=='LOCUST_weighted': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST_weighted - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST_weighted(self.filepath) #read the file
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from distribution_sources IDS - shot and run required\n",shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_beam_depo_IDS(self.shot,self.run)

        elif data_format=='TRANSP_fbm':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_fbm - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm(self.filepath) #read the file

        elif data_format=='TRANSP_fbm_gc':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_fbm_gc - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm_gc(self.filepath) #read the file

        elif data_format=='TRANSP_birth':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_birth - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_birth(self.filepath) #read the file 

        elif data_format=='TRANSP_birth_gc':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_birth_gc - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_birth_gc(self.filepath) #read the file 

        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST/IDS/TRANSP_fbm/TRANSP_fbm_gc/TRANSP_birth/TRANSP_birth_gc)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,equilibrium=None):
        """
        write beam_deposition to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_LOCUST(self.data,filepath)

        elif data_format=='LOCUST_weighted':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST_weighted - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_LOCUST_weighted(self.data,filepath)
         
        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to distribution_sources IDS - shot and run required\n",shot,run):
                dump_beam_depo_IDS(self.ID,self.data,shot,run)

        elif data_format=='ASCOT':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_ASCOT(self.data,filepath)

        elif data_format=='ASCOT_gc':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT_gc - filename and equilibrium required\n",filename,equilibrium):
                filepath=support.dir_input_files+filename
                dump_beam_depo_ASCOT_gc(self.data,filepath,equilibrium)
 
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST/LOCUST_weighted/IDS/ASCOT/ASCOT_gc)\n")

    def plot(self,some_equilibrium=False,grid=False,style='histogram',weight=True,number_bins=20,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,fill=True,ax=False,fig=False):
        """
        plots beam deposition

        notes:
            some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
            grid - grid-like object containing same 'axes' to bin against e.g. distribution_function object with ['R'] and ['Z'] data
            style - choose from scatter or histogram
            weight - toggle whether to include marker weights in histograms
            number_bins - set number of bins or levels
            axes - list of strings specifying which axes should be plotted
            LCFS - toggles whether plasma boundary is included (requires equilibrium arguement)
            limiters - toggles limiters on/off in 2D plots
            real_scale - sets r,z scale to real tokamak cross section
            colmap - set the colour map (use get_cmap names)
            fill - toggle contour fill on 2D plots
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        import scipy
        import numpy as np
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

        ndim=len(axes) #infer how many dimensions user wants to plot
        if ndim==1: #plot 1D histograms
            if weight:
                self_binned,self_binned_edges=np.histogram(self[axes[0]],bins=number_bins,weights=self['weight'])
            else:
                self_binned,self_binned_edges=np.histogram(self[axes[0]],bins=number_bins)
            self_binned_centres=(self_binned_edges[:-1]+self_binned_edges[1:])*0.5
            ax.plot(self_binned_centres,self_binned)
            ax.set_xlabel(axes[0])
            ax.set_title(self.ID)

        elif ndim==2: #plot 2D histograms

            if axes==['R','Z']: #check for commonly-used axes
                if real_scale is True: #set x and y plot limits to real scales
                    if some_equilibrium:
                        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')

            elif axes==['X','Y']:          
                if real_scale is True: 
                    if some_equilibrium:
                        ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
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
                self_binned_y,self_binned_x=np.meshgrid(self_binned_y,self_binned_x)
                
                if fill:
                    ax.set_facecolor(colmap(np.amin(self_binned)))
                    mesh=ax.pcolormesh(self_binned_x,self_binned_y,self_binned,cmap=colmap,vmin=np.amin(self_binned),vmax=np.amax(self_binned))
                else:
                    mesh=ax.contour(self_binned_x,self_binned_y,self_binned,levels=np.linspace(np.amin(self_binned),np.amax(self_binned),num=number_bins),cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(self_binned),vmax=np.amax(self_binned))
                    ax.clabel(mesh,inline=1,fontsize=10)

                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')

            elif style=='scatter':
                ax.scatter(self[axes[0]],self[axes[1]],color='red',marker='x',s=1)

            if axes==['R','Z']:
                if real_scale is True: #set x and y plot limits to real scales
                    if some_equilibrium:
                        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                if LCFS is True: #plot plasma boundary
                    ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],plot_style_LCFS) 
                if limiters is True: #add boundaries if desired
                    ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)

            elif axes==['X','Y']:
                if real_scale is True: #set x and y plot limits to real scales
                    if some_equilibrium:
                        ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                if LCFS is True: #plot plasma boundary
                    plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                    plasma_min_R=np.min(some_equilibrium['lcfs_r'])
                    ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS)
                    ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS)          
                if limiters is True: #add boundaries if desired
                    ax.set_xlim(-1.0*np.max(some_equilibrium['rlim']),np.max(some_equilibrium['rlim']))
                    ax.set_ylim(-1.0*np.max(some_equilibrium['rlim']),np.max(some_equilibrium['rlim']))
                    limiters_max_R=np.max(some_equilibrium['rlim'])
                    limiters_min_R=np.min(some_equilibrium['rlim'])
                    ax.plot(limiters_max_R*np.cos(np.linspace(0,2.0*pi,100)),limiters_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_limiters)
                    ax.plot(limiters_min_R*np.cos(np.linspace(0,2.0*pi,100)),limiters_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_limiters)           
            
            if ax_flag is True or fig_flag is True: #return the plot object
                if 'mesh' in locals():
                    return mesh

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
            ax.set_title(self.ID)
           
        elif ndim==3: #plot 3D scatter - assume X,Y,Z

            if style!='scatter':
                print("ERROR: plot_beam_deposition() can only plot scatter style in 3D!")
                return

            if ax_flag is False and len(axes)==3:
                ax = fig.gca(projection='3d')
            
            if LCFS: #plot periodic poloidal cross-sections in 3D
                for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                    x_points=some_equilibrium['lcfs_r']*np.cos(angle)
                    y_points=some_equilibrium['lcfs_r']*np.sin(angle)
                    z_points=some_equilibrium['lcfs_z']
                    ax.plot(x_points,y_points,zs=z_points,color=plot_style_LCFS)

            if limiters: #plot periodic poloidal cross-sections in 3D
                for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                    x_points=some_equilibrium['rlim']*np.cos(angle)
                    y_points=some_equilibrium['rlim']*np.sin(angle)
                    z_points=some_equilibrium['zlim']
                    ax.plot(x_points,y_points,zs=z_points,color=plot_style_limiters)

            if real_scale is True:
                ax.set_aspect('equal')
                if some_equilibrium:
                    ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D'])) 

            ax.scatter(self[axes[0]],self[axes[1]],self[axes[2]],color=colmap(np.random.uniform()),s=0.1)
        
        if ax_flag is False and fig_flag is False:
            plt.show() 
 
#################################
 
##################################################################
 
###################################################################################################