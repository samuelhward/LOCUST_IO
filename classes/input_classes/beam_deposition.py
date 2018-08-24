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
    from processing import utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from processing import process_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")
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
    from scipy.io import netcdf as ncdf
except:
    raise ImportError("ERROR: scipy.io.netcdf could not be imported!\nreturning\n")

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
 
pi=np.pi
amu=1.66053904e-27
mass_deuterium=2.0141017781*amu


################################################################## Beam_Deposition functions
 
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
        input_data['absorption_fraction']=float(lines[0])
        input_data['absorption_scaling']=float(lines[1])
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
 
def dump_beam_depo_LOCUST(output_data,filepath):
    """
    writes birth profile to LOCUST format - R Z phi V_R V_Z V_tor
     
    notes:

    """
 
    print("writing beam deposition to LOCUST")

    with open(filepath,'w') as file: #open file
 
        file.write("{}\n".format(utils.fortran_string(1.0,13))) #re-insert junk lines
        file.write("{}\n".format(utils.fortran_string(1.0,13)))
 
        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays

            file.write("{r}{phi}{z}{v_r}{v_tor}{v_z}\n".format(r=utils.fortran_string(output_data['R'][this_particle],14,6),phi=utils.fortran_string(output_data['phi'][this_particle],14,6),z=utils.fortran_string(output_data['Z'][this_particle],14,6),v_r=utils.fortran_string(output_data['V_R'][this_particle],14,6),v_tor=utils.fortran_string(output_data['V_tor'][this_particle],14,6),v_z=utils.fortran_string(output_data['V_Z'][this_particle],14,6)))
    
    print("finished writing beam deposition to LOCUST") 

def dump_beam_depo_LOCUST_weighted(output_data,filepath):
    """
    writes weighted birth profile to LOCUST format - R Z phi V_parallel V weight
     
    notes:
        assumes quantities are at the guiding centre
    """
 
    print("writing weighted beam deposition to LOCUST")

    if 'V_pitch' not in output_data:
        print("ERROR: dump_beam_depo_LOCUST_weighted() requires V_pitch!\n")
        return

    with open(filepath,'w') as file: #open file
 
        file.write("{}\n".format(utils.fortran_string(1.0,13))) #re-insert absorption fraction and scaling factor lines
        file.write("{}\n".format(utils.fortran_string(1.0,13)))

        V=np.sqrt(output_data['E']*2./mass_deuterium)
        V_parallel=V*output_data['V_pitch']
 
        for this_particle in range(output_data['R'].size): #iterate through all particles i.e. length of our dictionary's arrays

            file.write("{r}{phi}{z}{V_parallel}{V}{weight}\n".format(r=utils.fortran_string(output_data['R'][this_particle],14,6),phi=utils.fortran_string(output_data['phi'][this_particle],14,6),z=utils.fortran_string(output_data['Z'][this_particle],14,6),V_parallel=V_parallel,V=V,weight=utils.fortran_string(output_data['weight'][this_particle],14,6)))
    
    print("finished writing weighted beam deposition to LOCUST") 
 
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
 
def dump_beam_depo_IDS(ID,output_data,shot,run):
    """
    writes birth profile to a distribution_sources IDS
 
    notes:
    """
    
    print("writing beam deposition to IDS")

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

    file=ncdf.netcdf_file(filepath,'r')
    input_data={}

    input_data['R']=file.variables['bs_r_D_MCBEAM'].data*.01
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*pi/360
    input_data['Z']=file.variables['bs_z_D_MCBEAM'].data*.01
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
    
    file=ncdf.netcdf_file(filepath,'r')

    input_data={}

    input_data['R']=file.variables['bs_rgc_D_MCBEAM'].data*.01
    input_data['phi']=file.variables['bs_zeta_D_MCBEAM'].data*2.*pi/360
    input_data['Z']=file.variables['bs_zgc_D_MCBEAM'].data*.01
    input_data['E']=file.variables['bs_einj_D_MCBEAM'].data
    input_data['V_pitch']=-1.*file.variables['bs_xksid_D_MCBEAM'].data
    input_data['weight']=file.variables['bs_wght_D_MCBEAM'].data

    file.close()
    file.close()

    print("finished reading beam deposition from TRANSP birth CDF guiding centre format")

    return input_data

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
        
        file.write(" 3  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)\n")
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
            line+=utils.fortran_string(2,6,0,False) #mass and charge
            line+=utils.fortran_string(mass_deuterium,14,5)
            line+=utils.fortran_string(1,6,0,False)
            line+=utils.fortran_string(1.0,14,5)
            
            #line+=utils.fortran_string(pitch,13,5,False)
            
            line+=utils.fortran_string(360.0*(phi/(2.*pi)),18,9) #position
            line+=utils.fortran_string(R,18,9)
            line+=utils.fortran_string(Z,18,9)

            line+=utils.fortran_string(V_tor,18,9) #velocity
            line+=utils.fortran_string(V_R,18,9)
            line+=utils.fortran_string(V_Z,18,9)

            line+=utils.fortran_string(1.0,9,0,False) #origin
            line+=utils.fortran_string(weight,18,9) #weight
            line+=utils.fortran_string(i,10,0,False) #ID
            line+=utils.fortran_string(999.0,18,9) #Tmax

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
        
        file.write(" 3  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)\n")
        file.write("This is an example file for ASCOT4 particle input.\n")
        file.write("$HeadURL: https://solps-mdsplus.aug.ipp.mpg.de/repos/ASCOT/trunk/ascot4/input.particles $")
        file.write("$LastChangedRevision: 9738 $\n")
        file.write("\n")

        file.write(" {number_particles} # Number of particles\n".format(number_particles=output_data['R'].size))
        file.write("\n")

        file.write(" 18 # Number of different fields for each particle [10 first letters are significant]\n")
        file.write("Anum      - mass number of particle        (integer)\n")
        file.write("mass      - mass of the particle           (amu)\n")
        file.write("Znum      - charge number of particle      (integer)\n")
        file.write("charge    - charge of the particle         (elemental charge)\n")
        file.write("energy    - kinetic energy of particle     (eV)\n")
        #file.write("pitch     - pitch angle cosine of particle (vpar/vtot)\n")
        file.write("phi       - toroidal angle of particle     (deg)\n")
        file.write("R         - R of particle                  (m)\n")
        file.write("z         - z of particle                  (m)\n")
        file.write("vphi      - toroidal velocity of particle  (m/s)\n")
        file.write("vR        - radial velocity of particle    (m/s)\n")
        file.write("vz        - vertical velocity of particle  (m/s)\n")
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
            output_data['E']=0.5*mass_deuterium*(output_data['V_R']**2+output_data['V_tor']**2+output_data['V_Z']**2)

        #interpolate B field to particle locations with supplied equilibrium
        if 'B_field' not in equilibrium.data.keys(): #calculate B field if missing
            print("B_field not found in equilibrium - calculating!")
            if 'fpolrz' not in equilibrium.data.keys(): #calculate flux if missing
                print("fpolrz not found in equilibrium - calculating!")
                equilibrium.set(fpolrz=process_input.fpolrz_calc(equilibrium))
            equilibrium.set(B_field=process_input.B_calc(equilibrium))

        print("generating B_field interpolators")
        B_field_R_interpolator=utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,0]) #construct interpolators here
        B_field_tor_interpolator=utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,1])
        B_field_Z_interpolator=utils.interpolate_2D(equilibrium['R_1D'],equilibrium['Z_1D'],equilibrium['B_field'][:,:,2])
        print("finished generating B_field interpolators")


        #hard-code weight calculation for now
        beam_power=1.0
        energies_sum=np.sum(output_data['E'])
        weight=beam_power/energies_sum


        print("writing particle list to file")
        i=0 #counter for particle identifier
        for E,phi,R,Z,V_tor,V_R,V_Z in zip(output_data['E'],output_data['phi'],output_data['R'],output_data['Z'],output_data['V_tor'],output_data['V_R'],output_data['V_Z']): 
            
            line=''
            line+=utils.fortran_string(2,6,0,False) #mass and charge
            line+=utils.fortran_string(mass_deuterium,14,5)
            line+=utils.fortran_string(1,6,0,False)
            line+=utils.fortran_string(1.0,14,5)
            
            line+=utils.fortran_string(E,18,9) #energy
            #line+=utils.fortran_string(pitch,13,5,False)
            
            line+=utils.fortran_string(360.0*(phi/(2.*pi)),18,9) #position
            line+=utils.fortran_string(R,18,9)
            line+=utils.fortran_string(Z,18,9)

            line+=utils.fortran_string(V_tor,18,9) #velocity
            line+=utils.fortran_string(V_R,18,9)
            line+=utils.fortran_string(V_Z,18,9)

            line+=utils.fortran_string(1.0,9,0,False) #origin
            line+=utils.fortran_string(weight,18,9) #weight
            line+=utils.fortran_string(i,10,0,False) #ID
            line+=utils.fortran_string(999.0,18,9) #Tmax
            
            line+=utils.fortran_string(float(B_field_tor_interpolator(R,Z)),18,9) #B field
            line+=utils.fortran_string(float(B_field_R_interpolator(R,Z)),18,9) #must interpolate ad hoc to save memory
            line+=utils.fortran_string(float(B_field_Z_interpolator(R,Z)),18,9)
            line+="\n"
            
            file.write(line)
            i+=1

        file.write("#EOF\n")

    print("finished writing beam deposition to ASCOT guiding centre format")

################################################################## Beam_Deposition class
 
class Beam_Deposition(base_input.LOCUST_input):
    """
    class describing neutral beam deposition profile input for LOCUST
 
    inheritedfrom LOCUST_input:
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
 
        if utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_LOCUST(self.filepath) #read the file
         
        elif data_format=='IDS':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from distribution_sources IDS - shot and run data required\n",shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_beam_depo_IDS(self.shot,self.run)

        elif data_format=='TRANSP_fbm':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_fbm - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm(self.filepath) #read the file

        elif data_format=='TRANSP_fbm_gc':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_fbm_gc - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_fbm_gc(self.filepath) #read the file

        elif data_format=='TRANSP_birth':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_birth - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_beam_depo_TRANSP_birth(self.filepath) #read the file 

        elif data_format=='TRANSP_birth_gc':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from TRANSP_birth_gc - filename required\n",filename):
 
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
 
        if utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_LOCUST(self.data,filepath)

        elif data_format=='LOCUST_weighted':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST_weighted - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_LOCUST_weighted(self.data,filepath)
         
        elif data_format=='IDS':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to distribution_sources IDS - shot and run required\n",shot,run):
                dump_beam_depo_IDS(self.ID,self.data,shot,run)

        elif data_format=='ASCOT':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_beam_depo_ASCOT(self.data,filepath)

        elif data_format=='ASCOT_gc':
            if not utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT_gc - filename and equilibrium required\n",filename,equilibrium):
                filepath=support.dir_input_files+filename
                dump_beam_depo_ASCOT_gc(self.data,filepath,equilibrium)
 
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST/LOCUST_weighted/IDS/ASCOT/ASCOT_gc)\n")

 
 
#################################
 
##################################################################
 
###################################################################################################