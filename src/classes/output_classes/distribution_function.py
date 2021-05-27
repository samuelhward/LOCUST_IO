#distribution_function.py

'''
Samuel Ward
01/04/2018
----
class to handle LOCUST distribution function output data
---
usage:
    see README.md for usage

notes:         
---
'''


###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from output_classes import x") but best practice to import whole output_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
    import datetime
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
    import classes.base_output 
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/base_output.py could not be imported!\nreturning\n")
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


################################################################## Distribution_Function read functions

def read_distribution_function_LOCUST(filepath,**properties):
    """
    reads distribution function stored in unformatted fortran file

    notes:
        IDFTYP==1:F
                2:H
                3:X
                4:I
            else: #

        all n* variables are really n*-1 e.g. nV is nV-1 in LOCUST

        precision of LOCUST's output controlled by -DBP LOCUST flag (double precision hard-coded to be enabled by default in makefile optimisations)
        must keep an eye on the double/single GPU formats specified in LOCUST's prec_mod.f90 to ensure Dfn is read correctly from LOCUST binary format

        final Dfn is s^3/m^6
    """

    print("reading distribution function from LOCUST")

    filename=list(filepath.parts)[-1] #infer IDFTYP from first character of file name
    DFN_HEAD=filename[0] 
    if DFN_HEAD=='F': 
        IDFTYP=1
    elif DFN_HEAD=='H':    
        IDFTYP=2
    elif DFN_HEAD=='X':
        IDFTYP=3
    elif DFN_HEAD=='I':
        IDFTYP=4
    else:
        print('ERROR: cannot infer IDFTYPE from filename - {filename}\n'.format(filename))

    try:
        from scipy.io import FortranFile 
    except:
        raise ImportError("ERROR: read_distribution_function_LOCUST could not import scipy.io.FortranFile!\nreturning\n")
        return

    file=FortranFile(filepath,'r')
    input_data={} #initialise blank dictionary
    input_data['IDFTYP']=np.array(IDFTYP)

    if IDFTYP==3:

        #32-char checksum (32 bytes in Fortran, so essentially reading 8x32 bit floats)
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum
        
        if properties['EBASE']:
            #nEF-1
            input_data['nE']=file.read_ints() 
        else:
            #nV-1
            input_data['nV']=file.read_ints() #V

        #nPP-1
        input_data['nP_phi']=file.read_ints() #PPh
        #nmu-1
        input_data['nmu']=file.read_ints() #MU
        
        if properties['EBASE']:
            #dEFh
            input_data['dE']=file.read_reals(dtype=np.float32)         
        else:
            #dVFh
            input_data['dV']=file.read_reals(dtype=np.float32) 
        
        #dPPh
        input_data['dP_phi']=file.read_reals(dtype=np.float32) #dPPh    = P_phi1h - P_phi0h
        #dMUh
        input_data['dmu']=file.read_reals(dtype=np.float32)

        if properties['EBASE']:
            #EF+dEF/2 (nV long)
            input_data['E']=file.read_reals(dtype=np.float32)
        else:
            #VF+dVF/2 (nV long)
            input_data['V']=file.read_reals(dtype=np.float32)

        #PP+dPP/2 (nPP long)
        input_data['P_phi']=file.read_reals(dtype=np.float32) #PPhi
        #MU+MU/2 (nmu long)
        input_data['mu']=file.read_reals(dtype=np.float32)
        
        #for now these blocks do the same thing, but useful to leave as separate for now
        if properties['wtot']: #cumulative energy inventory (total energy injected so far)
            if properties['EBASE']:
                #Fh_norm (nEF by nPP by nmu)
                input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
            else:
                #Fh_norm (nV by nPP by nmu)
                input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
        else:
            #Fh (NOTE check dimensions e.g. nV by nPP by nmu)
            input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
    
        if True:#Fh_s:
            #Fh_s (nV by nPP by nmu) if Fh_s present
            input_data['dfn_s']=file.read_reals(dtype=np.float32) #Dfn. M.C. error
        if properties['J']:
            #Jh
            input_data['J']=file.read_reals(dtype=np.float64) #Jacobian for IDFTYP=3 
        if properties['J_s']:
            #Jh_s
            input_data['J_s']=file.read_ints(dtype=np.float64) #Jacobian error for IDFTYP=3
            if properties['cpu_time']:
                #cpuTime
                input_data['cpu_time']=file.read_reals(dtype=np.float64) #this may not be present in the file

        input_data['ddir']=np.array([1.]) #first dimension is just sign, so jacobian is unity
        input_data['dir']=np.array([-1,1.]) #first dimension is just sign, so jacobian is unity

    else:

        #32-char checksum
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum
        #zero indicates fast ions only
        input_data['0_1']=file.read_ints()

        if IDFTYP==4:
            #nPSIF-1
            input_data['npsi']=file.read_ints() #number of surface contours
            #nPOLF-1
            input_data['npol']=file.read_ints() #number of poloidal cells
        else:
            #nF-1
            input_data['nR']=file.read_ints() #R and Z dimension of the distribution function grid (equal here)
            input_data['nZ']=file.read_ints()

        if IDFTYP==1 or IDFTYP==4:
            #nL-1 
            input_data['nV_pitch']=file.read_ints() #Vphi/V cell boundaries
        else:
            #nPP-1
            input_data['nP_phi']=file.read_ints() #PPhi          

        if properties['EBASE']:
            #nE-1
            input_data['nE']=file.read_ints() #E
        else:
            #nV-1
            input_data['nV']=file.read_ints() #V        
        
        #nP-1
        input_data['nphase']=file.read_ints() #poloidal gyro-phase cell boundaries
        
        if properties['ITER']: #NOTE these blocks essentially do the same thing right now

            if properties['WIPE']:
                input_data['dfn']=[] #Final combined DFn. grid
                input_data['dfn_s']=[] #Dfn. M.C. error
                for line in range(int(input_data['nphase'])):
                    input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                for line in range(int(input_data['nphase'])):
                    input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                
            else:
                input_data['dfn']=[] #Final combined DFn. grid  
                input_data['dfn_s']=[] #Dfn. M.C. error
                for line in range(int(input_data['nphase'])): #Dfn. M.C. error
                    input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                for line in range(int(input_data['nphase'])):
                    input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
        
        else:
            input_data['dfn']=[] #Final combined DFn. grid   
            input_data['dfn_s']=[] #Dfn. M.C. error
            for line in range(int(input_data['nphase'])): #Dfn. M.C. error
                input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
            for line in range(int(input_data['nphase'])):
                input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long

        if input_data['IDFTYP']==2:
            input_data['dfn']=np.swapaxes(input_data['dfn'],2,3) #swap final order to ...r,z - means plotting functions can assume index order x,y
        
        nEQ=file.read_ints() 

        input_data['nR_1D']=np.array(nEQ[0]) #2D field grid R dimension
        input_data['nZ_1D']=np.array(nEQ[1]) #2D field grid Z dimension
   
        #LOCUST distribution function writing diverged in style after certain date - check for this
        day=int(filename[2:4])
        month=int(filename[5:7])
        year=int(filename[8:12])
        filedate=datetime.date(day=day,month=month,year=year)
        change_date=datetime.date(day=15,month=4,year=2019)
        if filedate>change_date:
            input_data['R_1D']=file.read_reals(dtype=np.float32) #Eqm. R grid. (nR_1D long) 
            input_data['Z_1D']=file.read_reals(dtype=np.float32) #Eqm. Z grid. (nZ_1D long)
            input_data['psirz']=file.read_reals(dtype=np.float32).reshape(input_data['nR_1D'],input_data['nZ_1D'],order='F') #PSIrz
        else:
            input_data['R_1D']=file.read_reals(dtype=np.float64) #Eqm. R grid. (nR_1D long) 
            input_data['Z_1D']=file.read_reals(dtype=np.float64) #Eqm. Z grid. (nZ_1D long)
            input_data['psirz']=file.read_reals(dtype=np.float64).reshape(input_data['nR_1D'],input_data['nZ_1D'],order='F') #PSIrz
   
        if IDFTYP==4:

            #dVOL (nPSIF by nPOLF)
            input_data['dVOL']=file.read_reals(dtype=np.float32) 
            #npn
            input_data['npn']=file.read_ints() #cell granularity for IDFTYP=4          
            #csb (nPSIF-1 by nPOLF by 2*npn by 2)
            input_data['csb']=file.read_reals(dtype=np.float32)
            #npolh (nPSIF long)
            input_data['npolh']=file.read_reals(dtype=np.float64)            
            #just the number 1.5 
            input_data['1.5']=file.read_reals(dtype=np.float32)

        else:

            #R+dR/2 (nF long)
            input_data['R']=file.read_reals(dtype=np.float32) #r space of dfn [m]
            #Z+dZ/2 (nF long)
            input_data['Z']=file.read_reals(dtype=np.float32) #z space of dfn [m]

        if IDFTYP==1 or IDFTYP==4:

            #L+dL/2 (nL long)
            input_data['V_pitch']=file.read_reals(dtype=np.float32) #pitch space of dfn

        else:

            #PP+dPP/2 (nPP long)
            input_data['P_phi']=file.read_reals(dtype=np.float32)

        if properties['EBASE']:
            #E+dE/2 (nE long)
            input_data['E']=file.read_reals(dtype=np.float32) #energy space of dfn
            input_data['E']/=constants.species_charge #convert energy to [eV]
            input_data['V']=np.array(np.sqrt(2.*input_data['E']*constants.species_charge/constants.species_mass)) #[m/s]
        else:
            #V+dV/2 (nV long)
            input_data['V']=file.read_reals(dtype=np.float32) #velocity space of dfn
            input_data['E']=np.array((0.5*constants.species_mass*input_data['V']**2)/constants.species_charge) #calculate energy [eV]
        
        #PG+dPG/2 (nP long)
        input_data['P']=file.read_reals(dtype=np.float32) #special dimension - simulation specific (e.g. gyrophase)

        input_data['Ab']=file.read_reals(dtype=np.float32) #fast ion masses
        input_data['Ai_1']=file.read_reals(dtype=np.float32) #first value of Ab
        input_data['Zb']=file.read_reals(dtype=np.float32) #trace particle Z
        input_data['Zi_1']=file.read_reals(dtype=np.float32) #first value of Zb
        input_data['Vsclh']=file.read_reals(dtype=np.float32) #vgrid upper bound
        input_data['Vnrm']=file.read_reals(dtype=np.float32) #normalising velocity
        input_data['icoll']=file.read_ints() #collisions (1=on 0=off)
        input_data['iscat']=file.read_ints() #scattering (1=on 0=off)
        input_data['idiff']=file.read_ints() #energy diffusion (1=on 0=off)
        input_data['iloss']=file.read_ints() #charge exchange losses (1=on 0=off)
        input_data['iterm']=file.read_ints() #terminate if ptcl. leaves plasma
        input_data['niter']=file.read_ints() #number of iterations for isym=1 simulation
        input_data['integrator']=file.read_ints() #integrator type
        input_data['npnt']=file.read_ints() #points per gyration
        input_data['one_1']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['one_2']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['one_3']=file.read_reals(dtype=np.float32) #the number 1.0
        input_data['999']=file.read_reals(dtype=np.float32) #the number 999.0
        input_data['0_2']=file.read_reals(dtype=np.float32) #the maximum integrator step size
        input_data['dt0']=file.read_reals(dtype=np.float32) #the number 0.0
        input_data['dt0']=file.read_reals(dtype=np.float32) #the number 0.0
        input_data['threads_per_block']=file.read_ints() #gpu threads per block 
        input_data['blocks_per_grid']=file.read_ints() #gpu blockers per grid

        if properties['TEST']:

            #Pdep/E0
            input_data['Pdep/E0']=file.read_reals(dtype=np.float32) #pdep is injected power
            #tau_s
            input_data['tau_s']=file.read_reals(dtype=np.float32) #zeroth order slowing down time
            #E0
            input_data['E0']=file.read_reals(dtype=np.float32) #energy (plasma frame)
            #EC
            input_data['EC']=file.read_reals(dtype=np.float32) #zeroth order critical energy
            #rho=Ai_1/(2*Ab)
            input_data['rho']=file.read_reals(dtype=np.float32)
            #siglg
            input_data['siglg']=file.read_reals(dtype=np.float32) #r.m.s. width of test src

            if properties['cpu_time']:
                input_data['cpu_time']=file.read_reals(dtype=np.float64)

        #extra derived data
        if input_data['nphase']>1:
            input_data['dP']=np.array(input_data['P'][1]-input_data['P'][0]) #special dimension bin width 
        else:
            input_data['dP']=np.array(2.*constants.pi)
        
    #some post processing
    if input_data['IDFTYP']==1:

        if properties['EBASE']:
            input_data['dfn_index']=np.array(['P','E','V_pitch','R','Z']) #reference for names of each dfn dimension
            for key in ['dfn','dfn_s']:
                input_data[key]=np.array(input_data[key]).reshape(int(input_data['nphase']),int(input_data['nE']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F')  
            input_data['nc']=len(input_data['dfn'])/input_data['nphase'] #nV*nV_pitch*nZ*nR (nZ, nR = nF)
        else:
            input_data['dfn_index']=np.array(['P','V','V_pitch','R','Z']) #reference for names of each dfn dimension
            for key in ['dfn','dfn_s']:
                input_data[key]=np.array(input_data[key]).reshape(int(input_data['nphase']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F')  
            input_data['nc']=len(input_data['dfn'])/input_data['nphase'] #nV*nV_pitch*nZ*nR (nZ, nR = nF)
        input_data['dfn']=np.swapaxes(input_data['dfn'],3,4) #swap final order to ...r,z - means plotting functions can assume index order x,y

        input_data['dfn']*=0.5 #convert from same pitch dimension as TRANSP (per dSolidAngle/4pi) to per unit pitch 
        input_data['dR']=np.array(input_data['R'][1]-input_data['R'][0]) #R bin width
        input_data['dZ']=np.array(input_data['Z'][1]-input_data['Z'][0]) #Z bin width
        input_data['dV_pitch']=np.array(input_data['V_pitch'][1]-input_data['V_pitch'][0]) #pitch bin width
        input_data['dV']=np.array(input_data['V'][1]-input_data['V'][0]) #velocity bin width
        input_data['dE']=np.array(np.abs(input_data['E'][1]-input_data['E'][0])) #energy bin width

    elif input_data['IDFTYP']==2:
        if properties['EBASE']:
            input_data['dfn_index']=np.array(['E','P_phi','R','Z']) #reference for names of each dfn dimension
        else:
            input_data['dfn_index']=np.array(['V','P_phi','R','Z']) #reference for names of each dfn dimension

    elif input_data['IDFTYP']==3:
        if properties['EBASE']:            
            for key in ['dfn','dfn_s','J','J_s']:
                input_data[key]=np.array(input_data[key]).reshape(2,int(input_data['nE']),int(input_data['nP_phi']),int(input_data['nmu']),order='F')  
            input_data['dfn_index']=np.array(['dir','E','P_phi','mu']) #reference for names of each dfn dimension
        else:
            input_data['dfn_index']=np.array(['dir','V','P_phi','mu']) #reference for names of each dfn dimension
            for key in ['dfn','dfn_s','J','J_s']:
                input_data[key]=np.array(input_data[key]).reshape(2,int(input_data['nV']),int(input_data['nP_phi']),int(input_data['nmu']),order='F')

    elif input_data['IDFTYP']==4:
        if properties['EBASE']:
            input_data['dfn_index']=np.array(['P','E','V_pitch','pol','psi']) #reference for names of each dfn dimension
        else:
            input_data['dfn_index']=np.array(['P','V','V_pitch','pol','psi']) #reference for names of each dfn dimension
     
    file.close()

    print("finished reading distribution function from LOCUST")
    return input_data

def read_distribution_function_ASCOT(filepath,**properties):
    """
    read distribution as written out in hdf5 file by ASCOT

    notes:

    """

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_distribution_function_ASCOT could not import h5py module!\n") 
        return

    with h5py.File(filepath,'r') as file:

        print("reading distribution function from ASCOT")

        input_data={}

        ''' 
        extract dfn to #[m^-3 eV^-1 dpitch^-1] format

        file['distributions/rzPitchEdist/ordinate'].value=[species,time,E,V_pitch,Z,R,ordinate]
        file['distributions/rzPitchEdist/abscissae']['dim1']=R
                                                 ...['dim2']=Z
                                                 ...['dim3']=V_pitch
                                                 ...['dim4']=E [J]
                                                 ...['dim5']=time
                                                 ...['dim6']=species'''
            
        input_data['R']=.5*(file['distributions/rzPitchEdist/abscissae/dim1'].value[1:]+file['distributions/rzPitchEdist/abscissae/dim1'].value[:-1]) #bin centres [m]
        input_data['Z']=.5*(file['distributions/rzPitchEdist/abscissae/dim2'].value[1:]+file['distributions/rzPitchEdist/abscissae/dim2'].value[:-1]) #[m]
        input_data['V_pitch']=.5*(file['distributions/rzPitchEdist/abscissae/dim3'].value[1:]+file['distributions/rzPitchEdist/abscissae/dim3'].value[:-1])
        input_data['E']=.5*(file['distributions/rzPitchEdist/abscissae/dim4'].value[1:]+file['distributions/rzPitchEdist/abscissae/dim4'].value[:-1])/constants.species_charge #[eV]
        input_data['V']=np.sqrt(2.*input_data['E']/constants.species_mass)

        input_data['dR']=np.abs(input_data['R'][1]-input_data['R'][0])
        input_data['dZ']=np.abs(input_data['Z'][1]-input_data['Z'][0])
        input_data['dV_pitch']=np.abs(input_data['V_pitch'][1]-input_data['V_pitch'][0])
        input_data['dE']=np.abs(input_data['E'][1]-input_data['E'][0])
        
        input_data['nR']=np.array(len(input_data['R']))
        input_data['nZ']=np.array(len(input_data['Z']))
        input_data['nV_pitch']=np.array(len(input_data['V_pitch']))
        input_data['nE']=np.array(len(input_data['E']))

        input_data['dfn']=file['distributions/rzPitchEdist/ordinate'].value #[m^-3 J^-1 dpitch^-1]
        input_data['dfn']=np.sum(input_data['dfn'],axis=0)
        input_data['dfn']=np.sum(input_data['dfn'],axis=0)
        input_data['dfn']=np.sum(input_data['dfn'],axis=-1)
        input_data['dfn']=np.swapaxes(input_data['dfn'],-1,-2) 
        input_data['dfn']=input_data['dfn'][np.newaxis,:] #insert dummy P dimension
        input_data['dfn']*=constants.species_charge #[m^-3 eV^-1 dpitch^-1]
        
        input_data['dfn_index']=np.array(['P','E','V_pitch','R','Z'])
        input_data['IDFTYP']=1
    
    print("finished reading distribution function from ASCOT")

    return input_data

def read_distribution_function_particle_list(source,**properties):
    """
    generate a distribution function from a pre-existing list of markers e.g. beam dsitribution or a loss list

    notes:
    """

    print("generating reading distribution function from particle list")

    input_data={}

    #do something with an nD histogram here

    print("finished generating reading distribution function from particle list")

    return input_data

def read_distribution_function_IDS(shot,run,**properties):
    """
    read distribution function from IDS

    notes:
        distribution function is already integrated down to R Z dimensions (so use transform=False in .plot etc.)
    """

    print("reading distribution function from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_distribution_function_IDS could not import IMAS module!\nreturning\n")
        return

    input_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    if 'username' not in properties: properties['username']=settings.username
    if 'imasdb' not in properties: properties['imasdb']=settings.imasdb
    if 'imas_version' not in properties: properties['imas_version']=settings.imas_version
    input_IDS.open_env(properties['username'],properties['imasdb'],properties['imas_version'])
    input_IDS.distributions.get() #open the file and get all the data from it
 
    input_data = {} #initialise blank dictionary to hold the data

    input_data['R_1D']=np.array(input_IDS.distributions.distribution[0].profiles_2d[0].grid.r)
    input_data['Z_1D']=np.array(input_IDS.distributions.distribution[0].profiles_2d[0].grid.z)
    input_data['dfn']=np.array(input_IDS.distributions.distribution[0].profiles_2d[0].density_fast)
    input_data['dfn_index']=np.array(['R','Z'])

    input_IDS.close()

    print("finished reading distribution function from IDS")

    return input_data


################################################################## Distribution_Function write functions

def dump_distribution_function_LOCUST(output_data,filepath,**properties): 
    """
    writes distribution function to LOCUST
    
    notes:
        
    """

    pass

def dump_distribution_function_TRANSP(output_data,filepath,**properties):
    """
    notes:
        check old LOCUST_IO code for functions which can generate new netCDF variables easily
    """
    pass


################################################################## Distribution_Function class

class Distribution_Function(classes.base_output.LOCUST_output):
    """
    class describing distribution function output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'distribution_function'
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
        the properties field for the distribution_function can contain a dictionary of LOCUST flags present which dictate how the dfn is written to file

        DoF for Fh (in order of array index i.e. my_dfn['dfn'][P,V/E,V_pitch,R,Z])
            IDFTYP=1
                P       - special, rare simulation specific (e.g. gyro phase dimension) - not used in EBASE mode
                V/E     - velocity dimension/energy dimension (dictated by EBASE)
                V_pitch - pitch dimension
                R       - r dimension of bin centres
                Z       - z dimension of bin centres
    """

    LOCUST_output_type='distribution_function'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read distribution function from file 

        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename                
                for variable, default_value in zip(['ITER','wtot','WIPE','TEST','EBASE','dfn_s','J','J_s','cpu_time'],[True,False,False,False,True,True,True,True,True]): #default properties settings
                    if variable not in properties:
                        properties[variable]=default_value
                self.properties={**properties}
                self.data=read_distribution_function_LOCUST(self.filepath,**properties)

        elif data_format=='ASCOT': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from ASCOT - filename required\n".format(self.ID),filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                for variable, default_value in zip(['EBASE'],[True]): #default properties settings
                    if variable not in properties:
                        properties[variable]=default_value
                self.properties={**properties}
                self.data=read_distribution_function_ASCOT(self.filepath,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from IDS - shot and run required\n".format(self.ID),shot,run):
                
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_distribution_function_IDS(self.shot,self.run,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/ASCOT/IDS)\n".format(self.ID))            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write distribution function to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_output_files / filename
                dump_distribution_function_LOCUST(self.data,filepath,**self.properties)
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))

    def plot(self,key='dfn',axes=['R','Z'],LCFS=False,limiters=False,gridlines=False,real_scale=False,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,transform=True,number_bins=20,fill=True,vminmax=None,label='',ax=False,fig=False):
        """
        plot the distribution function

        notes:
            if external figure or axes are supplied then, if possible, function returns plottable object for use with external colorbars etc 
            if user supplies full set of indices, code assumes those slices are dimension to plot over i.e. please crop before plotting
        args:
            key - select which data to plot
            axes - define plot axes in x,y order or as full list of indices/slices (see dfn_transform())
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - plot to Tokamak scale
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            line_style - set 1D line style
            transform - set to False if supplied dfn has already been cut down to correct dimensions
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            vminmax - set mesh Vmin/Vmax values
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        usage:
            plot_distribution_function(my_dfn) #basic default R,Z plot
            plot_distribution_function(my_dfn,axes=['E','V_pitch']) #basic pitch,energy plot
            plot_distribution_function(my_dfn,my_eq,axes=['R','Z'],LCFS=True,real_scale=True) #R,Z plotted with true scales and last closed flux surface from supplied equilibrium
            plot_distribution_function(my_dfn,my_eq,axes=['R','Z'],LCFS=True,real_scale=True,transform=False) #R,Z plot where my_dfn has already been cropped to correct dimensions
            plot_distribution_function(my_dfn,axes=[0,9,3,slice(None),slice(None)],vminmax=[0,1000],ax=my_ax,fig=my_fig) #R,Z plot at point 9,3 in E,pitch space without integrating and adding to my_ax on figure my_fig
        axes options:
            R,Z - integrate over pitch, gyrophase and velocity [m]^-3
            E,V_pitch - integrate over space and transform to [eV]^-1[dpitch]^-1 
            E - [eV]^-1
            R - [m]^-3 
            N - total #
            list of indices and slices
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,number_bins,vminmax,colmap_val=run_scripts.utils.literal_eval(axes,number_bins,vminmax,colmap_val)
        
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
            ax.plot(self[axes[0]],self[key],color=colmap(colmap_val),linestyle=line_style,label=label)
            ax.set_ylabel(key)

        #2D data
        elif self[key].ndim==2:

            X=self[axes[0]] #make a mesh
            Y=self[axes[1]]
            dx,dy=X[1]-X[0],Y[1]-Y[0]

            Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres
            Z=self[key]

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.min(Z)
                vmax=np.max(Z)
            
            #2D plot
            if fill is True:
                mesh=ax.contourf(X,Y,Z,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax)
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

            if LCFS:
                ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
            if limiters: #add boundaries if desired
                ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')

            if real_scale is True: #set x and y plot limits to real scales
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh
                
        #plot distribution function
        elif key=='dfn':
            
            #transform distribution function to the coordinates we want
            if transform is True:
                dfn_copy=self.transform(axes=axes) #user-supplied axes are checked for validity here
            else:
                dfn_copy=copy.deepcopy(self)

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.amin(dfn_copy[key])
                vmax=np.amax(dfn_copy[key])

            #check resulting dimensionality of distribution function
            if dfn_copy['dfn'].ndim==0: #user has given 0D dfn
                pass #XXX incomplete - should add scatter point
            elif dfn_copy['dfn'].ndim==1: #user chosen to plot 1D
                ax.plot(dfn_copy[axes[0]],dfn_copy[key],color=colmap(colmap_val),linestyle=line_style,label=label)
                ax.set_xlabel(axes[0])
                ax.set_ylabel(key)
            elif dfn_copy['dfn'].ndim==2: #user chosen to plot 2D

                if all(isinstance(axis,type('_')) for axis in axes): #user has supplied list of chars to denote axes
                    pass
                else: #user has supplied full list of indices to slice DFN -> need to determine convetional axes names  
                    axes=dfn_copy['dfn_index'][np.where([isinstance(axis,slice) for axis in axes])] #do this by assuming that user slices over dimensions they want to plot
                    #the above line works because dfn_index is a numpy array of strings - would break for lists

                if real_scale is True: #set x and y plot limits to real scales
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')

                X=dfn_copy[axes[0]] #make a mesh
                Y=dfn_copy[axes[1]]
                dx,dy=X[1]-X[0],Y[1]-Y[0]
                
                Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres

                if fill:
                    ax.set_facecolor(colmap(np.amin(dfn_copy[key])))
                    mesh=ax.pcolormesh(X,Y,dfn_copy[key],cmap=colmap,vmin=vmin,vmax=vmax)
                    #mesh=ax.contourf(X,Y,dfn_copy[key],levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                    '''for c in mesh.collections: #for use in contourf
                        c.set_edgecolor("face")'''
                else:
                    mesh=ax.contour(X,Y,dfn_copy[key],levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax)
                    if settings.plot_contour_labels:
                        ax.clabel(mesh,inline=1,fontsize=10)

                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')
                ax.set_xlabel(axes[0])
                ax.set_ylabel(axes[1])
                
                if real_scale is True: #set x and y plot limits to real scales
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                if LCFS: #plot plasma boundary
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')
                if gridlines:
                    #get bin edges
                    d_ax_0=dfn_copy[axes[0]][1]-dfn_copy[axes[0]][0]
                    d_ax_1=dfn_copy[axes[1]][1]-dfn_copy[axes[1]][0]
                    axes_0_edges=(dfn_copy[axes[0]][:-1]+dfn_copy[axes[0]][1:])*0.5 #get bin centres of bin centres i.e. bin edges
                    axes_1_edges=(dfn_copy[axes[1]][:-1]+dfn_copy[axes[1]][1:])*0.5
                    axes_0_edges=np.concatenate((axes_0_edges[0]-[d_ax_0],axes_0_edges,axes_0_edges[-1]+[d_ax_0])) #add outermost values to bin edges
                    axes_1_edges=np.concatenate((axes_1_edges[0]-[d_ax_1],axes_1_edges,axes_1_edges[-1]+[d_ax_1]))
                    for line_axis_0,line_axis_1 in zip(axes_0_edges,axes_1_edges):
                        ax.axvline(line_axis_0,color=settings.plot_colour_gridlines) #assume that axis 0 is X axis and axis 1 is y axis
                        ax.axhline(line_axis_1,color=settings.plot_colour_gridlines)

                if ax_flag is True or fig_flag is True: #return the plot object
                    return mesh

            else: #user has not supplied >2D dfn
                print("ERROR: distribution_function.plot() given >2D DFN - please reduce dimensionality")
                return 

        if ax_flag is False and fig_flag is False:
            plt.show()

    def transform(self,axes=['R','Z']):
        """
        transforms and integrates the distribution function according to pre-defined configurations 
        
        args:
            axes - the dimensions over which to transform the DFN to
        notes:
            remember dimensions of unedited dfn are my_dfn['dfn'][P,V/E,V_pitch,R,Z]
            assumes unedited dfn
            assumes the bin widths for a given dimension are constant (keep this as crop may add discontinuities in dimension axes)
            assumes toroidal symmetry (no toroidal dimension in dfn)
            if an array of indices is given, then slice the dfn accordingly and return without any integration
                note for an infinite slice, axes will need to contain slice() objects e.g. axes=[0,0,0,slice(None),slice(None)] for all R,Z values
            dimension P is meaningless in EBASE mode 

        axes options:
            IDFTYP==1
                R,Z - integrate over pitch, gyrophase and velocity [m]^-3
                E,V_pitch - integrate over space and transform to [eV]^-1[dpitch]^-1 
                E - [eV]^-1 
                R - [m]^-3
                Z - [m]^-3
                V_pitch - [dPitch]^-1  
                N - total # 
            list of indices and slices
        """

        dfn=copy.deepcopy(self) #make deep copy here since functions designed to repeatedly take fresh DFNs would otherwise permanently change it


        #general option
        if len(axes)==dfn['dfn'].ndim: #if user supplies all axes then slice WITHOUT integrating
            dfn['dfn']=dfn['dfn'][tuple(axes)]
            #XXX need to then reset dfn['nV'],dfn['R'],dfn['dfn_index'] etc data here?

        else:

            #begin list of specific options

            if dfn['IDFTYP']==1:

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

                    elif axes==['Z']:
                        dfn['dfn']*=dfn['dE']*dfn['dV_pitch'] #integrate
                        for counter in range(3): #sum over gyrophase, pitch and energy
                            dfn['dfn']=np.sum(dfn['dfn'],axis=0)
                        dfn['dfn']=np.sum(dfn['dfn'],axis=0) #sum over R

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
                        print("ERROR: dfn_transform given invalid axes argument: {axes} (ID={ID})".format(axes=str(axes),ID=self.ID))

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

            if dfn['IDFTYP']==3: #assume that each e.g. R dimension has corresponding dR quantity

                #assume rectilinear dimensions
                #XXX do not necessarily want to multiply by Jacobian here, since larger bins will measure larger densities
                #XXX the Jacobian is also completely numerical and not separated into the DFN dimensions 
                #XXX (i.e. like IDFTYP1 is where J in R dir is 2*pi*r, in Z dir is just 1 etc. these are all applicable separately) 
                dfn['dfn']=dfn['dfn']*dfn['J'] #XXX cannot use *= here due to overflow

                if axes==['N']:
                    axes_to_integrate_over=[[index,axis] for index,axis in enumerate(dfn['dfn_index'])]
                else:
                    axes_to_integrate_over=[[index,axis] for index,axis in enumerate(dfn['dfn_index']) if axis not in axes]                        

                for axes in reversed(axes_to_integrate_over):
                    dfn['dfn']*=dfn[f'd{axes[1]}']
                    for quantity in ['dfn','dfn_s']: dfn[quantity]=np.sum(dfn[quantity],axis=axes[0]) #sum
                    for quantity in ['J','J_s']: dfn[quantity]=np.sum(dfn[quantity],axis=axes[0]) #XXX is this correct thing to do?

            if len(axes)!=dfn['dfn'].ndim: #if user has not supplied all axes
                dfn['dfn_index']=np.array(axes)

        return dfn

    def crop(self,inside=True,**kwargs):
        """
        notes:
            warning! to work, dfn_index and 1D dfn axes must accurately reflect dfn which is still stored e.g. dfn['dfn'][r,z] must contain dfn['R'],dfn['Z'] and dfn['dfn_index']=['R','Z']
            always maintains shape of dfn
        args:
            kwargs - axes and their limits 
            inside - toggle whether to crop inside or outside supplied limits if 2D
        usage:
            new_dfn=crop(R=[1]) generates dfn at point closest to R=1
            new_dfn=crop(R=[0,1]) crops dfn between 0<R<1 (sets outside of this=0)
            new_dfn=crop(R=[1,0]) crops dfn between 1<R and R<0 (sets inside of this=0)
        """

        dfn=copy.deepcopy(self)

        if not inside:
            dfn['dfn'][:]=0

        for key,value in kwargs.items():
            if key not in dfn['dfn_index']:
                print("ERROR: crop supplied invalid axis name ({}) - see ['dfn_index'] for possible axes".format(key))    
            else:
                dimension_to_edit=dfn['dfn_index'].tolist().index(key) #figure out which dimension we are cropping over
                if len(value)==2: #user has supplied range - get new indices which satisfy range
                    if value[1]>value[0]: #range is a window - set outside to zero
                        i=np.any([(value[0]>dfn[key]),(dfn[key]>value[1])],axis=0)
                    else: #range is outside a window - set window to 0
                        i=np.all([(value[1]<dfn[key]),(dfn[key]<value[0])],axis=0)
                    
                    dfn['dfn']=np.moveaxis(dfn['dfn'],dimension_to_edit,0) #move desired axis of dfn array to front to crop
                    if inside: #invert if not inside
                        dfn['dfn'][i]=0 #crop dfn
                    else:
                        self['dfn']=np.moveaxis(self['dfn'],dimension_to_edit,0)
                        dfn['dfn'][i]=self['dfn'][i]
                        self['dfn']=np.moveaxis(self['dfn'],0,dimension_to_edit) 
                    dfn['dfn']=np.moveaxis(dfn['dfn'],0,dimension_to_edit) #move desired axis of dfn array back to original position                             
                
                elif len(value)==1: #user has supplied single value for nearest neighbour
                    i=[np.abs(dfn[key]-value[0]).argmin()]
                    dfn[key]=dfn[key][i] #crop 1D arrays accordingly
                    nkey='n{}'.format(key) #reset associated nkey values too e.g. reset nR if cropping R
                    dfn[nkey]=np.array(len(dfn[key]))
                    dfn['dfn']=np.moveaxis(dfn['dfn'],dimension_to_edit,0) #move desired axis of dfn array to front to crop
                    dfn['dfn']=dfn['dfn'][i] #crop dfn
                    dfn['dfn']=np.moveaxis(dfn['dfn'],0,dimension_to_edit) #move desired axis of dfn array back to original position             

        return dfn

#################################

##################################################################

###################################################################################################
