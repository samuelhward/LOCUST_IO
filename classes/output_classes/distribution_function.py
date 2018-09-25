#distribution_function.py

"""
Samuel Ward
01/04/2018
----
class to handle LOCUST distribution function output data
---
usage:
    see README.md for usage

notes:         
---
"""


###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from output_classes import x") but best practice to import whole output_classes module anyway
try:
    import numpy as np
    import copy
    import re
    from scipy.io import  FortranFile
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
    from classes import base_output 
except:
    raise ImportError("ERROR: base_output.py could not be imported!\nreturning\n")
    sys.exit(1) 
try:
    from classes import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)


np.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays
pi=np.pi
e_charge=1.60217662e-19 #define electron charge
mass_neutron=1.674929e-27 #define mass of neutron
amu=1.66053904e-27
mass_deuterium=2.0141017781*amu


################################################################## Distribution_Function functions

def read_distribution_function_LOCUST(filepath,ITER=True,wtot=False,WIPE=False,TEST=False,EBASE=False,dfn_s=True,Jh=True,Jh_s=True,cpu_time=True,**properties):
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

    filename=filepath.split('/')[-1] #infer IDFTYP from first character of file name
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
        print('ERROR: cannot infer IDFTYPE from filename\n')

    file=FortranFile(filepath,'r')
    input_data={} #initialise blank dictionary
    input_data['IDFTYP']=np.array(IDFTYP)

    if IDFTYP==3:

        #32-char checksum (32 bytes in Fortran, so essentially reading 8x32 bit floats)
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum
        
        if EBASE:
            #nEF-1
            input_data['nE']=file.read_ints() 
        else:
            #nV-1
            input_data['nV']=file.read_ints() #V

        #nPP-1
        input_data['nPP']=file.read_ints() #PPh
        #nMU-1
        input_data['nMU']=file.read_ints() #MU
        
        if EBASE:
            #dEFh
            input_data['dEh']=file.read_reals(dtype=np.float32)         
        else:
            #dVFh
            input_data['dV']=file.read_reals(dtype=np.float32) 
        
        #dPPh
        input_data['dPPh']=file.read_reals(dtype=np.float32) #dPPh    = P_phi1h - P_phi0h
        #dMUh
        input_data['dMUh']=file.read_reals(dtype=np.float32)

        if EBASE:
            #EF+dEF/2 (nV long)
            input_data['E']=file.read_reals(dtype=np.float32)
        else:
            #VF+dVF/2 (nV long)
            input_data['V']=file.read_reals(dtype=np.float32)

        #PP+dPP/2 (nPP long)
        input_data['PP']=file.read_reals(dtype=np.float32) #PPhi
        #MU+MU/2 (nMU long)
        input_data['MU']=file.read_reals(dtype=np.float32)
        
        #for now these blocks do the same thing, but useful to leave as separate for now
        if wtot: #cumulative energy inventory (total energy injected so far)
            if EBASE:
                #Fh_norm (nEF by nPP by nMU)
                input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
            else:
                #Fh_norm (nV by nPP by nMU)
                input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
        else:
            #Fh (NOTE check dimensions e.g. nV by nPP by nMU)
            input_data['dfn']=file.read_reals(dtype=np.float32) #Final combined DFn. grid
    
        if Fh_s:
            #Fh_s (nV by nPP by nMU) if Fh_s present
            input_data['dfn_s']=file.read_reals(dtype=np.float32) #Dfn. M.C. error
        if Jh:
            #Jh
            input_data['Jh']=file.read_reals(dtype=np.float64) #Jacobian for IDFTYP=3 
        if Jh_s:
            #Jh_s
            input_data['Jh_s']=file.read_ints(dtype=np.float64) #Jacobian error for IDFTYP=3
            if cpu_time:
                #cpuTime
                input_data['cpu_time']=file.read_reals(dtype=np.float64) #this may not be present in the file
    

        #extra derived data
        #input_data['dfn_index']=np.array(['E','V','P','MU']) #reference for names of each dfn dimension

    else:

        #32-char checksum
        input_data['EQBM_md5']=file.read_reals(dtype=np.float32) #equilibrium checksum
        #zero indicates fast ions only
        input_data['0_1']=file.read_ints()

        if IDFTYP==4:
            #nPSIF-1
            input_data['nPSIF']=file.read_ints() #number of surface contours
            #nPOLF-1
            input_data['nPOLF']=file.read_ints() #number of poloidal cells
        else:
            #nF-1
            input_data['nR']=file.read_ints() #R and Z dimension of the distribution function grid (equal here)
            input_data['nZ']=file.read_ints()

        if IDFTYP==1 or IDFTYP==4:
            #nL-1 
            input_data['nV_pitch']=file.read_ints() #Vphi/V cell boundaries
        else:
            #nPP-1
            input_data['nPP']=file.read_ints() #PPhi          

        if EBASE:
            #nE-1
            input_data['nE']=file.read_ints() #E
        else:
            #nV-1
            input_data['nV']=file.read_ints() #V        
        
        #nP-1
        input_data['nP']=file.read_ints() #poloidal gyro-phase cell boundaries
        
        if ITER: #NOTE these blocks essentially do the same thing right now

            if WIPE:
                input_data['dfn']=[] #Final combined DFn. grid
                input_data['dfn_s']=[] #Dfn. M.C. error
                for line in range(int(input_data['nP'])):
                    input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                for line in range(int(input_data['nP'])):
                    input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                
                input_data['dfn']=np.array(input_data['dfn']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F') 
                input_data['dfn_s']=np.array(input_data['dfn_s']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F')
                input_data['nc']=len(input_data['dfn'])/input_data['nP'] #nV*nV_pitch*nZ*nR (nZ = nR = nF)

            else:
                input_data['dfn']=[] #Final combined DFn. grid  
                input_data['dfn_s']=[] #Dfn. M.C. error
                for line in range(int(input_data['nP'])): #Dfn. M.C. error
                    input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
                for line in range(int(input_data['nP'])):
                    input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long

                input_data['dfn']=np.array(input_data['dfn']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F') 
                input_data['dfn_s']=np.array(input_data['dfn_s']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F')
                input_data['nc']=len(input_data['dfn'])/input_data['nP'] #nV*nV_pitch*nZ*nR (nZ, nR = nF)
        
        else:
            input_data['dfn']=[] #Final combined DFn. grid   
            input_data['dfn_s']=[] #Dfn. M.C. error
            for line in range(int(input_data['nP'])): #Dfn. M.C. error
                input_data['dfn'].extend(file.read_reals(dtype=np.float32)) #nP*nc long
            for line in range(int(input_data['nP'])):
                input_data['dfn_s'].extend(file.read_reals(dtype=np.float32)) #nP*nc long

            input_data['dfn']=np.array(input_data['dfn']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F') 
            input_data['dfn_s']=np.array(input_data['dfn_s']).reshape(int(input_data['nP']),int(input_data['nV']),int(input_data['nV_pitch']),int(input_data['nZ']),int(input_data['nR']),order='F')
            input_data['nc']=len(input_data['dfn'])/input_data['nP'] #nV*nV_pitch*nZ*nR (nZ, nR = nF)
        
        input_data['dfn']=np.swapaxes(input_data['dfn'],3,4) #swap final order to ...r,z - means plotting functions can assume index order x,y
        nEQ=file.read_ints() 
        input_data['nR_1D']=np.array(nEQ[0]) #2D field grid R dimension
        input_data['nZ_1D']=np.array(nEQ[1]) #2D field grid Z dimension
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
            input_data['R']=file.read_reals(dtype=np.float32) #r space of dfn
            #Z+dZ/2 (nF long)
            input_data['Z']=file.read_reals(dtype=np.float32) #z space of dfn

        if IDFTYP==1 or IDFTYP==4:

            #L+dL/2 (nL long)
            input_data['V_pitch']=file.read_reals(dtype=np.float32) #pitch space of dfn

        else:

            #PP+dPP/2 (nPP long)
            input_data['PP']=file.read_reals(dtype=np.float32)

        if EBASE:
            #E+dE/2 (nE long)
            input_data['E']=file.read_reals(dtype=np.float32) #energy space of dfn
        else:
            #V+dV/2 (nV long)
            input_data['V']=file.read_reals(dtype=np.float32) #velocity space of dfn
        
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

        if TEST:

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

            if cpu_time:
                input_data['cpu_time']=file.read_reals(dtype=np.float64)

        #extra derived data
        if input_data['nP']>1:
            input_data['dP']=np.array(input_data['P'][1]-input_data['P'][0]) #special dimension bin width 
        else:
            input_data['dP']=np.array(2.*pi)
        
        if EBASE:
            input_data['V']=np.array(np.sqrt(2.*input_data['E']*e_charge/mass_deuterium))
        else:
            input_data['E']=np.array((0.5*mass_deuterium*input_data['V']**2)/e_charge) #calculate energy [eV]
        
        input_data['dfn_index']=np.array(['P','V','V_pitch','R','Z']) #reference for names of each dfn dimension
   
        input_data['dR']=np.array(input_data['R'][1]-input_data['R'][0]) #R bin width
        input_data['dZ']=np.array(input_data['Z'][1]-input_data['Z'][0]) #Z bin width
        input_data['dV_pitch']=np.array(input_data['V_pitch'][1]-input_data['V_pitch'][0]) #pitch bin width
        input_data['dV']=np.array(input_data['V'][1]-input_data['V'][0]) #velocity bin width
        input_data['dE']=np.array(np.abs(input_data['E'][1]-input_data['E'][0])) #energy bin width

        
    file.close()

    print("finished reading distribution function from LOCUST")
    return input_data



def dump_distribution_function_LOCUST(output_data,filepath,**properties): 
    """
    writes distribution function to LOCUST
    
    notes:
        
    """

    pass



################################################################## Distribution_Function class

class Distribution_Function(base_output.LOCUST_output):
    """
    class describing distribution function output for LOCUST
    
    inheritedfrom LOCUST_output:
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
                P       - special, rare simulation specific (e.g. gyro phase dimension)
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

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading GEQDSKs

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties={**properties}
                self.data=read_distribution_function_LOCUST(self.filepath,**properties) #read the file
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write distribution function to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_distribution_function_LOCUST(self.data,filepath,**properties)
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST)\n")




#################################

##################################################################

###################################################################################################
