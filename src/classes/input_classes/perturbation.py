#perturbation.py
 
"""
Samuel Ward
29/07/2018
----
class to handle LOCUST perturbation input data
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
    import pathlib
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
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## Perturbation read functions
 
def read_perturbation_LOCUST(filepath,**properties):
    """
    reads perturbation stored in LOCUST format

    notes:
        assumes R is slower-varying dimension in file when inferring dimensions
    """

    print("reading LOCUST perturbation")

    with open(filepath,'r') as file:

        #initialise data
        input_data={}
        input_data['R_2D']=[]
        input_data['Z_2D']=[]
        input_data['dB_field_R_real']=[]
        input_data['dB_field_R_imag']=[]
        input_data['dB_field_Z_real']=[]
        input_data['dB_field_Z_imag']=[]
        input_data['dB_field_tor_real']=[]
        input_data['dB_field_tor_imag']=[]

        #read lazily
        for line in file:
            split_line=line.split()
            if len(split_line)==8:
                input_data['R_2D'].append(float(split_line[0]))
                input_data['Z_2D'].append(float(split_line[1]))
                input_data['dB_field_R_real'].append(float(split_line[2]))
                input_data['dB_field_R_imag'].append(float(split_line[3]))
                input_data['dB_field_Z_real'].append(float(split_line[4]))
                input_data['dB_field_Z_imag'].append(float(split_line[5]))
                input_data['dB_field_tor_real'].append(float(split_line[6]))
                input_data['dB_field_tor_imag'].append(float(split_line[7]))

        input_data['R_2D']=np.asarray(input_data['R_2D'])
        input_data['Z_2D']=np.asarray(input_data['Z_2D'])
        input_data['dB_field_R_real']=np.asarray(input_data['dB_field_R_real'])
        input_data['dB_field_R_imag']=np.asarray(input_data['dB_field_R_imag'])
        input_data['dB_field_Z_real']=np.asarray(input_data['dB_field_Z_real'])
        input_data['dB_field_Z_imag']=np.asarray(input_data['dB_field_Z_imag'])
        input_data['dB_field_tor_real']=np.asarray(input_data['dB_field_tor_real'])
        input_data['dB_field_tor_imag']=np.asarray(input_data['dB_field_tor_imag'])
  
        #infer the grid dimensions and axes
        if input_data['Z_2D'][0]==input_data['Z_2D'][1]: #Z is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(Z_dim,R_dim).T
            input_data['Z_2D']=input_data['Z_2D'].reshape(Z_dim,R_dim).T
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['dB_field_R_real']=input_data['dB_field_R_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_R_imag']=input_data['dB_field_R_imag'].reshape(Z_dim,R_dim).T
            input_data['dB_field_Z_real']=input_data['dB_field_Z_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_Z_imag']=input_data['dB_field_Z_imag'].reshape(Z_dim,R_dim).T
            input_data['dB_field_tor_real']=input_data['dB_field_tor_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_tor_imag']=input_data['dB_field_tor_imag'].reshape(Z_dim,R_dim).T
        else: #R is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(R_dim,Z_dim)
            input_data['Z_2D']=input_data['Z_2D'].reshape(R_dim,Z_dim)
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['dB_field_R_real']=input_data['dB_field_R_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_R_imag']=input_data['dB_field_R_imag'].reshape(R_dim,Z_dim)
            input_data['dB_field_Z_real']=input_data['dB_field_Z_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_Z_imag']=input_data['dB_field_Z_imag'].reshape(R_dim,Z_dim)
            input_data['dB_field_tor_real']=input_data['dB_field_tor_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_tor_imag']=input_data['dB_field_tor_imag'].reshape(R_dim,Z_dim)

    print("finished reading LOCUST perturbation")
    
    return input_data

def read_perturbation_LOCUST_field_data(filepath,**properties):
    """
    notes:
        reads from the file_data.out file produced by LOCUST BCHECK mode
    """

    print("reading LOCUST test field data")

    with open(filepath,'r') as file: #open file

        input_data={}

        lines=file.readlines()
        del(lines[0]) #delete headerline

        for quantity in ['R','phi','Z','time','B_field_R','B_field_tor','B_field_Z','dB_field_R','dB_field_tor','dB_field_Z','divB']:
            input_data[quantity]=[]

        for line in lines:
            for counter,quantity in enumerate(['R','phi','Z','time','B_field_R','B_field_tor','B_field_Z','dB_field_R','dB_field_tor','dB_field_Z','divB']):
            
                try:
                    input_data[quantity].append(float(line.split()[counter]))
                except:
                    input_data[quantity].append(0.0)

        for quantity in ['R','phi','Z','time','B_field_R','B_field_tor','B_field_Z','dB_field_R','dB_field_tor','dB_field_Z','divB']:
            input_data[quantity]=np.asarray(input_data[quantity])

    print("reading LOCUST test field data")

    return input_data

def read_perturbation_ASCOT_field_data(filepath,**properties):
    """
    read perturbation field data equivalent from ASCOT input particles file

    notes:
        since the input particles file essentially contains the magnetic field at different points where ions are born, this can be used to calibrate a 3D field
    """

    with open(filepath,'w') as file: #open file

        for line in file:
            if 'Number of particles' in line:
                number_particles=int(line.split()[0])
            if 'Number of different fields' in line:
                number_fields=int(line.split()[0])
                break

        fields=[] #this will hold the names of the quantites stored in the file - in order
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
        ascot_names_1=['Rprt','phiprt' ,'zprt','BR','Bphi','Bz'] #possible ASCOT fields - full orbit
        ascot_names_2=['R','phi' ,'z','BR','Bphi','Bz'] #possible ASCOT fields - guiding centre
        locust_names['R','phi','Z','B_field_R','B_field_tor','B_field_Z'] #corresponding LOCUST_IO fields that we want to retain
        for ascot_name_1,ascot_name_2,locust_io_name in zip(ascot_names_1,ascot_names_2,locust_io_names):
            if ascot_name_1 in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name_1])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])
            elif ascot_name_2 in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name_2])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])                
    
    return input_data

def read_perturbation_IDS_mhd_linear(shot,run,mode_number,**properties):
    """
    reads perturbation from an IMAS mhd_linear IDS

    args:
        shot - IDS shot number
        run -  IDS run number
        mode_number - toroidal mode number of perturbation to read

    notes:
        assumes coordinate 1 = R, coordinate 2 = phi, coordinate 3 = Z
        assumes reading from time_slice 0

    """ 

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_perturbation_IDS_mhd_linear could not import IMAS module!\nreturning\n")
        return

    print("reading perturbation from IMAS mhd_linear IDS")

    input_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    input_IDS.open_env(username,imasdb,'3')
    input_IDS.mhd_linear.get() #open the file and get all the data from it

    mode_index=None
    for counter,mode in enumerate(input_IDS.mhd_linear.time_slice[0].toroidal_mode): #determine where desired harmonic is stored in the IDS
        if mode_number==mode.n_tor:
            mode_index=counter

    if mode_index is None:
        print("ERROR: read_perturbation_IDS_mhd_linear could not find requested mode in IDS (shot - {shot}, run - {run}, n - {mode})!\nreturning\n!".format(shot=shot,run=run,mode=mode_number))
        return

    if input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.grid_type.index!=1:    
        print("WARNING: read_perturbation_IDS_mhd_linear detected non-rectangular grid geometry (shot - {shot}, run - {run}, n - {mode})!".format(shot=shot,run=run,mode=mode_number))

    input_data = {} #initialise blank dictionary to hold the data

    input_data['R_1D']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.grid.dim1)  
    input_data['Z_1D']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.grid.dim2)  
    input_data['dB_field_R_real']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate1.real)
    input_data['dB_field_R_imag']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate1.imaginary)
    input_data['dB_field_Z_real']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate2.real)
    input_data['dB_field_Z_imag']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate2.imaginary)
    input_data['dB_field_tor_real']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate3.real)
    input_data['dB_field_tor_imag']=np.array(input_IDS.mhd_linear.time_slice[0].toroidal_mode[mode_index].plasma.b_field_perturbed.coordinate3.imaginary)

    input_IDS.close()

    print("finished reading perturbation from IMAS mhd_linear IDS")

    return input_data

'''
def read_perturbation_JOREK(filepath,**properties):
    """
    notes:
    """

    try:
        import h5py
    except:
        print("ERROR: read_perturbation_JOREK could not import h5py!\n")
        return

    input_data={} #initialise data dictionary
    file = h5py.File(filepath, 'r')

    input_data['']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['sqrt(PSIn)']) 

    return input_data
'''

def read_perturbation_MARSF(filepath,**properties):
    """
    reads perturbation stored in MARSF format

    notes:
        assumes R is quickly-varying dimension in file when inferring dimensions
    """

    print("reading MARSF perturbation")

    with open(filepath,'r') as file:
                 
        #initialise data
        input_data={}
        input_data['R_2D']=[]
        input_data['Z_2D']=[]
        input_data['dB_field_R_real']=[]
        input_data['dB_field_R_imag']=[]
        input_data['dB_field_Z_real']=[]
        input_data['dB_field_Z_imag']=[]
        input_data['dB_field_tor_real']=[]
        input_data['dB_field_tor_imag']=[]

        #read lazily 
        #skip header lines
        counter=0
        number_headerlines=3
        for line in file:
            counter+=1
            if counter==number_headerlines:
                break

        for line in file:
            split_line=line.split()
            input_data['R_2D'].append(float(split_line[0]))
            input_data['Z_2D'].append(float(split_line[1]))
            input_data['dB_field_R_real'].append(float(split_line[2]))
            input_data['dB_field_R_imag'].append(float(split_line[3]))
            input_data['dB_field_Z_real'].append(float(split_line[4]))
            input_data['dB_field_Z_imag'].append(float(split_line[5]))
            input_data['dB_field_tor_real'].append(float(split_line[6]))
            input_data['dB_field_tor_imag'].append(float(split_line[7]))

        input_data['R_2D']=np.asarray(input_data['R_2D'])
        input_data['Z_2D']=np.asarray(input_data['Z_2D'])
        input_data['dB_field_R_real']=np.asarray(input_data['dB_field_R_real'])
        input_data['dB_field_R_imag']=np.asarray(input_data['dB_field_R_imag'])
        input_data['dB_field_Z_real']=np.asarray(input_data['dB_field_Z_real'])
        input_data['dB_field_Z_imag']=np.asarray(input_data['dB_field_Z_imag'])
        input_data['dB_field_tor_real']=np.asarray(input_data['dB_field_tor_real'])
        input_data['dB_field_tor_imag']=np.asarray(input_data['dB_field_tor_imag'])
    
        #infer the grid dimensions and axes
        if input_data['Z_2D'][0]==input_data['Z_2D'][1]: #Z is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(Z_dim,R_dim).T
            input_data['Z_2D']=input_data['Z_2D'].reshape(Z_dim,R_dim).T
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['dB_field_R_real']=input_data['dB_field_R_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_R_imag']=input_data['dB_field_R_imag'].reshape(Z_dim,R_dim).T
            input_data['dB_field_Z_real']=input_data['dB_field_Z_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_Z_imag']=input_data['dB_field_Z_imag'].reshape(Z_dim,R_dim).T
            input_data['dB_field_tor_real']=input_data['dB_field_tor_real'].reshape(Z_dim,R_dim).T
            input_data['dB_field_tor_imag']=input_data['dB_field_tor_imag'].reshape(Z_dim,R_dim).T
        else: #R is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(R_dim,Z_dim)
            input_data['Z_2D']=input_data['Z_2D'].reshape(R_dim,Z_dim)
            input_data['R_1D']=input_data['R_2D'][0,:].flatten()
            input_data['Z_1D']=input_data['Z_2D'][:,0].flatten()
            input_data['dB_field_R_real']=input_data['dB_field_R_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_R_imag']=input_data['dB_field_R_imag'].reshape(R_dim,Z_dim)
            input_data['dB_field_Z_real']=input_data['dB_field_Z_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_Z_imag']=input_data['dB_field_Z_imag'].reshape(R_dim,Z_dim)
            input_data['dB_field_tor_real']=input_data['dB_field_tor_real'].reshape(R_dim,Z_dim)
            input_data['dB_field_tor_imag']=input_data['dB_field_tor_imag'].reshape(R_dim,Z_dim)

    print("finished reading MARSF perturbation")
    
    return input_data

def read_perturbation_MARSF_bplas(filepath=pathlib.Path(''),response=True,ideal=False,phase=0,bcentr=1.75660107,rmaxis=1.70210874,nR_1D=400,nZ_1D=600):
    """
    read perturbation bplas files produced by MARSF for individual harmonics and coil sets 
    
    args:
       filepath - path to files in input_files/
       response - toggle whether vacuum field or with plasma response i.e. bplas_vac_upper/lower or bplas_ideal/resist_resp_upper/lower
       ideal - toggle whether resistive or ideal bplasma i.e. bplas_ideal_resp_upper/lower or bplas_resist_resp_upper/lower
       phase - phase shift between upper and lower rows (applied to upper coils) [radians]
       bcentr - vacuum toroidal magnetic field at rcentr
       rmaxis - R at magnetic axis (O-point)
       nR_1D - number of points in R
       nZ_1D - number of points in Z
    notes:
       adapted from David Ryan's scripts david.ryan@ukaea.uk
       response overrides ideal toggle setting
       assumes following default filenames within filepath directory:
          rmzm_geom
          rmzm_pest
          profeq
          bplas_vac_upper/bplas_ideal_resp_upper/bplas_resist_resp_upper
          bplas_vac_lower/bplas_ideal_resp_lower/bplas_resist_resp_lower
       e.g. of file name from MARS-F for individual harmonic=bplas_resist_resp_lower
       reading this way allows one to rotate coil sets with respect to each other
    """

    print("reading MARSF_bplas perturbation")

    class rzcoords():
 
        def __init__(self, path, nchi):
            rmzm=scipy.loadtxt(path)
            Nm0=int(rmzm[0,0]) #Num. poloidal harmonics for equilibrium quantities (not necessarily same as for perturbation quantities, but should be).
            Ns_plas=int(rmzm[0,1]) #Num. radial points in plasma
            Ns_vac=int(rmzm[0,2]) #Num. radial points in vacuum
            Ns=Ns_plas+Ns_vac
            R0EXP=rmzm[0,3]
            B0EXP=rmzm[1,3]
            s=rmzm[1:Ns+1, 0]
            RM=rmzm[Ns+1:,0]+1j*rmzm[Ns+1:,1]
            ZM=rmzm[Ns+1:,2]+1j*rmzm[Ns+1:,3]
            RM=RM.reshape((Nm0, Ns))
            RM=scipy.transpose(RM)
            ZM=ZM.reshape((Nm0, Ns))
            ZM=scipy.transpose(ZM)
            RM[:,1:]=2*RM[:,1:]
            ZM[:,1:]=2*ZM[:,1:]
            
            m=scipy.arange(0,Nm0,1)
            chi=scipy.linspace(-scipy.pi, scipy.pi,nchi)
            expmchi=scipy.exp(scipy.tensordot(m,chi,0)*1j)
            R=scipy.dot(RM[:,:Nm0],expmchi)
            Z=scipy.dot(ZM[:,:Nm0],expmchi)
            
            self.R=scipy.array(R.real) #R coordinates
            self.Z=scipy.array(Z.real) #Z coordinates
            self.Nchi=nchi
            self.Nm0=Nm0
            self.Ns_plas=Ns_plas       #number is s points in plasma
            self.Ns_vac=Ns_vac         #number of s points in vacuum
            self.Ns=Ns                 #total number of s points
            self.R0EXP=R0EXP           #normalisation length
            self.B0EXP=B0EXP           #normalisation magnetic field
            self.m=m                   #equilibrium poloidal harmonics
            self.chi=chi               #poloidal angle coordinate
            self.s=s                   #radial coordinate=sqrt(psi_pol)
 
    class jacobian():
 
        #decide what stuff is needed later, and add a self. in front of it.
        def __init__(self, rz):
            if not isinstance(rz, rzcoords):
                print("read_perturbation_MARSF_bplas.jacobian - must pass in coordinate system of type plotting_base.rzcoords")
                return

            self.Ns=rz.Ns #Used in jacobian.plot(), so pass in from rzcoords
            
            self.dRds=scipy.copy(rz.R) 
            self.dZds=scipy.copy(rz.R)
            self.dRdchi=scipy.copy(rz.R)
            self.dZdchi=scipy.copy(rz.R)
            self.jacobian=scipy.copy(rz.R)
    
            #this is for the vacuum region. these are overwritten for the plasma region
            #just having a number to denote the boundary index might be simpler in the future.
            #Vac_start variable should be all that's needed. II is way too complicated
            II_start=int(rz.Ns_plas)-1; II_end=len(rz.R[:,0])
            II2_start=int(rz.Ns_plas); II2_end=len(rz.R[:,0])
           
            s0=scipy.copy(rz.s[II_start:II_end]); R0=scipy.copy(rz.R[II_start:II_end, :])
            chi0=scipy.squeeze(scipy.copy(scipy.array(rz.chi))); Z0=scipy.copy(rz.Z[II_start:II_end, :])
           
            hs=0.5*(s0[1:]-s0[:-1]).min(); hs=min(hs,  2e-5)
            hchi=0.5*(chi0[1:]-chi0[:-1]).min(); hchi=min(hchi,  1e-4)
            s1=s0-hs; s2=s0+hs
            chi1=chi0-hchi;  chi2=chi0+hchi
           
            #compute dR/ds using R(s,chi)
            R1=scipy.zeros(scipy.shape(R0))
            R2=scipy.zeros(scipy.shape(R0))
            for i in range(rz.Nchi):
               R1[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s1[0], s0[-1]])(s1)
               R2[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s0[0], s2[-1]])(s2)
            self.dRds[II_start:II_end,:]=(R2-R1)/(2*hs)
    
            #compute dZ/ds using Z(s,chi) 
            Z1=scipy.zeros(scipy.shape(Z0))
            Z2=scipy.zeros(scipy.shape(Z0))
            for i in range(rz.Nchi):
               Z1[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s1[0], s0[-1]])(s1)
               Z2[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s0[0], s2[-1]])(s2)
            self.dZds[II_start:II_end,:]=(Z2-Z1)/(2*hs)
    
            # compute dR/dchi using R(s,chi) 
            R1=scipy.zeros(scipy.shape(R0))
            R2=scipy.zeros(scipy.shape(R0))
            for i in range(int(rz.Ns_vac)+1):
               R1[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
               R2[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi0[0], chi2[-1]])(chi2) 
               self.dRdchi[i+II_start,:]=(R2[i,:]-R1[i,:])/(2*hchi) 
       
            #compute dZ/dchi using Z(s,chi) 
            Z1=scipy.zeros(scipy.shape(Z0))
            Z2=scipy.zeros(scipy.shape(Z0))
            for i in range(int(rz.Ns_vac)+1):
               Z1[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
               Z2[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
            self.dZdchi[II_start:II_end,:]=(Z2-Z1)/(2*hchi)
       
            #Now do same calculations for plasma region
            II_start=0; II_end=rz.Ns_plas;
       
            s0=scipy.copy(rz.s[II_start:II_end]); R0=scipy.copy(rz.R[II_start:II_end, :])
            chi0=scipy.squeeze(scipy.copy(scipy.array(rz.chi))); Z0=scipy.copy(rz.Z[II_start:II_end, :])
       
            hs=0.5*(s0[1:]-s0[:-1]).min(); hs=min(hs,  2e-5)
            hchi=0.5*(chi0[1:]-chi0[:-1]).min(); hchi=min(hchi,  1e-4)
            s1=s0-hs; s2=s0+hs
            chi1=chi0-hchi;  chi2=chi0+hchi
       
            #compute dR/ds using R(s,chi) 
            R1=scipy.zeros(scipy.shape(R0))
            R2=scipy.zeros(scipy.shape(R0))
            for i in range(rz.Nchi):
                R1[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s1[0], s0[-1]])(s1)
                R2[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s0[0], s2[-1]])(s2)
            self.dRds[II_start:II_end,:]=(R2-R1)/(2*hs)
    
            #compute dZ/ds using Z(s,chi) 
            Z1=scipy.zeros(scipy.shape(Z0))
            Z2=scipy.zeros(scipy.shape(Z0))
            for i in range(rz.Nchi):
                Z1[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s1[0], s0[-1]])(s1)
                Z2[:,i]=scipy.interpolate.InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s0[0], s2[-1]])(s2)
                self.dZds[:II_end,i]=(Z2[:,i]-Z1[:,i])/(2*hs)
    
            # compute dR/dchi using R(s,chi)
            R1=scipy.zeros(scipy.shape(R0))
            R2=scipy.zeros(scipy.shape(R0))
            for i in range(int(rz.Ns_plas)):
               R1[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
               R2[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
            self.dRdchi[II_start:II_end,:]=(R2-R1)/(2*hchi)
    
            #compute dZ/dchi using Z(s,chi) 
            Z1=scipy.zeros(scipy.shape(Z0))
            Z2=scipy.zeros(scipy.shape(Z0))
            for i in range(int(rz.Ns_plas)):
               Z1[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
               Z2[i,:]=scipy.interpolate.InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
            self.dZdchi[II_start:II_end,:]=(Z2-Z1)/(2*hchi)
    
            G11=scipy.square(self.dRds)+scipy.square(self.dZds)
            G12=scipy.multiply(self.dRds, self.dRdchi)+scipy.multiply(self.dZds, self.dZdchi)
            G22=scipy.square(self.dRdchi)+scipy.square(self.dZdchi)
            G22[0,:]=G22[1,:]
            G33=scipy.square(rz.R)
    
            #Metrics elements
            self.G11=G11
            self.G12=G12
            self.G22=G22
            self.G33=G33
    
            self.jacobian=(-self.dRdchi*self.dZds+self.dRds*self.dZdchi)*rz.R
            self.jacobian[0,:]=self.jacobian[1,:]
    
    class bplasma():
 
        def __init__(self, path, rz, jc):
        
            self.path=path
            bplasma=scipy.loadtxt(self.path)
            Nm1=int(bplasma[0,0]) #Number of perturbation poloidal harmonics (should be same as equilibrium harmonics)
            self.bm1=bplasma[Nm1+1:, 0]+1j*bplasma[Nm1+1:, 1]
            self.bm2=bplasma[Nm1+1:, 2]+1j*bplasma[Nm1+1:, 3]
            self.bm3=bplasma[Nm1+1:, 4]+1j*bplasma[Nm1+1:, 5]
        
            self.bm1=self.bm1.reshape((Nm1, rz.Ns))
            self.bm2=self.bm2.reshape((Nm1, rz.Ns))
            self.bm3=self.bm3.reshape((Nm1, rz.Ns))
        
            #bm2 and bm3 are defined at half int points. 
            #3 ways to recompute at int points, see MacReadBPLASMA.m, lines 109-124
            #For now, simplest implemented. Assume spline_B23==2.
            for i in range(len(self.bm2[:,1])):
                self.bm2[i, 1:]=self.bm2[i, :-1]
                self.bm3[i, 1:]=self.bm3[i, :-1]
        
            m=scipy.array(bplasma[1:Nm1+1,0])
        
            expmchi=scipy.exp(scipy.tensordot(m,rz.chi,0)*1j)
        
            self.b1=scipy.dot(self.bm1.T,expmchi)
            self.b2=scipy.dot(self.bm2.T,expmchi)
            self.b3=scipy.dot(self.bm3.T,expmchi)
        
            self.bn=self.b1/scipy.sqrt(jc.G22*jc.G33)
        
            self.m=m 
        
            self.Br=scipy.divide(scipy.multiply(self.b1, jc.dRds)+scipy.multiply(self.b2, jc.dRdchi), jc.jacobian)
            self.Bz=scipy.divide(scipy.multiply(self.b1, jc.dZds)+scipy.multiply(self.b2,jc.dZdchi), jc.jacobian)
            self.Bphi=scipy.divide(scipy.multiply(self.b3, rz.R), jc.jacobian)
        
            self.Br[0,:]=self.Br[1,:]
            self.Bz[0,:]=self.Bz[1,:]
            self.Bphi[0:2,:]=self.Bphi[3,:]
        
            self.AbsB=scipy.sqrt(scipy.square(scipy.absolute(self.Br))+scipy.square(scipy.absolute(self.Bz))+scipy.square(scipy.absolute(self.Bphi)))
        
            self.Nm1=Nm1
        
            self.rz=rz
            self.jc=jc
       
    #start of function
 
    import scipy
    import scipy.interpolate

    rmzm_geom_path=filepath / 'rmzm_geom'
    rmzm_pest_path=filepath / 'rmzm_pest'
    profeq_path=filepath / 'profeq'

    if response:
        if ideal: #ideal plasma response
          bplas_u_path=filepath / 'bplas_ideal_resp_upper'
          bplas_l_path=filepath / 'bplas_ideal_resp_lower'
        else: #resistive plasma response
            bplas_u_path=filepath / 'bplas_resist_resp_upper'
            bplas_l_path=filepath / 'bplas_resist_resp_lower'
    else:
        bplas_u_path=filepath / 'bplas_vac_upper'
        bplas_l_path=filepath / 'bplas_vac_lower'
 
    #make ascot input from B field
    nchi=2400
  
    R_min=0.5
    R_max=2.5
  
    Z_min=-1.5
    Z_max=1.5
  
    numR=nR_1D
    numZ=nZ_1D
  
    rz_geom=rzcoords(rmzm_geom_path, nchi)
    jc_geom=jacobian(rz_geom)
    bplas_u_geom=bplasma(bplas_u_path, rz_geom,jc_geom)
    bplas_l_geom=bplasma(bplas_l_path, rz_geom,jc_geom)
  
    BN=(bplas_u_geom.bn*scipy.exp(1j*phase)+bplas_l_geom.bn)*bcentr
    BR=(bplas_u_geom.Br*scipy.exp(1j*phase)+bplas_l_geom.Br)*bcentr
    BZ=(bplas_u_geom.Bz*scipy.exp(1j*phase)+bplas_l_geom.Bz)*bcentr
    BP=(bplas_u_geom.Bphi*scipy.exp(1j*phase)+bplas_l_geom.Bphi)*bcentr
  
    B=scipy.sqrt(scipy.square(scipy.absolute(BR))+scipy.square(scipy.absolute(BZ))+scipy.square(scipy.absolute(BP)))
  
    R1=rz_geom.R*rmaxis
    Z1=rz_geom.Z*rmaxis
  
    B1=B
    BR1=BR
    BZ1=BZ
    BP1=BP
  
    R_rect=scipy.linspace(R_min, R_max, numR)
    Z_rect=scipy.linspace(Z_min, Z_max, numZ)
  
    R_grid, Z_grid=scipy.meshgrid(R_rect,Z_rect)
  
    BR_rect=scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BR1.ravel(), (R_grid, Z_grid), method='linear')
    BR_rect=BR_rect.reshape((numZ,numR))
  
    BZ_rect=scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BZ1.ravel(), (R_grid, Z_grid), method='linear')
    BZ_rect=BZ_rect.reshape((numZ,numR))
  
    BP_rect=scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BP1.ravel(), (R_grid, Z_grid), method='linear')
    BP_rect=BP_rect.reshape((numZ,numR))
  
    B_rect=scipy.sqrt(scipy.absolute(BR_rect)**2+scipy.absolute(BZ_rect)**2+scipy.absolute(BP_rect)**2)
 
    input_data={}
    input_data['R_2D']=np.array(R_grid,ndmin=2).swapaxes(0,1)
    input_data['Z_2D']=np.array(Z_grid,ndmin=2).swapaxes(0,1)
    input_data['R_1D']=np.squeeze(np.array(input_data['R_2D'][:,0],ndmin=2))
    input_data['Z_1D']=np.squeeze(np.array(input_data['Z_2D'][0,:],ndmin=2))
    input_data['dB_field_R_real']=np.array(BR_rect.real,ndmin=2).swapaxes(0,1)
    input_data['dB_field_R_imag']=np.array(BR_rect.imag,ndmin=2).swapaxes(0,1)
    input_data['dB_field_Z_real']=np.array(BZ_rect.real,ndmin=2).swapaxes(0,1)
    input_data['dB_field_Z_imag']=np.array(BZ_rect.imag,ndmin=2).swapaxes(0,1)
    input_data['dB_field_tor_real']=np.array(BP_rect.real,ndmin=2).swapaxes(0,1)
    input_data['dB_field_tor_imag']=np.array(BP_rect.imag,ndmin=2).swapaxes(0,1)

    print("finished reading MARSF_bplas perturbation")

    return input_data
    
################################################################## Perturbation write functions
 
def dump_perturbation_LOCUST(output_data,filepath,**properties):
    """
    writes perturbation to LOCUST format

    notes:
        dumps data with quickly-varying Z like usual LOCUST input
    """
 
    print("writing LOCUST perturbation")

    output_data['Z_2D'],output_data['R_2D']=np.meshgrid(output_data['Z_1D'],output_data['R_1D'])

    with open(filepath,'w') as file: #open file
        
        quantities=['R_2D','Z_2D','dB_field_R_real','dB_field_R_imag','dB_field_Z_real','dB_field_Z_imag','dB_field_tor_real','dB_field_tor_imag']

        for row in np.array([output_data[quantity].flatten() for quantity in quantities]).T:
            line=''
            for number in row:
                line+=processing.utils.fortran_string(number_out=number,length=18,decimals=10,exponential=True)
            line+=' \n' 
            file.write(line)

    print("finished writing LOCUST perturbation")

def dump_perturbation_point_data_LOCUST(output_data,filepath='point_data.inp',BCHECK=1,**properties):
    """
    generates the point_data.inp file for checking magnetic perturbations using LOCUST -DBCHECK

    args:
        BCHECK - coordinate format setting for LOCUST field checking (1=RPhiZ,2=XYZ)  
    notes:
        uses R_point_data, phi_point_data, Z_point_data and time_point_data arrays stored in perturbation class 
    """

    print("writing point_inp.dat test points")

    with open(filepath,'w') as file: #open file

        if BCHECK==1:
            for R,Phi,Z,time in zip(output_data['R_point_data'],output_data['phi_point_data'],output_data['Z_point_data'],output_data['time_point_data']):
                line=' '
                line+=processing.utils.fortran_string(R,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(Phi,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(Z,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(time,11,6,exponential=False)
                line+='  '
                file.write('{}\n'.format(line))

        elif BCHECK==2:
            for X,Y,Z,time in zip(output_data['X_point_data'],output_data['Y_point_data'],output_data['Z_point_data'],output_data['time_point_data']):
                line=' '
                line+=processing.utils.fortran_string(X,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(Y,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(Z,11,6,exponential=False)
                line+=' '
                line+=processing.utils.fortran_string(time,11,6,exponential=False)
                line+='  '
                file.write('{}\n'.format(line))

    print("finished writing point_inp.dat test points")

def dump_perturbation_IDS_mhd_linear(ID,output_data,shot,run,mode_number,**properties):
    """
    dumps perturbation to an IMAS mhd_linear IDS

    args:
        shot - IDS shot number
        run -  IDS run number
        mode_number - toroidal mode number of perturbation to write

    notes:
        assumes coordinate 1 = R, coordinate 2 = phi, coordinate 3 = Z
        assumes writing to time_slice 0

    """ 

    try:
        import imas 
    except:
        raise ImportError("ERROR: dump_perturbation_IDS_mhd_linear could not import IMAS module!\nreturning\n")
        return

    print("dumping perturbation from IMAS mhd_linear IDS")

    output_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    output_IDS.open_env(username,imasdb,'3')
    output_IDS.mhd_linear.get() #open the file and get all the data from it

    output_IDS.mhd_linear.ids_properties.comment=ID #write out identification
    output_IDS.mhd_linear.code.name="LOCUST_IO"
    if settings.commit_hash_default_LOCUST_IO: output_IDS.mhd_linear.code.commit=str(settings.commit_hash_default_LOCUST_IO)
    output_IDS.mhd_linear.code.version=support.LOCUST_IO_version
    output_IDS.mhd_linear.ids_properties.homogeneous_time=1   #must set homogeneous_time variable
    output_IDS.mhd_linear.time=np.array([0.0]) #define timebase
    output_IDS.mhd_linear.time_slice.resize(1) #create first time slice

    output_IDS.mhd_linear.time_slice[0].toroidal_mode.resize(len(output_IDS.mhd_linear.time_slice[0].toroidal_mode)+1) #add a perturbation

    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].n_tor=mode_number #set mode number 
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.grid_type.index=1 #define geometry
   
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.grid.dim1=output_data['R_1D']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.grid.dim2=output_data['Z_1D']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate1.real=output_data['dB_field_R_real']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate1.imaginary=output_data['dB_field_R_imag']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate2.real=output_data['dB_field_Z_real']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate2.imaginary=output_data['dB_field_Z_imag']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate3.real=output_data['dB_field_tor_real']
    output_IDS.mhd_linear.time_slice[0].toroidal_mode[-1].plasma.b_field_perturbed.coordinate3.imaginary=output_data['dB_field_tor_imag']

    output_IDS.mhd_linear.put()
    output_IDS.close()

    print("finished dumping perturbation from IMAS mhd_linear IDS")

def dump_perturbation_POCA(output_data,filepath,**properties):
    """
    writes perturbation to POCA format

    notes:
        dumps data with quickly-varying Z like usual LOCUST input
    """
 
    print("writing POCA perturbation")

    with open(filepath,'w') as file: #open file

        file.write(' IPEC_BRZPHI: Total perturbed field\n')
        file.write('\n')
        file.write('   nr =   {nr}  nz =   {nz}\n'.format(nr=len(output_data['R_1D']),nz=len(output_data['Z_1D'])))
        file.write('\n')
        file.write('  l               r               z       real(b_r)       imag(b_r)       real(b_z)       imag(b_z)     real(b_phi)     imag(b_phi) \n')

        quantities=['R_2D','Z_2D','dB_field_R_real','dB_field_R_imag','dB_field_Z_real','dB_field_Z_imag','dB_field_tor_real','dB_field_tor_imag']

        for row in np.array([output_data[quantity].flatten() for quantity in quantities]).T:
            line='  '
            line+='0' #index indicating whether the grid is inside or outside of separatrix - set just to 1 here
            for number in row:
                line+=processing.utils.fortran_string(number_out=number,length=16,decimals=8,exponential=True)
            line+=' \n' 
            file.write(line)

    print("finished writing POCA perturbation")

################################################################## perturbation class
 
class Perturbation(classes.base_input.LOCUST_input):
    """
    class describing magnetic field perturbation for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'perturbation'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in input_files folder
        self.filepath               full path of file in input_files folder  
        self.shot                   shot number
        self.run                    run number
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in input_files folder
 
    notes:
        most methods adapted for fourier-expanded method, for a field expressed as a sum over poloidal harmonics for each toroidal harmonic
    """
 
    LOCUST_input_type='perturbation'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,mode_number=None,**properties):
        """
        read perturbation from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                if mode_number is not None: self.mode_number=int(mode_number)
                self.properties={**properties}
                self.data=read_perturbation_LOCUST(self.filepath,**properties)

        elif data_format=='LOCUST_field_data': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_field_data - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                if mode_number is not None: self.mode_number=int(mode_number)
                self.properties={**properties}
                self.data=read_perturbation_LOCUST_field_data(self.filepath,**properties)

        elif data_format=='ASCOT_field_data': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_field_data - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                if mode_number is not None: self.mode_number=int(mode_number)
                self.properties={**properties}
                self.data=read_perturbation_ASCOT_field_data(self.filepath,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from mhd_linear IDS - shot, run and mode_number required\n".format(self.ID),shot,run,mode_number):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                if mode_number is not None: self.mode_number=int(mode_number)
                self.properties={**properties}
                self.data=read_perturbation_IDS_mhd_linear(self.shot,self.run,self.mode_number,**properties)

        elif data_format=='MARSF': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from MARSF - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                self.mode_number=mode_number
                self.properties={**properties}
                self.data=read_perturbation_MARSF(self.filepath,**properties)

        elif data_format=='MARSF_bplas': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from MARSF_bplas - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files / filename
                if mode_number is not None: self.mode_number=int(mode_number)
                self.properties={**properties}
                self.data=read_perturbation_MARSF_bplas(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/LOCUST_field_data/ASCOT_field_data/IDS/MARSF/MARSF_bplas)\n".format(self.ID))            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,BCHECK=1,**properties):
        """
        write perturbation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run (ID = {})".format(self.ID))
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_perturbation_LOCUST(self.data,filepath,**properties)

        elif data_format=='point_data':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to point_data.inp - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_perturbation_point_data_LOCUST(self.data,filepath,BCHECK,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to IDS - shot and run required\n".format(self.ID),shot,run):
                dump_perturbation_IDS_mhd_linear(self.ID,self.data,shot,run,self.mode_number,**properties)

        elif data_format=='POCA':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to POCA - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files / filename
                dump_perturbation_POCA(self.data,filepath,**properties)

        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/point_data/IDS/POCA)\n".format(self.ID))

    def plot(self,key='dB_field_R_real',axes=['R','Z'],LCFS=False,limiters=False,number_bins=20,fill=True,vminmax=None,i3dr=-1,phase=0.,colmap=cmap_default,colmap_val=np.random.uniform(),gridlines=False,ax=False,fig=False):
        """
        plots a perturbation
        
        notes:
            
        args:
            key - selects which data in perturbation to plot
            axes - list of strings specifying which axes should be plotted (current options only XY or RZ)
            LCFS - toggles plasma boundary on/off in 2D plots (requires equilibrium argument)
            limiters - object which contains limiter data rlim and zlim
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            vminmax - set mesh Vmin/Vmax values
            i3dr - flip definition of phi (+1 anti-clockwise, -1 clockwise)
            phase - global field phase shift (of field origin with respect to locust origin) (radians, anti-clockwise)
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            gridlines - toggle gridlines on plot
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,number_bins,vminmax,i3dr,phase,colmap_val=run_scripts.utils.literal_eval(axes,number_bins,vminmax,i3dr,phase,colmap_val)

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
        if key in self.data.keys():
            if self[key].ndim==0:
                print([key])
                return
        
        #>0D data is plottable
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            polar=True if axes==['X','Y'] else False
            ax = fig.add_subplot(111,polar=polar)
        
        ax.grid(False) if gridlines is False else ax.grid(True)
        if gridlines is False: ax.set_yticklabels([])

        ax.set_title(self.ID)

        #1D data
        if key in self.data.keys():
            if self[key].ndim==1:
                ax.plot(self[key],color=colmap(colmap_val))
                ax.set_ylabel(key)

        #2D data
        if len(axes)==2:
            if axes==['R','Z']:

                R=self['R_1D'] #make a mesh
                Z=self['Z_1D'] 
                Z,R=np.meshgrid(Z,R) #swap since things are defined r,z 
                values=self[key] #2D array (nR_1D,nZ_1D) of poloidal flux
     
                if vminmax:
                    vmin=vminmax[0]
                    vmax=vminmax[1]
                else:
                    vmin=np.amin(values)
                    vmax=np.amax(values)

                #2D plot
                if fill is True:
                    mesh=ax.contourf(R,Z,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    for c in mesh.collections: #for use in contourf
                        c.set_edgecolor("face")
                else:
                    mesh=ax.contour(R,Z,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    if plot_contour_labels:
                        ax.clabel(mesh,inline=1,fontsize=10)
                    
                #mesh=ax.pcolormesh(R,Z,values,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(values),vmax=np.amax(values))

                #3D plot
                #ax=ax.axes(projection='3d')
                #ax.view_init(elev=90, azim=None) #rotate the camera
                #ax.plot_surface(R,Z,values,rstride=1,cstride=1,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(values),vmax=np.amax(values))
                
                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')
                ax.set_aspect('equal')
                ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
                ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D']))
                ax.set_xlabel('R [m]')
                ax.set_ylabel('Z [m]')

                if LCFS:
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],plot_style_LCFS) 
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],plot_style_limiters) 

            elif axes==['X','Y']:

                R=self['R_1D'] #make a mesh
                phi=np.linspace(0.,2.*np.pi*self.mode_number,100) 
                nR,nphi=len(R),len(phi)
                phi,R=np.meshgrid(phi,R)
                R_flat,phi_flat=R.flatten(),phi.flatten()
                Z_flat=np.zeros(len(phi_flat))

                dB_R,dB_tor,dB_Z=self.evaluate(R_flat,phi_flat,Z_flat,mode_number=self.mode_number,i3dr=i3dr,phase=phase)

                if key=='dB_field_R':
                    values=dB_R
                elif key=='dB_field_tor':
                    values=dB_tor
                elif key=='dB_field_Z':
                    values=dB_Z
                elif key=='dB_field':
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                else:
                    values=self[key]
                
                values=values.reshape(nR,nphi)

                if vminmax:
                    vmin=vminmax[0]
                    vmax=vminmax[1]
                else:
                    vmin=np.amin(values)
                    vmax=np.amax(values)

                #2D plot
                if fill is True:
                    mesh=ax.contourf(phi,R,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    for c in mesh.collections: #for use in contourf
                        c.set_edgecolor("face")
                else:
                    mesh=ax.contour(phi,R,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    if plot_contour_labels:
                        ax.clabel(mesh,inline=1,fontsize=10)
                    
                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')

                if LCFS: #plot plasma boundary
                    plasma_max_R=np.max(LCFS['lcfs_r'])
                    plasma_min_R=np.min(LCFS['lcfs_r'])
                    ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(plasma_max_R,100),plot_style_LCFS)
                    ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(plasma_min_R,100),plot_style_LCFS)
                    ax.set_rmax(1.1*plasma_max_R)
 
                if limiters: #add boundaries if desired
                    limiters_max_R=np.max(limiters['rlim'])
                    limiters_min_R=np.min(limiters['rlim'])
                    ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(limiters_max_R,100),plot_style_limiters)
                    ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(limiters_min_R,100),plot_style_limiters)
                    ax.set_rmax(1.1*limiters_max_R)

                ax.set_rmin(0.0)

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()


    def evaluate(self,R,phi,Z,mode_number,i3dr=1,phase=0):
        """
        returns the three components of perturbation field at a point in the plasma 
        
        args:
            R - list of R coordinates
            phi - list of phi coordinates
            Z - list of Z coordinates
            mode_number - mode number of this toroidal harmonic
            i3dr - flip definition of phi (+1 anti-clockwise, -1 clockwise)
            phase - global field phase shift (of field origin with respect to locust origin) (radians, anti-clockwise)
        notes:

        usage:
            dB_R,dB_tor,dB_Z=my_equilibrium.B_calc_point(R=[1,2,3],phi=[0,0,0],Z=[1,2,3])
        """

        print("perturbation_calc_point generating B_field interpolators")
        dB_field_R_real_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_R_real']) #construct interpolators here
        dB_field_R_imag_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_R_imag']) 
        dB_field_tor_real_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_tor_real'])
        dB_field_tor_imag_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_tor_imag'])
        dB_field_Z_real_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_Z_real'])
        dB_field_Z_imag_interpolator=processing.utils.interpolate_2D(self['R_1D'],self['Z_1D'],self['dB_field_Z_imag'])
        print("perturbation_calc_point finished generating B_field interpolators")

        dB_R=[]
        dB_tor=[]
        dB_Z=[]  

        for R_point,phi_point,Z_point in zip(R,phi,Z):
            dB_field_R_real=float(dB_field_R_real_interpolator(R_point,Z_point))
            dB_field_R_imag=float(dB_field_R_imag_interpolator(R_point,Z_point))
            dB_field_tor_real=float(dB_field_tor_real_interpolator(R_point,Z_point))
            dB_field_tor_imag=float(dB_field_tor_imag_interpolator(R_point,Z_point))
            dB_field_Z_real=float(dB_field_Z_real_interpolator(R_point,Z_point))
            dB_field_Z_imag=float(dB_field_Z_imag_interpolator(R_point,Z_point))

            cp = np.cos(mode_number*phi_point*i3dr-phase)
            sp = np.sin(mode_number*phi_point*i3dr-phase)

            dB_R.extend([dB_field_R_real*cp - dB_field_R_imag*sp])
            dB_tor.extend([(dB_field_tor_real*cp - dB_field_tor_imag*sp)*i3dr])
            dB_Z.extend([dB_field_Z_real*cp - dB_field_Z_imag*sp])

        dB_R=np.asarray(dB_R)
        dB_tor=np.asarray(dB_tor)
        dB_Z=np.asarray(dB_Z)

        return dB_R,dB_tor,dB_Z

    def plot_components(self,R,Z,phi,phase=0,i3dr=1,LCFS=False,limiters=False,number_bins=50,vminmax=None,absolute=False,colmap=cmap_default,colmap_val=np.random.uniform(),ax_array=False,fig=False):
        """
        generates plot of perturbation components for field checking

        args:
            R - R coordinate at which 3D field is expanded toroidally
            Z - R coordinate at which 3D field is expanded toroidally
            phi - toroidal point at which to display components 
            phase - global field phase shift (of field origin with respect to locust origin) (radians, anti-clockwise)
            i3dr - flip definition of phi (+1 anti-clockwise, -1 clockwise)
            LCFS - toggles plasma boundary on/off in 2D plots
            limiters - toggles limiters on/off in 2D plots
            number_bins - set number of bins or levels
            vminmax - set mesh Vmin/Vmax values 
            absolute - plot absolute value of perturbation
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            ax_array - take input array of axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        notes:
            user must either supply both fig and ax_array or none
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        R,Z,phi,phase,i3dr,number_bins,vminmax,colmap_val=run_scripts.utils.literal_eval(R,Z,phi,phase,i3dr,number_bins,vminmax,colmap_val)

        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib import cm
        colmap=matplotlib.cm.get_cmap('plasma') #set default colourmap

        if ax_array is False:
            ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
        else:
            ax_flag=True

        if fig is False:
            fig_flag=False
        else:
            fig_flag=True

        if fig_flag is False and ax_flag is False:
            fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)
        else:
            ax1,ax2,ax3,ax4=ax_array
            ax4.cla()

        number_points_toroidal=10000
        R_toroidal=np.full(number_points_toroidal,R) #these are the points to evaluate the field at when we take a single point and expand toroidally
        Z_toroidal=np.full(number_points_toroidal,Z)
        phi_toroidal=np.linspace(0.,2.*constants.pi,number_points_toroidal)

        R_poloidal_dim=len(self['R_1D'])
        Z_poloidal_dim=len(self['Z_1D'])
        R_poloidal=np.linspace(np.min(self['R_1D']),np.max(self['R_1D']),R_poloidal_dim) #these are the points to evaluate the field at when we look at a single plane
        Z_poloidal=np.linspace(np.min(self['Z_1D']),np.max(self['Z_1D']),Z_poloidal_dim)
        R_poloidal,Z_poloidal=np.meshgrid(R_poloidal,Z_poloidal) 
        R_poloidal,Z_poloidal=R_poloidal.flatten(),Z_poloidal.flatten()
        phi_poloidal=np.full(len(R_poloidal),phi)

        dB_toroidal=self.evaluate(R=R_toroidal,phi=phi_toroidal,Z=Z_toroidal,phase=phase,i3dr=i3dr,mode_number=self.mode_number) #evaluate toroidally
        dB_poloidal=self.evaluate(R=R_poloidal,phi=phi_poloidal,Z=Z_poloidal,phase=phase,i3dr=i3dr,mode_number=self.mode_number) #evaluate poloidally

        if absolute:
            dB_toroidal=np.abs(dB_toroidal)
            dB_poloidal=np.abs(dB_poloidal)

        for ax,component_toroidal,component_poloidal,component_name,colour in zip([ax1,ax2,ax3],dB_toroidal,dB_poloidal,['dB_field_R','dB_field_tor','dB_field_Z'],['k','b','r']):
    
            component_poloidal=component_poloidal.reshape(Z_poloidal_dim,R_poloidal_dim) #plot the poloidal variation
            Z_poloidal,R_poloidal=Z_poloidal.reshape(Z_poloidal_dim,R_poloidal_dim),R_poloidal.reshape(Z_poloidal_dim,R_poloidal_dim)
            if not vminmax:
                vminmax=[np.amin(component_poloidal),np.amax(component_poloidal)]
            mesh=ax.pcolormesh(R_poloidal,Z_poloidal,component_poloidal,cmap=colmap,vmin=np.amin(vminmax),vmax=np.amax(vminmax))
            #fig.colorbar(mesh,ax=ax,orientation='vertical')
            
            if absolute: #if plotting absolute value, add tag to axis labels
                component_name='abs( ' + component_name + ' )'

            ax.set_title(component_name)
            ax.set_aspect('equal')

            ax4.plot(phi_toroidal,component_toroidal,colour) #plot the toroidal variation

            if LCFS:
                ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],plot_style_LCFS) #add a LCFS

            if limiters:
                ax.plot(limiters['lcfs_r'],limiters['lcfs_z'],plot_style_limiters) #add a LCFS

        legend=['dB_field_R','dB_field_tor','dB_field_Z']
        if absolute:
            legend=['abs( ' + leg + ' )' for leg in legend]
        ax4.legend(legend)
        ax4.set_title(self.ID)

        if ax_flag is False and fig_flag is False:
            plt.show()

#################################
 
##################################################################
 
###################################################################################################