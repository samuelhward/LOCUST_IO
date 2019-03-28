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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
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
        input_data['B_field_R_real']=[]
        input_data['B_field_R_imag']=[]
        input_data['B_field_Z_real']=[]
        input_data['B_field_Z_imag']=[]
        input_data['B_field_tor_real']=[]
        input_data['B_field_tor_imag']=[]

        #read lazily
        for line in file:
            split_line=line.split()
            input_data['R_2D'].append(float(split_line[0]))
            input_data['Z_2D'].append(float(split_line[1]))
            input_data['B_field_R_real'].append(float(split_line[2]))
            input_data['B_field_R_imag'].append(float(split_line[3]))
            input_data['B_field_Z_real'].append(float(split_line[4]))
            input_data['B_field_Z_imag'].append(float(split_line[5]))
            input_data['B_field_tor_real'].append(float(split_line[6]))
            input_data['B_field_tor_imag'].append(float(split_line[7]))

        input_data['R_2D']=np.asarray(input_data['R_2D'])
        input_data['Z_2D']=np.asarray(input_data['Z_2D'])
        input_data['B_field_R_real']=np.asarray(input_data['B_field_R_real'])
        input_data['B_field_R_imag']=np.asarray(input_data['B_field_R_imag'])
        input_data['B_field_Z_real']=np.asarray(input_data['B_field_Z_real'])
        input_data['B_field_Z_imag']=np.asarray(input_data['B_field_Z_imag'])
        input_data['B_field_tor_real']=np.asarray(input_data['B_field_tor_real'])
        input_data['B_field_tor_imag']=np.asarray(input_data['B_field_tor_imag'])
  
        #infer the grid dimensions and axes
        if input_data['Z_2D'][0]==input_data['Z_2D'][1]: #Z is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(Z_dim,R_dim).T
            input_data['Z_2D']=input_data['Z_2D'].reshape(Z_dim,R_dim).T
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['B_field_R_real']=input_data['B_field_R_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_R_imag']=input_data['B_field_R_imag'].reshape(Z_dim,R_dim).T
            input_data['B_field_Z_real']=input_data['B_field_Z_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_Z_imag']=input_data['B_field_Z_imag'].reshape(Z_dim,R_dim).T
            input_data['B_field_tor_real']=input_data['B_field_tor_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_tor_imag']=input_data['B_field_tor_imag'].reshape(Z_dim,R_dim).T
        else: #R is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(R_dim,Z_dim)
            input_data['Z_2D']=input_data['Z_2D'].reshape(R_dim,Z_dim)
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['B_field_R_real']=input_data['B_field_R_real'].reshape(R_dim,Z_dim)
            input_data['B_field_R_imag']=input_data['B_field_R_imag'].reshape(R_dim,Z_dim)
            input_data['B_field_Z_real']=input_data['B_field_Z_real'].reshape(R_dim,Z_dim)
            input_data['B_field_Z_imag']=input_data['B_field_Z_imag'].reshape(R_dim,Z_dim)
            input_data['B_field_tor_real']=input_data['B_field_tor_real'].reshape(R_dim,Z_dim)
            input_data['B_field_tor_imag']=input_data['B_field_tor_imag'].reshape(R_dim,Z_dim)

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

        for quantity in ['R','phi','Z','time','B_field_R_mag','B_field_tor_mag','B_field_Z_mag','dB_field_R','dB_field_tor','dB_field_Z','divB']:
            input_data[quantity]=[]

        for line in lines:
            for counter,quantity in enumerate(['R','phi','Z','time','B_field_R_mag','B_field_tor_mag','B_field_Z_mag','dB_field_R','dB_field_tor','dB_field_Z','divB']):
            
                try:
                    input_data[quantity].append(float(line.split()[counter]))
                except:
                    input_data[quantity].append(0.0)

        for quantity in ['R','phi','Z','time','B_field_R_mag','B_field_tor_mag','B_field_Z_mag','dB_field_R','dB_field_tor','dB_field_Z','divB']:
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
        locust_names['R','phi','Z','B_field_R_mag','B_field_tor_mag','B_field_Z_mag'] #corresponding LOCUST_IO fields that we want to retain
        for ascot_name_1,ascot_name_2,locust_io_name in zip(ascot_names_1,ascot_names_2,locust_io_names):
            if ascot_name_1 in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name_1])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])
            elif ascot_name_2 in raw_data.keys():
                input_data[locust_io_name]=copy.deepcopy(raw_data[ascot_name_2])
                input_data[locust_io_name]=np.asarray(input_data[locust_io_name])                
    
    return input_data

def read_perturbation_IDS(shot,run,**properties):
    """
    notes:

    """ 

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_perturbation_IDS could not import IMAS module!\nreturning\n")
        return

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
        input_data['B_field_R_real']=[]
        input_data['B_field_R_imag']=[]
        input_data['B_field_Z_real']=[]
        input_data['B_field_Z_imag']=[]
        input_data['B_field_tor_real']=[]
        input_data['B_field_tor_imag']=[]

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
            input_data['B_field_R_real'].append(float(split_line[2]))
            input_data['B_field_R_imag'].append(float(split_line[3]))
            input_data['B_field_Z_real'].append(float(split_line[4]))
            input_data['B_field_Z_imag'].append(float(split_line[5]))
            input_data['B_field_tor_real'].append(float(split_line[6]))
            input_data['B_field_tor_imag'].append(float(split_line[7]))

        input_data['R_2D']=np.asarray(input_data['R_2D'])
        input_data['Z_2D']=np.asarray(input_data['Z_2D'])
        input_data['B_field_R_real']=np.asarray(input_data['B_field_R_real'])
        input_data['B_field_R_imag']=np.asarray(input_data['B_field_R_imag'])
        input_data['B_field_Z_real']=np.asarray(input_data['B_field_Z_real'])
        input_data['B_field_Z_imag']=np.asarray(input_data['B_field_Z_imag'])
        input_data['B_field_tor_real']=np.asarray(input_data['B_field_tor_real'])
        input_data['B_field_tor_imag']=np.asarray(input_data['B_field_tor_imag'])
    
        #infer the grid dimensions and axes
        if input_data['Z_2D'][0]==input_data['Z_2D'][1]: #Z is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(Z_dim,R_dim).T
            input_data['Z_2D']=input_data['Z_2D'].reshape(Z_dim,R_dim).T
            input_data['R_1D']=input_data['R_2D'][:,0]
            input_data['Z_1D']=input_data['Z_2D'][0,:]
            input_data['B_field_R_real']=input_data['B_field_R_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_R_imag']=input_data['B_field_R_imag'].reshape(Z_dim,R_dim).T
            input_data['B_field_Z_real']=input_data['B_field_Z_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_Z_imag']=input_data['B_field_Z_imag'].reshape(Z_dim,R_dim).T
            input_data['B_field_tor_real']=input_data['B_field_tor_real'].reshape(Z_dim,R_dim).T
            input_data['B_field_tor_imag']=input_data['B_field_tor_imag'].reshape(Z_dim,R_dim).T
        else: #R is slowly-varying
            R_dim=int(np.where(input_data['Z_2D']==input_data['Z_2D'][0])[0].size)
            Z_dim=int(np.where(input_data['R_2D']==input_data['R_2D'][0])[0].size)
            input_data['R_2D']=input_data['R_2D'].reshape(R_dim,Z_dim)
            input_data['Z_2D']=input_data['Z_2D'].reshape(R_dim,Z_dim)
            input_data['R_1D']=input_data['R_2D'][0,:]
            input_data['Z_1D']=input_data['Z_2D'][:,0]
            input_data['B_field_R_real']=input_data['B_field_R_real'].reshape(R_dim,Z_dim)
            input_data['B_field_R_imag']=input_data['B_field_R_imag'].reshape(R_dim,Z_dim)
            input_data['B_field_Z_real']=input_data['B_field_Z_real'].reshape(R_dim,Z_dim)
            input_data['B_field_Z_imag']=input_data['B_field_Z_imag'].reshape(R_dim,Z_dim)
            input_data['B_field_tor_real']=input_data['B_field_tor_real'].reshape(R_dim,Z_dim)
            input_data['B_field_tor_imag']=input_data['B_field_tor_imag'].reshape(R_dim,Z_dim)
        
    print("finished reading MARSF perturbation")
    
    return input_data
    
################################################################## Perturbation write functions
 
def dump_perturbation_LOCUST(output_data,filepath,**properties):
    """
    writes perturbation to LOCUST format

    notes:
        dumps data with quickly-varying Z like usual LOCUST input
    """
 
    print("writing LOCUST perturbation")

    with open(filepath,'w') as file: #open file
        
        quantities=['R_2D','Z_2D','B_field_R_real','B_field_R_imag','B_field_Z_real','B_field_Z_imag','B_field_tor_real','B_field_tor_imag']
        np.savetxt(filepath,np.array([output_data[quantity].flatten() for quantity in quantities]).T,fmt='%1.10E',delimiter=' ')

    print("finished writing LOCUST perturbation")

def dump_perturbation_point_data_LOCUST(output_data,filepath='point_data.inp',BCHECK=1,**properties):
    """
    generates the point_data.inp file for checking magnetic perturbations using LOCUST -DBCHECK

    args:
        BCHECK - coordinate format setting for LOCUST field checking (1=RPhiZ,2=XYZ)  
    notes:
        uses R_point_data, phi_point_data, Z_point_data and time_point_data arrays stored in perturbation class as point_data.inp coordinates 
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
    """
 
    LOCUST_input_type='perturbation'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
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
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_perturbation_LOCUST(self.filepath,**properties)

        elif data_format=='LOCUST_field_data': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from LOCUST_field_data - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties={**properties}
                self.data=read_perturbation_LOCUST_field_data(self.filepath,**properties)

        elif data_format=='ASCOT_field_data': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from ASCOT_field_data - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties={**properties}
                self.data=read_perturbation_ASCOT_field_data(self.filepath,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from magnetics IDS - shot and run required\n".format(self.ID),shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_perturbation_IDS(self.shot,self.run,**properties)

        elif data_format=='MARSF': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot read_data() from MARSF - filename required\n".format(self.ID),filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_perturbation_MARSF(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/LOCUST_field_data/ASCOT_field_data/IDS/MARSF)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,BCHECK=1,**properties):
        """
        write perturbation to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)

        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
         
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files+filename
                dump_perturbation_LOCUST(self.data,filepath,**properties)

        elif data_format=='point_data':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: {} cannot dump_data() to point_data.inp - filename required\n".format(self.ID),filename):
                filepath=support.dir_input_files+filename
                dump_perturbation_point_data_LOCUST(self.data,filepath,BCHECK,**properties)

        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/point_data)\n")

    def plot(self,key='B_field_R_real',LCFS=False,limiters=False,number_bins=20,fill=True,colmap=cmap_default,ax=False,fig=False):
        """
        plots a perturbation
        
        notes:
            
        args:
            key - selects which data in perturbation to plot
            LCFS - toggles plasma boundary on/off in 2D plots (requires equilibrium arguement)
            limiters - object which contains limiter data rlim and zlim
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            colmap - set the colour map (use get_cmap names)
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

        #0D data
        if self[key].ndim==0:
            print([key])
            return
        
        #>0D data is plottable
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)
        ax.set_title(self.ID)

        #1D data
        if self[key].ndim==1:
            ax.plot(self[key])
            ax.set_ylabel(key)

        #2D data
        elif self[key].ndim==2:

            X=self['R_1D'] #make a mesh
            Y=self['Z_1D'] 
            Y,X=np.meshgrid(Y,X) #swap since things are defined r,z 
            Z=self[key] #2D array (nR_1D,nZ_1D) of poloidal flux
            
            #2D plot
            if fill is True:
                mesh=ax.contourf(X,Y,Z,levels=np.linspace(np.amin(Z),np.amax(Z),num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
                for c in mesh.collections: #for use in contourf
                    c.set_edgecolor("face")
            else:
                mesh=ax.contour(X,Y,Z,levels=np.linspace(np.amin(Z),np.amax(Z),num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
                if plot_contour_labels:
                    ax.clabel(mesh,inline=1,fontsize=10)
                
            #mesh=ax.pcolormesh(X,Y,Z,colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))

            #3D plot
            #ax=ax.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
            
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')
            ax.set_aspect('equal')
            ax.set_xlim(np.min(self['R_1D']),np.max(self['R_1D']))
            ax.set_ylim(np.min(self['Z_1D']),np.max(self['Z_1D']))
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

            if LCFS:
                ax.plot(self['lcfs_r'],self['lcfs_z'],plot_style_LCFS) 
            if limiters: #add boundaries if desired
                ax.plot(self['rlim'],self['zlim'],plot_style_limiters) 

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()

            
#################################
 
##################################################################
 
###################################################################################################