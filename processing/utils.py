#utils.py
 
"""
Samuel Ward
02/11/2017
----
supporting functions for LOCUST input/output classes
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
    import os
    import re
    import time
    import itertools
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from scipy.io import netcdf as ncdf
    import scipy.interpolate
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import h5py
except:
    print("WARNING: h5py could not be imported!\n") 
try:
    from classes import support
except:
    raise ImportError("ERROR: LOCUST_IO/classes/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from processing import plot_output
except:
    raise ImportError("ERROR: LOCUST_IO/processing/plot_output.py could not be imported!\nreturning\n")


np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
plot_style_LCFS='m-' #set plot style for LCFS
plot_style_limiters='w-' #set plot style for limiters

pi=np.pi
e_charge=1.60217662e-19 #define electron charge
mass_neutron=1.674929e-27 #define mass of neutron
amu=1.66053904e-27
mass_deuterium=2.0141017781*amu



################################################################################################### Supporting functions
 
def none_check(ID,LOCUST_input_type,error_message,*args):
    """
    generic function for checking if a None value appears in *args
 
    notes:
        message should be something specific to the section of code which called none_check
    """
    if all(arg is not None for arg in args):
        return False
    else:
        print("WARNING: none_check returned True (LOCUST_input_type={LOCUST_input_type}, ID={ID}): {message}".format(LOCUST_input_type=LOCUST_input_type,ID=ID,message=error_message))
        return True
         
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
 
def get_next(obj):
    """
    generic object iterator
 
    notes:
        every time get_next() is called, it will advance the counter one
    """
    pyVer = sys.version_info[0]
    if pyVer == 2:
        return obj.next() #NOTE how exactly does .next() work? will it accept tok in toklist values from file_numbers() and allow file_numbers to then read in the next bunch (since yield will pause the While True loop?)
    else:
        return next(obj)
     
def dict_set(data,**kwargs):
    """
    generalised function (based upon safe_set) to set values in dictionary 'data' after checking source data exists 
    
    notes:
    usage:
        dict_set(data,some_key=some_source,some_other_key=[1,2,3,4]) to set multiple values simultaneously
        dict_set(data,**{'some_key':100,'some_other_key':200}) equally
    """
    keys=kwargs.keys()
    values=kwargs.values()
    
    for key,value in zip(keys,values): #loop through kwargs
        if key is not None and value is not None: #use 'is' since we may want to set data to 0
            data[key]=value

def safe_set(target,source):
    """
    generalised function to set a target value to source if it exists
    
    notes:
    """
    if source is not None: #use 'is' since we may want to set data to 0
        target=source

 
def fortran_string(number_out,length,decimals=None,exponential=True):
    """
    produces a string stream of specified length, padded with space at the start

    notes:
        supposed to be a quick fix to fortran format descriptors
    args:
        number_out - variable holding number to write out 
        length - total number of characters for output string stream including white space 
        decimals - required for non-ints, how many decimal places the number is written to 
        exponential - true for output to exponential format, otherwise formatted like float
    """

    if decimals is not None:
        if exponential is True:
            output_string='{:.{decimals}e}'.format(number_out,decimals=decimals)        
        else:
            output_string='{:.{decimals}f}'.format(number_out,decimals=decimals)        
    else: #assume integer if decimals not specified
        output_string=str(number_out)

    number_spaces=length-len(output_string) #caclulate amount of padding needed from length of the number_out
    if number_spaces>=0:
        string_stream='{}{}'.format(' '*number_spaces,output_string) #construct the string stream
        return string_stream
    else:
        print('ERROR: fortran_string() too many decimal places for requested string length!')


def sort_arrays(main_array,*args):
    """
    sort an arbitrary number of arrays in parallel (*args) according to main_array

    notes:
        assumes all arrays are numpy arrays
    usage:
        to sort ascending according to array x, 
        x_sorted,y_sorted,z_sorted,...=sort_arrays(x,y,z,...)
    """

    sorted_indices=main_array.argsort() #get the order of sorted indices
    
    main_array=main_array[sorted_indices] #convert all the arrays
    returned_arrays=[main_array] #add arrays to list
    for array in args:
        array=array[sorted_indices]
        returned_arrays.append(array)

    return returned_arrays

def angle_pol(R_major_width,R,Z):
    """
    returns poloidal angle

    notes:
    """

    R_minor=R-R_major_width
    angle=np.arctan2(Z,R_minor)

    return angle

def interpolate_2D(X_axis,Y_axis,Z_grid,type='RBS',rect_grid=True):
    """
    generate a 2D grid interpolator
    notes:
        keep as separate functions so can freely swap out interpolation method
        RBF - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
            - high memory overhead, most accurate
        RBS - https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html#scipy.interpolate.RectBivariateSpline
              https://scipython.com/book/chapter-8-scipy/examples/two-dimensional-interpolation-with-scipyinterpolaterectbivariatespline/
    args:
        X_axis - 1D x-axis
        Y_axis - 1D y-axis
        Z_grid - 2D z-axis
        rect_grid - toggle whether X_axis, Y_axis are regular rectangular grid edges or arbitrary cordinates 

    usage:
        my_interpolator=interpolate_2D(X_axis,Y_axis,data)
        interpolated_value=my_interpolator(x,y)
    """
    
    if type=='RBF':
        if rect_grid:
            Y_grid,X_grid=np.meshgrid(Y_axis,X_axis) #swap since things are defined r,z 
        else:
            Y_grid,X_grid=Y_axis,X_axis
        interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function='multiquadric',smooth=0)
    
    elif type=='RBS':
        interpolator=scipy.interpolate.RectBivariateSpline(X_axis,Y_axis,Z_grid) #normally order is other way in RBS but I have swapped my axes

    return interpolator

def interpolate_1D(X_axis,Y_axis):
    """
    generate a 1D line interpolator

    notes:
        keep as separate functions so can freely swap out interpolation method
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    """

    interpolator=scipy.interpolate.Rbf(X_axis,Y_axis,function='cubic',smooth=0)

    return interpolator

def RphiZ_to_XYZ(R,phi):
    """
    converts R,phi positions to X,Y 
    """

    X=R*np.sin(phi)
    Y=R*np.cos(phi)

    return X,Y

def dump_profiles_ASCOT(filename,temperature_i,temperature_e,density_i,density_e,rotation_toroidal):
    """
    dumps collection of kinetic profiles to ASCOT input.plasma_1d format

    notes:
        profiles must be mapped to same poloidal flux axis
        currently allows single ion species - in future, pass list of profile objects and iterate
        if rotation profile added to LOCUST data in future, will need to update this function accordingly
    args:
        filename - output filename
        temperature_e - electron temperature object (eV)
        temperature_i - ion temperature object (eV)
        density_e - electron density object (#/m^3)
        density_i - ion density object (#/m^3)
        rotation_toroidal - array-like toroidal rotation (rad/s)
    """

    print("dumping profiles to ASCOT format")
 
    filepath=support.dir_input_files+filename
    with open(filepath,'w') as file:

        file.write("# Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation\n")
        file.write("# range must cover [0,1] of normalised poloidal rho. It can exceed 1. \n")
        file.write("# 18Jan08 for testing (first 3 lines are comment lines)\n")
        file.write(str(temperature_e['flux_pol_norm'].size)+"   "+"1"+" # Nrad,Nion\n")
        file.write("1           # ion Znum\n")
        file.write("1           # ion Amass\n")
        file.write("1 1         # collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons\n")
        file.write("    RHO (pol)       Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)\n")

        flux_pol_norm_sqrt=np.sqrt(np.abs(temperature_e['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,temperature_e,density_e,rotation_toroidal,temperature_i,density_i=sort_arrays(flux_pol_norm_sqrt,temperature_e['T'],density_e['n'],rotation_toroidal,temperature_i['T'],density_i['n']) #check order

        for RHO,Te,Ne,Vtor_I,Ti1,Ni1 in zip(flux_pol_norm_sqrt,temperature_e,density_e,rotation_toroidal,temperature_i,density_i): 
            line=fortran_string(RHO,16,7)+fortran_string(Te,16,7)+fortran_string(Ne,16,7)+fortran_string(Vtor_I,15,7)+fortran_string(Ti1,17,7)+fortran_string(Ni1,15,7)+"\n"
            file.write(line)

    print("finished dumping profiles to ASCOT format")

def dump_wall_ASCOT(filename,output_data):
    """
    dumps 2D wall outline to ASCOT input.wall_2d format
    
    notes:
        currently sets all divertor flags = 0 i.e. treats everything as 'wall'
    args:
        filename - output filename
        output_data - data structure holding wall with same wall variable names as GEQDSK equilibrium i.e. rlim,zlim

    """

    print("dumping wall to ASCOT format")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file:

        file.write("{number_points} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n".format(number_points=int(output_data['rlim'].size)))
        
        for r,z in zip(output_data['rlim'],output_data['zlim']):
            line=fortran_string(r,16,7)+fortran_string(z,16,7)+fortran_string(0.0,4,0,False)+"\n"
            file.write(line)

    print("finished dumping wall to ASCOT format")

def dump_rotation_MARSF(filename,output_data):
    """
    writes rotation profile to MARSF Mogui ASCII format 

    notes
        assumes structure similar number_density or temperature, with normalised poloidal flux given against rotation
        re-usable if rotation class written for LOCUST_IO
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    args:
        filename - output filename
        output_data - data structure holding 'rotation' against normalised poloidal flux
    """

    print("writing rotation to MARSF mogui")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file: #open file

        flux_pol_norm_sqrt=np.sqrt(np.abs(output_data['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,rotation=sort_arrays(flux_pol_norm_sqrt,output_data['rotation']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(flux_pol_norm_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(flux_pol_norm_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm_sqrt}{rotation}\n".format(flux_pol_norm_sqrt=fortran_string(flux_pol_norm_sqrt[point],24,18),rotation=fortran_string(rotation[point],24,18)))

    print("finished writing rotation to MARSF mogui")


def dump_coil_currents_MARSF(filename,output_data):
    """
    writes coil currents to MARSF Mogui ASCII format

    notes:
        output_data structure assumed to be dict holding 'upper_coil_currents' and 'lower_coil_currents' arrays
        MARSF mogui written by David Ryan
    args:  
        filename - output filename
        output_data - data structure holding coil currents
    """

    print("writing coil currents to MARSF mogui")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file: #open file

        some_ID='34835_2000'
        file.write('{}\n'.format(some_ID)) #write headerline
        
        upper_coil_currents=''
        lower_coil_currents=''

        for coil in output_data['upper_coil_currents']:
            upper_coil_currents+=str(coil)+', '

        for coil in output_data['lower_coil_currents']:
            lower_coil_currents+=str(coil)+', '

        upper_coil_currents=upper_coil_currents[0:-2]
        lower_coil_currents=lower_coil_currents[0:-2]

        upper_coil_currents+='\n'
        lower_coil_currents+='\n'

        file.write(upper_coil_currents) #write upper coil currents
        file.write(lower_coil_currents) #write lower coil currents

    print("finished writing coil currents to MARSF mogui")



def ASCOT_run_gen(run_file='ascot4.cmd',initialdir=None,output_file='ascot.out',max_proc=50,min_proc=25,error_file='ascot.err',executable='test_ascot'):
    """
    writes out freia batch file for ASCOT run

    notes:
    """

    if initialdir is None:
        initialdir=os.path.dirname(os.path.abspath(__file__)) #use cwd
    
    with open(run_file,'w') as file:
        file.write('# @ input = dev/null/\n')
        file.write('# @ initialdir = {initialdir}\n'.format(initialdir=initialdir))
        file.write('# @ output = {output_file}\n'.format(output_file=output_file))
        file.write('# @ error = {error_file}\n'.format(error_file=error_file))
        file.write('# @ jobtype = openmpi\n')
        file.write('# @ max_processors = {max_proc}\n'.format(max_proc=max_proc))
        file.write('# @ min_processors = {min_proc}\n'.format(min_proc=min_proc))
        file.write('# @ queue\n\n')
        file.write('date\n')
        file.write('mpirun -np $NSLOTS {executable} -output ascot_freia_$JOB_ID\n'.format(executable=executable))
        file.write('date\n')











################################################################################################### Supporting classes

################################################################## TRANSP

class TRANSP_output:
    """
    base class for generic CDF file produced by TRANSP
    
    notes:
        TRANSP_output and child classes usually have different methodsto LOCUST_IO objects and function differently
        as a result, the classes tend to purely use their own methods to combat confusion  
        please see individual class/method docstrings for help
    """

    def __init__(self,filename):
        """
        constructor for TRANSP_output class
        
        notes: 

        """

        self.filename=filename
        self.data={}
        self.read_data(self.filename)

    def __getitem__(self,key):
        """
        special method to set member data via []

        notes:
            returns the nested netCDF object associated with 'key', not the raw data
        usage:
            access member data via my_ncdf[key].data, attributes via my_ncdf[key]._attributes etc
            
        """
    
        return self.data[key]

    def __setitem__(self,key,value):
        """
        special method to set member data via []
        """

        self.data[key]=value
    
    def read_data(self,filepath):
        """
        read data from TRANSP output netCDF file into LOCUST_IO-like data structure

        notes:
            designed to be overloaded by child classes of TRANSP_output with specific netCDF data structures
        """

        pass

    '''
    def new_var(self,key,typ,**dimensions):
        """
        interface for creating new netCDF variable

        notes:
            if dimension is supplied that does not already exist then new dimension is created
            final variable shape does not currently reflect order of **dimensions arguements
        args:
            key - name of variable to create
            typ - dtype of variable
            dimensions - kwarg for dimension names and corresponding lengths
        useage:
            my_cdf.new_var('new_variable','float',x=100,y=50)
        """

        for dimension in dimensions.keys():
            if dimension not in self.data.dimensions:
                self.data.dimensions[dimension]=dimensions[dimension]
            elif dimensions[dimension]!=self.data.dimensions[dimension]:
                print("ERROR - new_var(): specified dimension {dimension} already exists with value={value}!\n".format(dimension=dimension,value=dimensions[dimension]))
                return 

        if key in self.data.variables.keys():
            print("{filename}: overwriting {key}".format(filename=self.filename,key=key))
        self.data.createVariable(key,typ,dimensions)
    '''

class TRANSP_output_FI(TRANSP_output):
    """
    class to describe RunNumberCode_FI.cdf fast ion output netCDF files generated by get_fbm 

    notes:
        because distribution function collected on non-LOCUST grid must specify custom class methods here
        
    """

    def read_data(self,filename):
        """
        read data from TRANSP output netCDF file into LOCUST_IO-like data structure
        
        notes:
            reads all data into nested dictionary before renaming and reshaping some data according to LOCUST_IO variable names  
            my_TRANSP_output_FI[key]['data'] returns data from netCDF
            my_TRANSP_output_FI[key]['units'] returns any associated units if any 
            my_TRANSP_output_FI[key]['description'] returns any associated description if any
        """

        print("reading TRANSP fast ion distribution function")

        self.filename=filename
        self.filepath=support.dir_output_files+filename
        file=ncdf.netcdf_file(self.filepath,'r') #open netCDF file

        for key in file.variables.keys(): #read everything into a nested dictionary
            self[key]={}
            try: #this try except structure is more suitable than using hasattr
                self[key]['units']=copy.deepcopy(file.variables[key].units.decode('utf-8'))
            except:
                self[key]['units']=None
            try:
                self[key]['description']=copy.deepcopy(file.variables[key].long_name.decode('utf-8'))
            except:
                self[key]['description']=None
            try:   
                self[key]['data']=copy.deepcopy(file.variables[key].data)
            except:
                self[key]['data']=None

        del(file) #close file

        #map some data to more LOCUST_IO-style variable names here
        self['F_D_NBI']['data']*=.5e6 #get dfn in terms of dPitch^-1 m^-3 
        self['F_D_NBI']['units']='#/m^3/eV/dPitchBin'
        self['R2D']['data']*=.01
        self['R2D']['units']='[metres]'
        self['Z2D']['data']*=.01
        self['Z2D']['units']='[metres]'

        print("finished reading TRANSP fast ion distribution function")


    def dfn_integrate(self,space=True,pitch=True,energy=True):
        """
        integrate the fast ion distribution function over specified dimensions

        notes:
            assumes regular pitch and energy grids
            dimensions=[space,pitch,energy]
        args:
            space - toggle to integrate over space
            pitch - toggle to integrate over pitch
            energy - toggle to integrate over energy 
        """

        print("integrating TRANSP fast ion distribution function")

        self['F_D_NBI_int']={}
        self['F_D_NBI_int']['data']=copy.deepcopy(self['F_D_NBI']['data'])

        sum_indices=[] #figure out which indices we need to sum over
        if energy: #apply Jacobian and integrate each dimension
            dE=np.abs(self['E_D_NBI']['data'][1]-self['E_D_NBI']['data'][0]) #energy bin width
            self['F_D_NBI_int']['data']*=dE
            sum_indices.append(2)
        if pitch:
            dP=np.abs(self['A_D_NBI']['data'][1]-self['A_D_NBI']['data'][0]) #pitch bin width
            self['F_D_NBI_int']['data']*=dP
            sum_indices.append(1)
        if space: #integrate over various dimensions here
            for counter,volume_element in enumerate(self['BMVOL']['data']):
                self['F_D_NBI_int']['data'][counter,:,:]*=volume_element
            sum_indices.append(0)

        for sum_index in sum_indices: #sum over unwanted dimensions - must be descending order to avoid changing desired index to sum over
            self['F_D_NBI_int']['data']=np.sum(self['F_D_NBI_int']['data'],axis=sum_index)

        print("finished integrating TRANSP fast ion distribution function")

    def dfn_plot(self,some_equilibrium=None,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,ax=False,fig=False,**kwargs):
        """
        plot the distribution function

        notes:
            functions differently to LOCUST_IO's plot_distribution_function() (see args list below)
            assumes distribution function has already been integrated to specified 'axes' option
        args:
            some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
            axes - array specifying plot type from list of options below
            LCFS - show plasma boundary outline (requires equilibrium arguement)
            limiters - toggles limiters on/off in 2D plot
            real_scale - plot to Tokamak scale
            colmap - select desired colourmap
            ax - external axis object
            fig - external figure object
        axes options:
            R,Z - integrate over pitch and velocity [m]^-3
            E,V_pitch - integrate over space and transform to [eV]^-1[dpitch]^-1  
            E,time - [eV]^-1 over multiple timesteps (supply list of additional objects in **kwargs e.g. ...fig=False, TRANSP_output_FI_list=[FI_CDF_Time1,FI_CDF_Time2])
            E - [eV]^-1
            R - [m]^-3 
            N - total #
        """
        
        if not ax:
            ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
        else:
            ax_flag=True
        if not fig:
            fig_flag=False
        else:
            fig_flag=True
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)   

        #add specific options for plotting here
        
        if axes==['R','Z']:
            #R,Z=np.meshgrid(self['R2D']['data'],self['Z2D']['data'])

            X=np.linspace(np.min(self['R2D']['data']),np.max(self['R2D']['data']),np.sqrt(len(self['R2D']['data']))) #need to interpolate since irregular grid
            Y=np.linspace(np.min(self['Z2D']['data']),np.max(self['Z2D']['data']),np.sqrt(len(self['Z2D']['data'])))
            X,Y=np.meshgrid(X,Y)
            interpolator=interpolate_2D(self['Z2D']['data'],self['R2D']['data'],self['F_D_NBI_int']['data'],type='RBF',rect_grid=False)
            new_dfn=interpolator(Y,X)
            ax.set_facecolor(colmap(np.amin(new_dfn)))
            mesh=ax.pcolormesh(X,Y,new_dfn,cmap=colmap,vmin=np.amin(new_dfn),vmax=np.amax(new_dfn))
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

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
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')

        elif axes==['E','V_pitch']: 
            E,V_pitch=np.meshgrid(self['E_D_NBI']['data'],self['A_D_NBI']['data']) #X,Y this way because F_D_NBI dimension ordering
            ax.set_facecolor(colmap(np.amin(self['F_D_NBI_int']['data'])))
            mesh=ax.pcolormesh(E,V_pitch,self['F_D_NBI_int']['data'],cmap=colmap,vmin=np.amin(self['F_D_NBI_int']['data']),vmax=np.amax(self['F_D_NBI_int']['data']))            
            ax.set_xlabel('energy [eV]')
            ax.set_ylabel('pitch (v_parallel/v)')
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')            

        elif axes==['E','time']: #use the kwargs to accept TRANSP_output_FI_list, a list of other integrated TRANSP_output_FI objects to plot in a single mesh plot
            if 'TRANSP_output_FI_list' in kwargs.keys():

                E_list=[self['F_D_NBI_int']['data']]
                time_list=[self['TIME']['data']]
                
                for FI_object in kwargs['TRANSP_output_FI_list']:
                    E_list.append(FI_object['F_D_NBI_int']['data'])
                    time_list.append(FI_object['TIME']['data'])

                F_D_NBI_int_all=np.array(E_list,ndmin=2).T #combine all objects in list into one
                time,E=np.meshgrid(np.array(time_list),self['E_D_NBI']['data']) #automatically sorts in ascending time since pcolormesh does not require increasing meshgrid
                mesh=ax.pcolormesh(E,time,F_D_NBI_int_all,cmap=cmap_default,vmin=np.amin(F_D_NBI_int_all),vmax=np.amax(F_D_NBI_int_all))
                ax.set_xlabel('energy [eV]')                
                ax.set_ylabel('time [s]')
                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')
            
            else:
                print("ERROR: 'TRANSP_output_FI_list' not supplied in TRANSP_output_FI.dfn_plot(...**kwargs) - please see docstring for TRANSP_output_FI.dfn_plot()\n")

        #elif axes==['E']:

        #elif axes==['R']:

        else:
            print("ERROR: TRANSP_output_FI.dfn_plot() given unknown axes option - please check dfn_plot() docstring\n")

        if ax_flag is False and fig_flag is False:
            plt.show()

################################################################## ASCOT

class ASCOT_output:
    """
    class for encapsulating ASCOT HDF5 output file

    notes:
        very hacky
        mimics a LOCUST_IO object - use pull_data and methods like dfn_transform to then access standard LOCUST_IO functions
        my_output.file['key/path/to/data'].values will return leaf-level tree data from HDF5 file
    example:
        my_output=ASCOT_output('ID',some_filename)
        my_output.pull_data('distribution_function')
        my_output.dfn_plot(axes=['R','Z'],real_scale=True) 
    """

    def __init__(self,ID,filename=None):
        """
        constructor for TRANSP_output class
        
        notes: 

        """

        self.ID=ID
        self.data={}
        if filename: #if supplied new filename, overwrite previous filepath
            self.filename=filename
            self.filepath=support.dir_output_files+filename

    def __getitem__(self,key):
        """
        special method to get member data via []

        notes:
            this does not return the data at leaf level, this must be done using .values attribute
            could adapt this to neatly print the whole tree with '\t * #recursion levels'
        usage:
        """
    
        return self.data[key]

    def __setitem__(self,key,value):
        """
        special method to set member data via []
        """

        self.data[key]=value

    def look(self,key=None):
        """
        prints file sub-branches from branch 'key'

        notes:
            can navigate tree with '/' like file directory
        usage:
            my_ASCOT_output.look() #print top level
            my_ASCOT_output.look('bfield')
        """

        if not self.file:
            print("no file member data found by .look()! use .file_open()")
            return

        if key:
            if 'keys' in dir(self.file[key]):
                for var in self.file[key].keys():
                    print(var)
            else:
                print(self.file[key].value)
        else:
            for var in self.file.keys():
                print(var)

    def file_open(self,filename=None):
        """
        re-open/open new HDF5 file
        """

        if filename: #if supplied new filename, overwrite previous filepath
            self.filename=filename
            self.filepath=support.dir_output_files+filename
        self.file=h5py.File(self.filepath,'r') 

    def file_close(self):
        """
        close annoying HDF5 file handle

        notes:
        """

        del(self.file)

    def pull_data(self,datatype='distribution_function'):
        """
        deep copy and pull data from HDF5 file 

        notes:
            extracts and maps out the output HDF5 file to LOCUST_IO variable names where possible
            see individual datatype branches for more detail
            NOTE work in progress

        args:
            datatype - read data from ASCOT output file and create structure resembling this LOCUST_IO datatype    
        """

        self.datatype=datatype
        self.file_open()

        if datatype=='distribution_function':
            ''' 
            extract dfn to #[m^-3 eV^-1 dpitch^-1] format

            file['distributions/rzPitchEdist/ordinate'].value=[species?,time?,E,V_pitch,Z,R,ordinate]
            file['distributions/rzPitchEdist/abscissae']['dim1']=R
                                                     ...['dim2']=Z
                                                     ...['dim3']=V_pitch
                                                     ...['dim4']=E [J]
                                                     ...['dim5']=time?
                                                     ...['dim6']=species?'''
            
            self['R']=.5*(self.file['distributions/rzPitchEdist/abscissae/dim1'].value[1:]+self.file['distributions/rzPitchEdist/abscissae/dim1'].value[:-1]) #bin centres [m]
            self['Z']=.5*(self.file['distributions/rzPitchEdist/abscissae/dim2'].value[1:]+self.file['distributions/rzPitchEdist/abscissae/dim2'].value[:-1]) #[m]
            self['V_pitch']=.5*(self.file['distributions/rzPitchEdist/abscissae/dim3'].value[1:]+self.file['distributions/rzPitchEdist/abscissae/dim3'].value[:-1])
            self['E']=.5*(self.file['distributions/rzPitchEdist/abscissae/dim4'].value[1:]+self.file['distributions/rzPitchEdist/abscissae/dim4'].value[:-1])/e_charge #[eV]
            self['V']=np.sqrt(2.*self['E']/mass_neutron)

            self['dR']=np.abs(self['R'][1]-self['R'][0])
            self['dZ']=np.abs(self['Z'][1]-self['Z'][0])
            self['dV_pitch']=np.abs(self['V_pitch'][1]-self['V_pitch'][0])
            self['dE']=np.abs(self['E'][1]-self['E'][0])
            
            self['nR']=np.array(len(self['R']))
            self['nZ']=np.array(len(self['Z']))
            self['nV_pitch']=np.array(len(self['V_pitch']))
            self['nE']=np.array(len(self['E']))

            self['dfn']=self.file['distributions/rzPitchEdist/ordinate'].value #[m^-3 J^-1]
            self['dfn']=np.sum(self['dfn'],axis=0)
            self['dfn']=np.sum(self['dfn'],axis=0)
            self['dfn']=np.sum(self['dfn'],axis=-1)
            self['dfn']=np.swapaxes(self['dfn'],-1,-2) 
            self['dfn']*=e_charge/self['dV_pitch'] #[m^-3 eV^-1 dpitch^-1] XXX check this
            
            self['dfn_index']=np.array(['dummy','E','V_pitch','R','Z'])

        elif datatype=='equilibrium':
            pass #unfinished
        



        self.file_close()

    def dfn_transform(self,axes=['R','Z']):
        """
        transforms and integrates the distribution function according to pre-defined configurations 
        
        args:
            axes - the dimensions over which to transform the DFN to
        notes:
            'overloads' processing.process_output.dfn_transform
            overwrites dfn generated by pull_data() - call pull_data() again to reset
            remember dimensions of unedited dfn are [E,V_pitch,R,Z] #[m^-3 eV^-1 dpitch^-1]
            assumes unedited dfn
            assumes the bin widths for a given dimension are constant
            assumes toroidal symmetry (no toroidal dimension in dfn)
            assumes user has execute self.pull_data('distribution_function')
            if an array of indices is given, then slice the dfn accordingly and return without any integration
                note for an infinite slice, axes will need to contain slice() objects e.g. axes=[0,0,slice(None),slice(None)] for all R,Z values

        axes options:
            R,Z - integrate over pitch and energy [m]^-3
            E,V_pitch - integrate over space and transform to [eV]^-1[dpitch]^-1 
            E - [eV]^-1
            R - [m]^-3
            N - total #
        """

        #begin list of specific options

        if axes==['R','Z']:
            self.data['dfn']*=self.data['dE']*self.data['dV_pitch'] #integrate
            for counter in range(2): #sum
                self.data['dfn']=np.sum(self.data['dfn'],axis=0) #XXX IS A FACTOR OF 2 NEEDED WHEN INTEGRATING OVER PITCH LIKE TRANSP?

        elif axes==['E','V_pitch']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(self.data['R'])):
                self.data['dfn'][:,:,r,:]*=self.data['R'][r]*2.*pi*self.data['dR']*self.data['dZ']
            #then need to integrate over the unwanted coordinates
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #over Z
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #over R

        elif axes==['E']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(self.data['R'])):
                self.data['dfn'][:,:,r,:]*=self.data['R'][r]*2.*pi*self.data['dR']*self.data['dZ']
            self.data['dfn']*=self.data['dV_pitch'] #integrate over pitch
            
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #sum over Z
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #sum over R
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #sum over V_pitch

        elif axes==['R']:
            self.data['dfn']*=self.data['dE']*self.data['dV_pitch'] #integrate
            for counter in range(2): #sum over pitch and energy
                self.data['dfn']=np.sum(self.data['dfn'],axis=0)
            self.data['dfn']=np.sum(self.data['dfn'],axis=-1) #sum over Z

        elif axes==['N']:
            #applying full Jacobian and integrate over toroidal angle
            for r in range(len(self.data['R'])):
                self.data['dfn'][:,:,r,:]*=self.data['R'][r]*2.*pi*self.data['dR']*self.data['dZ']*self.data['dV_pitch']*self.data['dE']
            for all_axes in range(self.data['dfn'].ndim): #sum over all dimensions
                self.data['dfn']=np.sum(self.data['dfn'],axis=0) 
        

        #general option
        
        elif len(axes)==self.data['dfn'].ndim: #if user supplies all axes then slice WITHOUT integrating
            self.data['dfn']=self.data['dfn'][tuple(axes)]
        else:
            print("ERROR: dfn_transform given invalid axes arguement: "+str(axes))

    def dfn_plot(self,some_equilibrium=None,key='dfn',axes=['R','Z'],LCFS=False,real_scale=False,colmap=cmap_default,transform=True,ax=False,fig=False):
        """
        wrapper to plot_distribution_function

        notes:
        """

        if transform: 
            self.dfn_transform(axes=axes)
        plot_output.plot_distribution_function(self,some_equilibrium,key,axes,LCFS,real_scale,colmap,False,ax,fig) #call standard plot_distribution function but with LOCUST_IO version of transform disabled
        self.pull_data(datatype='distribution_function') #re-pull data that has been overwritten by dfn_transform



#################################
 
##################################################################
 
###################################################################################################