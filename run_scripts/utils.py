#processing.utils.py
 
"""
Samuel Ward
24/08/2018
----
supporting functions for LOCUST run scripts
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
    import os
    import copy
    import subprocess
    import numpy as np
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning")
    sys.exit(1)

try:
    import processing.utils 
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
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


###################################################################################################
#Main Code


################################################################################################### TRANSP classes and functions

def TRANSP_get_fbm_FI_CDF(run_ID,number_files,particle_position=True,guiding_centre=True,device='d3d'):
    """
    notes:
        looks for files in output_files
    args:
        run_ID - full transp run number e.g. 157418S01
        number_files - total number of .DATA# files to extract CDF from 
        particle_position - toggle whether to generate set of CDFs at particle positions
        guiding_centre - toggle whether to generate set of CDFs at guiding centres
        device - device code for machine under study
    """

    project_dir=os.getcwd()
    os.chdir(support.dir_output_files) #change working directory to output files briefly due to bug in CCFE get_fbm installation
    
    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(run_ID,'_fi_',file_ID,'.cdf') #current output file
        output_filename_gc='{}{}{}{}'.format(run_ID,'_fi_',file_ID,'_gc.cdf')

        if guiding_centre:
            fbm_input = """
                        {run_ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        c
                        """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device)

            print("writing TRANSP FI netCDF file {}".format(output_filename_gc))
            proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out, err = proc.communicate(input=fbm_input.encode('utf-8'))
            try:
                os.rename(output_filename,output_filename_gc) #add '_gc' label onto file
            except:
                print("WARNING: TRANSP_get_fbm_FI_CDF could not rename {} to {}".format(output_filename,output_filename_gc))
                    
            print("finished writing TRANSP FI netCDF file {}".format(output_filename_gc))

        if particle_position:
            fbm_input = """
                        {run_ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        p
                        """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device)

            print("writing TRANSP FI netCDF file {}".format(output_filename))
            proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out, err = proc.communicate(input=fbm_input.encode('utf-8'))
            print("finished writing TRANSP FI netCDF file {}".format(output_filename))

    os.chdir(project_dir) #change back to original working directory    

def TRANSP_get_fbm_FI_birth_deposition(run_ID,number_files,device='d3d'):
    """
    notes:
        looks for files in output_files
        assumes:
            all beams
            all energy components
            dumping at particle location
            dumping in x y z space
            uses Akima Hermite spline interpolation
            random seed = 1
            sample size = 1,000,000
    args:
        run_ID - full transp run number e.g. 157418S01
        number_files - total number of .DATA# files to extract CDF from 
        device - device code for machine under study
    """

    project_dir=os.getcwd()
    os.chdir(support.dir_output_files) #change working directory to output files briefly due to bug in CCFE get_fbm installation
    
    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(run_ID,'_fdep_nb_en_',file_ID,'.out') #current output file

        fbm_input = """
                    {run_ID}
                    {path}
                    {file_ID}
                    t
                    {device}
                    1
                    d
                    n
                    0
                    0
                    p
                    b
                    c
                    Y
                    1
                    1000000
                    {output_filename}
                    """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device,output_filename=output_filename)

        print("writing TRANSP FI random sample birth deposition {}".format(output_filename))
        proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        out, err = proc.communicate(input=fbm_input.encode('utf-8'))
        print("finished writing TRANSP FI random sample birth deposition {}".format(output_filename))

    os.chdir(project_dir) #change back to original working directory  


class TRANSP_output:
    """
    base class for generic CDF file produced by TRANSP
    
    notes:
        TRANSP_output and child classes usually have different methodsto LOCUST_IO objects and function differently
        as a result, the classes tend to purely use their own methods to combat confusion  
        please see individual class/method docstrings for help
    """

    def __init__(self,ID,filename):
        """
        constructor for TRANSP_output class
        
        notes: 

        """

        self.ID=ID
        self.filename=filename
        self.filepath=support.dir_output_files+filename
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

    def look(self):
        """
        print class information and data
        """

        print("\n-----------------------")
        print("ID - {ID}".format(ID=self.ID))  

        if hasattr(self,'filename'):
            print("Input Filename - {filename}".format(filename=self.filename))

        if hasattr(self,'properties') and self.properties:
            print("Properties:".format(properties=self.properties))
            for key in self.properties:
                if self.properties[key]: #do not print if the data is empty
                    print("{key} - {value}".format(key=key,value=self.properties[key])) 
        
        print("|")
        print("|")

        if hasattr(self,'data') and self.data:
            for key in self.data:
                if type(self.data[key])==type({}): #check for dicts since these mess things up
                    if self.data[key]: #if dict then check if dict is empty
                        print(key+":")
                        for sub_key in self.data[key]: 
                            print("     {sub_key} - {value}".format(sub_key=sub_key,value=self.data[key][sub_key]))
                elif not self.data[key].size==0: #if not a dict assume a numpy array
                    print("{key} - {value}".format(key=key,value=self.data[key]))
                
        print("-----------------------\n")
    
    def read_data(self,filename):
        """
        read data from TRANSP output netCDF file into LOCUST_IO-like data structure

        notes:
            designed to be overloaded by child classes of TRANSP_output with specific netCDF data structures
        """

        try:
            from scipy.io import netcdf as ncdf 
        except:
            raise ImportError("ERROR: TRANSP_output.read_data could not import netcdf module!\nreturning\n")
            return

        self.filename=filename
        self.filepath=support.dir_output_files+filename
        file=ncdf.netcdf_file(self.filepath,'r') #open netCDF file

        for key in file.variables.keys(): #read everything into a nested dictionary
            try: #this try except structure is more suitable than using hasattr
                self[key]=copy.deepcopy(file.variables[key].data)
            except:
                self[key]=None

        del(file) #close file

class TRANSP_output_FI(TRANSP_output):
    """
    class to describe RunNumberCode_FI.cdf fast ion output netCDF files generated by the get_fbm program 

    notes:
        because distribution function collected on non-LOCUST grid must specify custom class methods here
        
    """

    def read_data(self,filename):
        """
        read data from TRANSP output netCDF file into LOCUST_IO-like data structure
        
        notes:
            reads all data into nested dictionary before renaming and reshaping some data according to LOCUST_IO variable names  
            my_TRANSP_output_FI[key] returns data from netCDF
        """

        print("reading TRANSP fast ion distribution function")

        try:
            from scipy.io import netcdf as ncdf 
        except:
            raise ImportError("ERROR: TRANSP_output_FI.read_data could not import netcdf module!\nreturning\n")
            return

        self.filename=filename
        self.filepath=support.dir_output_files+filename
        file=ncdf.netcdf_file(self.filepath,'r') #open netCDF file

        for key in file.variables.keys(): #read everything into a nested dictionary
            try: #this try except structure is more suitable than using hasattr
                self[key]=copy.deepcopy(file.variables[key].data)
            except:
                self[key]=None

        del(file) #close file

        #map some data to more LOCUST_IO-style variable names here
        self['F_D_NBI']*=.5e6 #get dfn in terms of dPitch^-1 m^-3 eV^-1
        self['R2D']*=.01 #convert to metres
        self['Z2D']*=.01 #convert to metres
        self['BMVOL']*=1.e-6 #convert to m^3
        self['dVOL']=self.data.pop('BMVOL')
        self['dE']=np.abs(self['E_D_NBI'][1]-self['E_D_NBI'][0]) #energy bin width
        self['dV_pitch']=np.abs(self['A_D_NBI'][1]-self['A_D_NBI'][0])
        self['dfn']=self.data.pop('F_D_NBI')
        self['E']=self.data.pop('E_D_NBI')
        self['V_pitch']=self.data.pop('A_D_NBI')
        self['dfn_index']=['RZ','V_pitch','E']

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

        dfn_copy=copy.deepcopy(self)

        sum_indices=[] #figure out which indices we need to sum over
        if energy: #apply Jacobian and integrate each dimension
            dfn_copy['dfn']*=dfn_copy['dE']
            sum_indices.append(2)
        if pitch:
            dfn_copy['dfn']*=dfn_copy['dV_pitch']
            sum_indices.append(1)
        if space: #integrate over various dimensions here
            for counter,volume_element in enumerate(dfn_copy['dVOL']):
                dfn_copy['dfn'][counter,:,:]*=volume_element
            sum_indices.append(0)

        sum_indices.sort(reverse=True) #must be descending order to avoid changing desired index to sum over
        for sum_index in sum_indices: #sum over unwanted dimensions
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=sum_index)

        print("finished integrating TRANSP fast ion distribution function")

        return dfn_copy

    def dfn_plot(self,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,number_bins=20,fill=True,vminmax=None,ax=False,fig=False,**kwargs):
        """
        plot the distribution function

        notes:
            functions differently to LOCUST_IO's plot_distribution_function() (see args list below)
            assumes distribution function has already been integrated to specified 'axes' option
        args:
            axes - array specifying plot type from list of options below
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - plot to Tokamak scale
            colmap - select desired colourmap
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            vminmax - set mesh Vmin/Vmax values
            ax - external axis object
            fig - external figure object
        axes options:
            R,Z - assumes integrated over pitch and velocity [m]^-3
            E,V_pitch - assumes integrated over space and transform to [eV]^-1[dpitch]^-1  
            E,time - [eV]^-1 over multiple timesteps (supply list of additional objects in **kwargs e.g. ...fig=False, TRANSP_output_FI_list=[FI_CDF_Time1,FI_CDF_Time2])
            E - [eV]^-1
            R - [m]^-3 x
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

        ax.set_title(self.ID) #set title to object's ID descriptor

        
        #add specific options for plotting here
        if axes==['R','Z']:

            dfn_copy=self.dfn_integrate(space=False)

            R=np.linspace(np.min(dfn_copy['R2D']),np.max(dfn_copy['R2D']),np.sqrt(len(dfn_copy['R2D']))) #need to interpolate since irregular grid
            Z=np.linspace(np.min(dfn_copy['Z2D']),np.max(dfn_copy['Z2D']),np.sqrt(len(dfn_copy['Z2D'])))
            R,Z=np.meshgrid(R,Z)
            interpolator=processing.utils.interpolate_2D(dfn_copy['Z2D'],dfn_copy['R2D'],dfn_copy['dfn'],type='RBF',rect_grid=False)
            new_dfn=interpolator(Z,R)

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.amin(new_dfn['dfn'])
                vmax=np.amax(new_dfn['dfn'])

            if fill:
                ax.set_facecolor(colmap(np.amin(new_dfn)))
                mesh=ax.pcolormesh(R,Z,new_dfn,cmap=colmap,vmin=vmin,vmax=vmax)
            else:
                mesh=ax.contour(R,Z,new_dfn,levels=np.linspace(vmin,vmax,num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                #ax.clabel(mesh,inline=1,fontsize=10)

            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

            if real_scale is True: #set x and y plot limi`ts to real scales
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            if LCFS: #plot plasma boundary
                ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],plot_style_LCFS) 
            if limiters: #add boundaries if desired
                ax.plot(limiters['rlim'],limiters['zlim'],plot_style_limiters)
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh

        elif axes==['E','V_pitch']: 

            dfn_copy=self.dfn_integrate(pitch=False,energy=False)

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.amin(dfn_copy['dfn'])
                vmax=np.amax(dfn_copy['dfn'])

            E,V_pitch=np.meshgrid(dfn_copy['E'],dfn_copy['V_pitch']) #X,Y this way because dfn dimension ordering

            if fill:
                ax.set_facecolor(colmap(np.amin(dfn_copy['dfn'])))
                mesh=ax.pcolormesh(E,V_pitch,dfn_copy['dfn'],cmap=colmap,vmin=vmin,vmax=vmax)            
            else:
                mesh=ax.contour(E,V_pitch,dfn_copy['dfn'],levels=np.linspace(vmin,vmax,num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                #ax.clabel(mesh,inline=1,fontsize=10)

            ax.set_xlabel('energy [eV]')
            ax.set_ylabel('pitch (v_parallel/v)')
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh            

        elif axes==['E','time']: #use the kwargs to accept TRANSP_output_FI_list, a list of other integrated TRANSP_output_FI objects to plot in a single mesh plot
            
            dfn_copy=self.dfn_integrate(energy=False)

            if 'TRANSP_output_FI_list' in kwargs.keys():

                E_list=[dfn_copy['dfn']]
                time_list=[dfn_copy['TIME']]
                
                for FI_object in kwargs['TRANSP_output_FI_list']:
                    E_list.append(FI_object['dfn'])
                    time_list.append(FI_object['TIME'])

                dfn_int_all=np.array(E_list,ndmin=2) #combine all objects in list into one
                E,time=np.meshgrid(dfn_copy['E'],np.array(time_list)) #automatically sorts in ascending time since pcolormesh does not require increasing meshgrid

                if vminmax:
                    vmin=vminmax[0]
                    vmax=vminmax[1]
                else:
                    vmin=np.amin(dfn_int_all)
                    vmax=np.amax(dfn_int_all)

                if fill:
                    ax.set_facecolor(colmap(np.amin(dfn_int_all)))
                    mesh=ax.pcolormesh(time,E,dfn_int_all,cmap=colmap,vmin=vmin,vmax=vmax)           
                else:
                    mesh=ax.contour(time,E,dfn_int_all,levels=np.linspace(vmin,vmax,num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    #ax.clabel(mesh,inline=1,fontsize=10)

                ax.set_xlabel('energy [eV]')                
                ax.set_ylabel('time [s]')
                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')
                    
                if ax_flag is True or fig_flag is True: #return the plot object
                    return mesh  

            else:
                print("ERROR: 'TRANSP_output_FI_list' not supplied in TRANSP_output_FI.dfn_plot(...**kwargs) - please see docstring for TRANSP_output_FI.dfn_plot()\n")

        #elif axes==['E']:

        #elif axes==['R']:

        elif len(axes)==self['dfn'].ndim: #assume user wants to plot energy pitch at point in real space

            dfn_copy=self['dfn'][tuple(axes)]
            E,V_pitch=np.meshgrid(self['E'],self['V_pitch']) #X,Y this way because dfn dimension ordering

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.amin(dfn_copy)
                vmax=np.amax(dfn_copy)

            if fill:          
                ax.set_facecolor(colmap(np.amin(dfn_copy)))
                mesh=ax.pcolormesh(E,V_pitch,dfn_copy,cmap=colmap,vmin=vmin,vmax=vmax)         
            else:
                mesh=ax.contour(E,V_pitch,dfn_copy,levels=np.linspace(vmin,vmax,num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                #ax.clabel(mesh,inline=1,fontsize=10)

            ax.set_xlabel('energy [eV]')
            ax.set_ylabel('pitch (v_parallel/v)')
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')

            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh 
                            
        else:
            print("ERROR: TRANSP_output_FI.dfn_plot() given unknown axes option - please check dfn_plot() docstring\n")

        if ax_flag is False and fig_flag is False:
            plt.show()

################################################################################################### ASCOT classes and functions

def dump_run_file_ASCOT(run_file='ascot4.cmd',initialdir=None,output_file='ascot.out',max_proc=50,min_proc=25,error_file='ascot.err',executable='test_ascot',tag=''):
    """
    writes out freia batch file for ASCOT run
    args:
        run_file - 
        initialdir - 
        output_file - 
        max_proc - 
        min_proc - 
        error_file - 
        executable - 
        tag - optional identifier tag for each set of run files produced
    notes:
    """

    filepath=support.dir_input_files+run_file+tag
    if initialdir is None:
        initialdir=os.path.dirname(os.path.abspath(__file__)) #use cwd
    
    with open(filepath,'w') as file:
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
        file.write("    RHO (pol)       Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)\n") #rho_pol is defined as sqrt((psi - psi_axis)/(psi_sep - psi_axis)), where psi is the poloidal magnetic flux, and psi_axis/psi_sep are its value evaluated on-axis/at the separatrix
 
        flux_pol_norm_sqrt=np.sqrt(np.abs(temperature_e['flux_pol_norm'])) #calculate profiles vs sqrt(flux_pol)
        flux_pol_norm_sqrt,te,ne,rot,ti,ni=processing.utils.sort_arrays(flux_pol_norm_sqrt,temperature_e['T'],density_e['n'],rotation_toroidal,temperature_i['T'],density_i['n']) #check order

        for RHO,Te,Ne,Vtor_I,Ti1,Ni1 in zip(flux_pol_norm_sqrt,te,ne,rot,ti,ni): 
            line=processing.utils.fortran_string(RHO,16,7)+processing.utils.fortran_string(Te,16,7)+processing.utils.fortran_string(Ne,16,7)+processing.utils.fortran_string(Vtor_I,15,7)+processing.utils.fortran_string(Ti1,17,7)+processing.utils.fortran_string(Ni1,15,7)+"\n"
            file.write(line)

    print("finished dumping profiles to ASCOT format")

def dump_wall_ASCOT(filename,wall):
    """
    dumps 2D wall outline to ASCOT input.wall_2d format
    
    notes:
        currently sets all divertor flags = 0 i.e. treats everything as 'wall'
    args:
        filename - output filename
        wall - data structure holding wall with same wall variable names as GEQDSK equilibrium e.g. rlim,zlim - i.e. a GEQDSK equilibrium object 

    """

    print("dumping wall to ASCOT format")

    filepath=support.dir_input_files+filename

    with open(filepath,'w') as file:

        file.write("{number_points} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n".format(number_points=int(wall['rlim'].size)))
        
        for r,z in zip(wall['rlim'],wall['zlim']):
            line=processing.utils.fortran_string(r,16,7)+processing.utils.fortran_string(z,16,7)+processing.utils.fortran_string(0.0,4,0,False)+"\n"
            file.write(line)

    print("finished dumping wall to ASCOT format")

def ASCOT_run_gen(temperature_i,temperature_e,density_i,density_e,rotation_toroidal,equilibrium,beam_deposition,wall,guiding_centre=True,tag=''):
    """
    generates full run input data for ASCOT

    notes:
        dumps everything into input_files folder
        uses default ASCOT filenames
        some inputs from list below are missing and must be generated by hand for now
        ASCOT run requires:
            - ascot4.cmd = freia batch submission file
            - input.magn_bkg = B field (currently generated with matlab scripts)
            - input.magn_header = B field (currently generated with matlab scripts)
            - input.options = run options namelist
            - input.particles = particle birth list
            - input.plasma_1d = kinetic profile inputs
            - input.wall_2d = limiter contour outline
            - binary = executable
    args:
        temperature_i - ion temperature object (eV)
        temperature_e - electron temperature object (eV)
        density_i - ion density object (#/m^3)
        density_e - electron density object (#/m^3)
        rotation_toroidal - array-like toroidal rotation (rad/s)
        equilibrium - equilibrium object 
        beam_deposition - beam deposition object
        wall - wall object
        guiding_centre - toggle dumping birth list at guiding-centre or particle position
        tag - optional identifier tag for each set of run files produced
    """

    print("ASCOT_run_gen creating ASCOT inputs")

    dump_run_file_ASCOT(initialdir=support.dir_input_files,tag=tag) #generate run file
    dump_profiles_ASCOT(filename='input.plasma_1d'+tag,temperature_i=temperature_i,temperature_e=temperature_e,density_i=density_i,density_e=density_e,rotation_toroidal=rotation_toroidal)
    dump_wall_ASCOT(filename='input.wall_2d'+tag,wall=wall)
    if guiding_centre:
        beam_deposition.dump_data(data_format='ASCOT_gc',filename='input.particles'+tag,equilibrium=equilibrium)
    else:
        beam_deposition.dump_data(data_format='ASCOT',filename='input.particles'+tag,equilibrium=equilibrium)

    print("ASCOT_run_gen finished")

class ASCOT_output:
    """
    class for encapsulating ASCOT HDF5 output file

    notes:
        mimics a LOCUST_IO object - use pull_data and methods like dfn_transform to then access standard LOCUST_IO functions
        my_output.file['key/path/to/data'].values will return leaf-level tree data from HDF5 file
    example:
        my_output=ASCOT_output('ID',some_filename,some_datatype)
        my_output.dfn_plot(axes=['R','Z'],real_scale=True) 
    """

    def __init__(self,ID,filename,datatype):
        """
        constructor for TRANSP_output class
        
        notes: 
        """

        self.ID=ID
        self.data={}
        self.filename=filename
        self.filepath=support.dir_output_files+filename
        self.datatype=datatype
        self.read_data(self.filename,self.datatype)

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

        notes:
        """

        self.data[key]=value

    def look(self):
        """
        print class information and data
        """

        print("\n-----------------------")
        print("ID - {ID}".format(ID=self.ID))  
        
        if hasattr(self,'filename'):
            print("Input Filename - {filename}".format(filename=self.filename))

        if hasattr(self,'properties') and self.properties:
            print("Properties:".format(properties=self.properties))
            for key in self.properties:
                if self.properties[key]: #do not print if the data is empty
                    print("{key} - {value}".format(key=key,value=self.properties[key])) 
        
        print("|")
        print("|")

        if hasattr(self,'data') and self.data:
            for key in self.data:
                if type(self.data[key])==type({}): #check for dicts since these mess things up
                    if self.data[key]: #if dict then check if dict is empty
                        print(key+":")
                        for sub_key in self.data[key]: 
                            print("     {sub_key} - {value}".format(sub_key=sub_key,value=self.data[key][sub_key]))
                elif not self.data[key].size==0: #if not a dict assume a numpy array
                    print("{key} - {value}".format(key=key,value=self.data[key]))
                
        print("-----------------------\n")

    def CDFlook(self,key=None):
        """
        prints file sub-branches from branch 'key'

        notes:
            can navigate tree with '/' like file directory
        usage:
            my_ASCOT_output.CDFlook() #print top level
            my_ASCOT_output.CDFlook('bfield')
        """

        if not self.file:
            print("no file member data found by .CDFlook()! use .file_open()")
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

    def file_open(self,filepath=None):
        """
        re-open/open new HDF5 file
        """
        
        try:
            import h5py
        except:
            raise ImportError("ERROR: ASCOT_output.file_open() could not import h5py module!\nreturning\n")
            return

        self.file=h5py.File(filepath,'r') 

    def file_close(self):
        """
        close annoying HDF5 file handle

        notes:
        """

        del(self.file)

    def read_data(self,filename=None,datatype='distribution_function'):
        """
        deep copy and pull data from HDF5 file 

        notes:
            extracts and maps out the output HDF5 file to LOCUST_IO variable names where possible
            see individual datatype branches for more detail
            NOTE work in progress

        args:
            datatype - read data from ASCOT output file and create structure resembling this LOCUST_IO datatype    
        """

        self.datatype=datatype #overwrite datatype member data
        if filename: #if supplied new filename, overwrite previous filepath
            self.filename=filename
            self.filepath=support.dir_output_files+filename
        self.file_open(self.filepath)

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
            self['V']=np.sqrt(2.*self['E']/mass_deuterium)

            self['dR']=np.abs(self['R'][1]-self['R'][0])
            self['dZ']=np.abs(self['Z'][1]-self['Z'][0])
            self['dV_pitch']=np.abs(self['V_pitch'][1]-self['V_pitch'][0])
            self['dE']=np.abs(self['E'][1]-self['E'][0])
            
            self['nR']=np.array(len(self['R']))
            self['nZ']=np.array(len(self['Z']))
            self['nV_pitch']=np.array(len(self['V_pitch']))
            self['nE']=np.array(len(self['E']))

            self['dfn']=self.file['distributions/rzPitchEdist/ordinate'].value #[m^-3 J^-1 dpitch^-1]
            self['dfn']=np.sum(self['dfn'],axis=0)
            self['dfn']=np.sum(self['dfn'],axis=0)
            self['dfn']=np.sum(self['dfn'],axis=-1)
            self['dfn']=np.swapaxes(self['dfn'],-1,-2) 
            self['dfn']*=e_charge #[m^-3 eV^-1 dpitch^-1]
            
            self['dfn_index']=np.array(['E','V_pitch','R','Z'])

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

        dfn_copy=copy.deepcopy(self)

        #begin list of specific options

        if axes==['R','Z']:
            dfn_copy['dfn']*=dfn_copy['dE']*dfn_copy['dV_pitch'] #integrate
            for counter in range(2): #sum
                dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=0)

        elif axes==['E','V_pitch']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn_copy['R'])):
                dfn_copy['dfn'][:,:,r,:]*=dfn_copy['R'][r]*2.*pi*dfn_copy['dR']*dfn_copy['dZ']
            #then need to integrate over the unwanted coordinates
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #over Z
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #over R

        elif axes==['E']:
            #applying real space Jacobian and integrate over toroidal angle
            for r in range(len(dfn_copy['R'])):
                dfn_copy['dfn'][:,:,r,:]*=dfn_copy['R'][r]*2.*pi*dfn_copy['dR']*dfn_copy['dZ']
            dfn_copy['dfn']*=dfn_copy['dV_pitch'] #integrate over pitch
            
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #sum over Z
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #sum over R
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #sum over V_pitch

        elif axes==['R']:
            dfn_copy['dfn']*=dfn_copy['dE']*dfn_copy['dV_pitch'] #integrate
            for counter in range(2): #sum over pitch and energy
                dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=0)
            dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=-1) #sum over Z

        elif axes==['N']:
            #applying full Jacobian and integrate over toroidal angle
            for r in range(len(dfn_copy['R'])):
                dfn_copy['dfn'][:,:,r,:]*=dfn_copy['R'][r]*2.*pi*dfn_copy['dR']*dfn_copy['dZ']*dfn_copy['dV_pitch']*dfn_copy['dE']
            for all_axes in range(dfn_copy['dfn'].ndim): #sum over all dimensions
                dfn_copy['dfn']=np.sum(dfn_copy['dfn'],axis=0) 
        
        #general option
        
        elif len(axes)==dfn_copy['dfn'].ndim: #if user supplies all axes then slice WITHOUT integrating
            dfn_copy['dfn']=dfn_copy['dfn'][tuple(axes)]
        else:
            print("ERROR: dfn_transform given invalid axes arguement: "+str(axes))

        return dfn_copy

    def dfn_plot(self,key='dfn',axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,transform=True,number_bins=20,fill=True,vminmax=None,ax=False,fig=False):
        """
        wrapper to plot_distribution_function

        notes:
        """

        #XXX - THIS IS FUDGE CODE COPIED FROM DFN.PLOT METHOD

        def plot_distribution_function(self,key='dfn',axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,transform=True,number_bins=20,fill=True,vminmax=vminmax,ax=False,fig=False):
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
                transform - set to False if supplied dfn has already been cut down to correct dimensions
                number_bins - set number of bins or levels
                fill - toggle contour fill on 2D plots
                vminmax - set mesh Vmin/Vmax values
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
                print(self[key])
                return
            
            #>0D data is plottable
            if fig_flag is False:
                fig = plt.figure() #if user has not externally supplied figure, generate
            
            if ax_flag is False: #if user has not externally supplied axes, generate them
                ax = fig.add_subplot(111)

            #1D data
            if self[key].ndim==1:
                ax.plot(self[axes[0]],self[key])
                ax.set_ylabel(key)
                ax.set_title(self.ID)

            #plot distribution function
            elif key=='dfn':
                
                #transform distribution function to the coordinates we want
                dfn_copy=copy.deepcopy(self)
                if transform is True:
                    dfn_copy=processing.process_output.dfn_transform(self,axes=axes) #user-supplied axes are checked for validity here

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
                    ax.plot(self[key])
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
                    Y,X=np.meshgrid(Y,X) #dfn is r,z so need to swap order here

                    if fill:
                        ax.set_facecolor(colmap(np.amin(dfn_copy[key])))
                        mesh=ax.pcolormesh(X,Y,dfn_copy[key],cmap=colmap,vmin=vmin,vmax=vmax)
                        #mesh=ax.contourf(X,Y,dfn_copy[key],levels=np.linspace(np.amin(dfn_copy[key]),np.amax(dfn_copy[key]),num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                        '''for c in mesh.collections: #for use in contourf
                            c.set_edgecolor("face")'''
                    else:
                        mesh=ax.contour(X,Y,dfn_copy[key],levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                        #ax.clabel(mesh,inline=1,fontsize=10)

                    if fig_flag is False:    
                        fig.colorbar(mesh,ax=ax,orientation='horizontal')
                    ax.set_xlabel(axes[0])
                    ax.set_ylabel(axes[1])
                    ax.set_title(self.ID)
                    
                    if real_scale is True: #set x and y plot limits to real scales
                        ax.set_aspect('equal')
                    else:
                        ax.set_aspect('auto')
                    if LCFS: #plot plasma boundary
                        ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],plot_style_LCFS) 
                    if limiters: #add boundaries if desired
                        ax.plot(limiters['rlim'],limiters['zlim'],plot_style_limiters)
                    
                    if ax_flag is True or fig_flag is True: #return the plot object
                        return mesh

                else: #user has not supplied >2D dfn
                    print("ERROR: plot_distribution_function given >2D DFN - please reduce dimensionality")
                    return 

            if ax_flag is False and fig_flag is False:
                plt.show()

        #XXX - THIS IS FUDGE CODE COPIED FROM DFN.PLOT METHOD

        if transform: 
            dfn_copy=self.dfn_transform(axes=axes)
        return plot_distribution_function(dfn_copy,key=key,axes=axes,LCFS=LCFS,limiters=limiters,real_scale=real_scale,colmap=colmap,transform=False,number_bins=number_bins,fill=fill,ax=ax,fig=fig) #call standard plot_distribution function but with LOCUST_IO version of transform disabled



################################################################################################### MARSF classes and functions

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
        flux_pol_norm_sqrt,rotation=processing.utils.sort_arrays(flux_pol_norm_sqrt,output_data['rotation']) #check order
 
        file.write("{length} {some_number}\n".format(length=int(flux_pol_norm_sqrt.size),some_number=1)) #re-insert line containing length
        
        for point in range(flux_pol_norm_sqrt.size): #iterate through all points i.e. length of our dictionary's arrays
            file.write("{flux_pol_norm_sqrt}{rotation}\n".format(flux_pol_norm_sqrt=processing.utils.fortran_string(flux_pol_norm_sqrt[point],24,18),rotation=processing.utils.fortran_string(rotation[point],24,18)))

    print("finished writing rotation to MARSF mogui")

'''
def read_coil_currents_pyuda(shot_number):
    """
    notes:
    args:
        shot_number - MAST shot number for signal
    """

    import pyuda
    client=pyuda.Client(shot_number)
'''

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


################################################################################################### LOCUST classes and functions

class FINT_LOCUST:
    """
    hacky class to read in LOCUST FINT data 

    notes:
    """
    
    def __init__(self,ID,filename='FINT.dat'): 
        """
        just reads in data from supplied filename
        """
        
        self.ID=ID
        self.data={}
        self.filename=filename
        self.filepath=support.dir_output_files+filename
        self.read_data(self.filepath)

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

    def read_data(self,filepath=None):
        """
        read data from filepath
        """

        with open(filepath,'r') as file:

            lines=file.readlines() #reading entire file anyway so grab all at once
            if not lines: #check to see if the file opened
                raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+filepath)


            self['nE']=float(lines[0])
            del(lines[0])
        
            self['E']=[] #start reading the energy spacings
            numbers_read=0
            while numbers_read<self['nE']:
                for line in lines:
                    split_line=line.split()
                    for number in split_line:
                        self['E'].extend([float(number)])
                        numbers_read+=1
                    del(lines[0])
                    break #must break out of loop here to re-check while condition

            self['E_norm']=float(lines[0])
            del(lines[0])
            self['Pdep']=float(lines[0])
            del(lines[0])
            
            #rest of file now structured the same in loops
            self['time']=[]
            self['PFC_power']=[]
            self['dfn']=[]
            while lines: 
                self['time'].extend([float(lines[0])])
                del(lines[0])
                self['PFC_power'].extend([float(lines[0])])
                del(lines[0])
                numbers_read=0
                while numbers_read<self['nE']:
                    for line in lines:
                        split_line=line.split()
                        for number in split_line:
                            self['dfn'].extend([float(number)])
                            numbers_read+=1
                        del(lines[0])
                        break #must break out of loop here to re-check while condition

            self['nE']=np.asarray(self['nE'])
            self['E']=np.asarray(self['E'])
            self['E_norm']=np.asarray(self['E_norm'])
            self['Pdep']=np.asarray(self['Pdep'])
            self['time']=np.asarray(self['time'])
            self['PFC_power']=np.asarray(self['PFC_power'])
            self['dfn']=np.asarray(self['dfn'])
            self['dfn']=self['dfn'].reshape(int(len(self['time'])),int(self['nE'])) #XXX check order of this

    def dfn_plot(self,axes=['E','time'],colmap=cmap_default,number_bins=20,ax=False,fig=False):
        """
        plot dfn vs time
        
        note:
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

        E,time=np.meshgrid(self['E'],self['time'])

        if fill:
            ax.set_facecolor(colmap(np.amin(self['dfn'])))
            mesh=ax.pcolormesh(time,E,self['dfn'],cmap=colmap,vmin=np.amin(self['dfn']),vmax=np.amax(self['dfn']))
        else:
            mesh=ax.contour(time,E,self['dfn'],levels=np.linspace(np.amin(self['dfn']),np.amax(self['dfn']),num=number_bins),colours=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(self['dfn']),vmax=np.amax(self['dfn']))
            #ax.clabel(mesh,inline=1,fontsize=10)

        ax.set_xlabel('time [s]')
        ax.set_ylabel('energy [eV]')                
        ax.set_title(self.ID)
        
        if fig_flag is False:    
            fig.colorbar(mesh,ax=ax,orientation='horizontal')

        if ax_flag is True or fig_flag is True: #return the plot object
            return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()

################################################################################################### misc functions



#################################
 
##################################################################
 
###################################################################################################
