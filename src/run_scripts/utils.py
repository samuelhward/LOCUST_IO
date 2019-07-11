#run_scripts.utils.py
 
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
    import pathlib
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
    import classes.input_classes.equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/equilibrium.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    import classes.input_classes.temperature
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/temperature.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    import classes.input_classes.number_density
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/number_density.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    import classes.input_classes.beam_deposition
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/beam_deposition.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.wall
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/wall.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
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

def TRANSP_get_fbm_FI_CDF(run_ID,shot_number,number_files,particle_position=True,guiding_centre=True,path_TRANSP='',device='d3d'):
    """
    notes:
        looks for files in output_files
    args:
        run_ID - TRANSP run_ID e.g. W01
        shot_number - TRANSP shot number e.g. 29034
        number_files - total number of .DATA# files to extract CDF from 
        particle_position - toggle whether to generate set of CDFs at particle positions
        guiding_centre - toggle whether to generate set of CDFs at guiding centres
        path_TRANSP - path to TRANSP files in input_files dir (input_files/path_TRANSP...)
        device - device code for machine under study
    """

    ID=str(shot_number)+str(run_ID)

    project_dir=os.getcwd()
    target_dir=os.path.join(support.dir_output_files,path_TRANSP)
    os.chdir(target_dir) #change working directory to output files briefly due to bug in CCFE get_fbm installation

    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(ID,'_fi_',file_ID,'.cdf') #current output file
        output_filename_gc='{}{}{}{}'.format(ID,'_fi_',file_ID,'_gc.cdf')

        if guiding_centre:
            fbm_input = """
                        {ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        c
                        """.format(ID=ID,path='q',file_ID=file_ID,device=device)

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
                        {ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        p
                        """.format(ID=ID,path='q',file_ID=file_ID,device=device)

            print("writing TRANSP FI netCDF file {}".format(output_filename))
            proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out, err = proc.communicate(input=fbm_input.encode('utf-8'))
            print("finished writing TRANSP FI netCDF file {}".format(output_filename))

    os.chdir(project_dir) #change back to original working directory    

def TRANSP_get_fbm_FI_birth_deposition(run_ID,shot_number,number_files,path_TRANSP='',device='d3d'):
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
        run_ID - TRANSP run_ID e.g. W01
        shot_number - TRANSP shot number e.g. 29034
        number_files - total number of .DATA# files to extract CDF from 
        path_TRANSP - path to TRANSP files in input_files dir (input_files/path_TRANSP...)
        device - device code for machine under study
    """

    ID=str(shot_number)+str(run_ID)

    project_dir=os.getcwd()
    target_dir=os.path.join(support.dir_output_files,path_TRANSP)
    os.chdir(target_dir) #change working directory to output files briefly due to bug in CCFE get_fbm installation
    
    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(ID,'_fdep_nb_en_',file_ID,'.out') #current output file

        fbm_input = """
                    {ID}
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
                    """.format(ID=ID,path='q',file_ID=file_ID,device=device,output_filename=output_filename)

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
        self.filepath=support.dir_output_files / filename
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
        self.filepath=support.dir_output_files / filename
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
        self.filepath=support.dir_output_files / filename
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
        self['dfn_index']=np.array(['RZ','V_pitch','E'])

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
            try tricontourf for TRANSP plotting instead of interpolation? unstructured grid contours!
                https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.tricontour.html
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
            V_pitch - [dPitch]^-1
            R - [m]^-3 x
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
        ax.set_title(self.ID)   
        
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
                vmin=np.amin(new_dfn)
                vmax=np.amax(new_dfn)

            if fill:
                ax.set_facecolor(colmap(np.amin(new_dfn)))
                mesh=ax.pcolormesh(R,Z,new_dfn,cmap=colmap,vmin=vmin,vmax=vmax)
            else:
                mesh=ax.contour(R,Z,new_dfn,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
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
                mesh=ax.contour(E,V_pitch,dfn_copy['dfn'],levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
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
                    mesh=ax.contour(time,E,dfn_int_all,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
                    #ax.clabel(mesh,inline=1,fontsize=10)

                ax.set_xlabel('energy [eV]')                
                ax.set_ylabel('time [s]')
                if fig_flag is False:    
                    fig.colorbar(mesh,ax=ax,orientation='horizontal')
                    
                if ax_flag is True or fig_flag is True: #return the plot object
                    return mesh  

            else:
                print("ERROR: 'TRANSP_output_FI_list' not supplied in TRANSP_output_FI.dfn_plot(...**kwargs) - please see docstring for TRANSP_output_FI.dfn_plot()\n")

        elif axes==['E']: #integrate over all volume and plot as a function of energy in #/eV

            dfn_copy=self.dfn_integrate(energy=False)
            ax.plot(dfn_copy[axes[0]],dfn_copy['dfn'],color=colmap(np.random.uniform()))
            ax.set_xlabel('energy [eV]')
            ax.set_ylabel('density [#/eV]')

        elif axes==['V_pitch']:

            dfn_copy=self.dfn_integrate(pitch=False)
            ax.plot(dfn_copy[axes[0]],dfn_copy['dfn'],color=colmap(np.random.uniform()))
            ax.set_xlabel('pitch [V||/V]')
            ax.set_ylabel('density [#/dPitch]')

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
                mesh=ax.contour(E,V_pitch,dfn_copy,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=vmin,vmax=vmax)
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

def read_inputs_TRANSP(run_ID,shot_number,input_path=pathlib.Path(''),beam_depo_GC=True,beam_depo_number=None,GEQDSKFIX=0):
    """
    reads full input_data from TRANSP run

    notes:
        assumes some default filenames based on the TRANSP run_ID
        equilibrium must be stored in format shot_number+run_ID.GEQDSK
        assumes _birth.cdf files are in input_files/ despite being TRANSP outputs
        XXX reading UFILE data (e.g. for kinetic profiles) will currently break if more than one timestep is stored
        XXX only reads first beam deposition file
    args:
        run_ID - TRANSP run_ID e.g. W01
        shot_number - TRANSP shot number e.g. 29034
        input_path - path to target in input_files dir (input_files/path/)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        beam_depo_number - integer number of beam depo file to read elif None then read all available beam depositions and combine into single object 
        GEQDSKFIX - LOCUST-equivalent flag to optionally flip fields in GEQDSK
    """
    
    print("read_inputs_TRANSP()")

    filepath_temperature_i_1=input_path / 'OMF'+shot_number+'.TIO' 
    filepath_temperature_i_2=input_path / 'OMF'+shot_number+'.TI2'    
    filepath_temperature_e_1=input_path / 'OMF'+shot_number+'.TEL' 
    filepath_temperature_e_2=input_path / 'OMF'+shot_number+'.TER' 
    filepath_number_density_e_1=input_path / 'OMF'+shot_number+'.NEL' 
    filepath_number_density_e_2=input_path / 'OMF'+shot_number+'.NER' 
    filepath_equilibrium=input_path / 'g'+shot_number
    if beam_depo_number:
        filepaths_beam_deposition=list(input_path / shot_number+run_ID+'_birth.cdf{}'.format(str(beam_depo_number))) #list is just one file long
    else:
        filepaths_beam_deposition=list(pathlib.Path(support.dir_input_files / input_path).glob('*_birth.cdf*')) #find all birth CDF files in supplied input_path directory  
    filepath_wall=input_path / 'OMF'+shot_number+'.LIM'

    if beam_depo_GC:
        data_format_beam_depo='TRANSP_birth_gc'
    else:
        data_format_beam_depo='TRANSP_birth'

    try:
        temperature_i=classes.input_classes.temperature.Temperature(ID=shot_number+run_ID,data_format='UFILE',species='ions',filename=filepath_temperature_i_1)
    except:
        try:
            temperature_i=classes.input_classes.temperature.Temperature(ID=shot_number+run_ID,data_format='UFILE',species='ions',filename=filepath_temperature_i_2)            
        except:
            print("WARNING: read_inputs_TRANSP() could not read ion temperature from LOCUST_IO/input_files/{} - returning None".format(str(filepath_temperature_i_1)+" or "+str(filepath_temperature_i_2)))
            temperature_i=None
    try:
        temperature_e=classes.input_classes.temperature.Temperature(ID=shot_number+run_ID,data_format='UFILE',species='electrons',filename=filepath_temperature_e_1)
    except:
        try:
            temperature_e=classes.input_classes.temperature.Temperature(ID=shot_number+run_ID,data_format='UFILE',species='electrons',filename=filepath_temperature_e_2)
        except:
            print("WARNING: read_inputs_TRANSP() could not read electron temperature from LOCUST_IO/input_files/{} - returning None".format(str(filepath_temperature_e_1)+" or "+str(filepath_temperature_e_2)))
            temperature_e=None
    try:
        density_e=classes.input_classes.number_density.Number_Density(ID=shot_number+run_ID,data_format='UFILE',species='electrons',filename=filepath_number_density_e_1)
    except:
        try:
            density_e=classes.input_classes.number_density.Number_Density(ID=shot_number+run_ID,data_format='UFILE',species='electrons',filename=filepath_number_density_e_2)            
        except:
            print("WARNING: read_inputs_TRANSP() could not read electron density from LOCUST_IO/input_files/{} - returning None".format(str(filepath_number_density_e_1)+" or "+str(filepath_number_density_e_2)))
            density_e=None
    try:
        equilibrium=classes.input_classes.equilibrium.Equilibrium(ID=shot_number+run_ID,data_format='GEQDSK',filename=filepath_equilibrium,GEQDSKFIX=GEQDSKFIX)
    except:
        print("WARNING: read_inputs_TRANSP() could not read equilibrium from LOCUST_IO/input_files/{} - returning None".format(filepath_equilibrium))
        equilibrium=None

    beam_deposition=classes.input_classes.beam_deposition.Beam_Deposition(ID=shot_number+run_ID) #generate blank beam deposition which we will append to using .combine
    for filepath_beam_deposition in filepaths_beam_deposition:
        try:
            beam_deposition.combine(classes.input_classes.beam_deposition.Beam_Deposition(ID=shot_number+run_ID,data_format=data_format_beam_depo,filename=filepath_beam_deposition))
        except:
            print("WARNING: read_inputs_TRANSP() could not read beam deposition from LOCUST_IO/input_files/{} - returning None".format(filepath_beam_deposition))

    try:
        wall=classes.input_classes.wall.Wall(ID=shot_number+run_ID,data_format='UFILE',filename=filepath_wall)
    except:
        print("WARNING: read_inputs_TRANSP() could not read wall from LOCUST_IO/input_files/{} - returning None".format(filepath_wall))
        wall=None

    print("finished read_inputs_TRANSP()")

    return temperature_i,temperature_e,density_e,equilibrium,beam_deposition,wall

################################################################################################### ASCOT classes and functions

def dump_run_file_ASCOT(run_file='ascot4.cmd',initialdir=None,output_file='ascot.out',max_proc=50,min_proc=25,error_file='ascot.err',executable='test_ascot',input_path=pathlib.Path(''),tag='',user='sward'):
    """
    writes out freia batch file for ASCOT run
    args:
        run_file - batch file name 
        initialdir - directory holding ASCOT run input, run file is dumped here
        output_file - direct STDOUT to this file  
        max_proc - maximum number of processors
        min_proc - minimum number of processors
        error_file - direct STDERR to this file
        executable - filename of executable binary
        input_path - path to target in input_files dir (input_files/path/)
        tag - optional identifier tag for each set of run files produced
        user - username for email notifications
    notes:
    """

    print("dumping ASCOT run file")

    filepath=support.dir_input_files / input_path / str(run_file+tag)
    if initialdir is None:
        initialdir=support.dir_input_files / input_path #use write location as default
    
    with open(filepath,'w') as file:
        file.write('# @ input = dev/null/\n')
        file.write('# @ initialdir = {initialdir}\n'.format(initialdir=initialdir))
        file.write('# @ output = {output_file}\n'.format(output_file=output_file))
        file.write('# @ error = {error_file}\n'.format(error_file=error_file))
        file.write('# @ jobtype = openmpi\n')
        file.write('# @ max_processors = {max_proc}\n'.format(max_proc=max_proc))
        file.write('# @ min_processors = {min_proc}\n'.format(min_proc=min_proc))
        file.write('# @ notify_user = {user}\n'.format(user=user))
        file.write('# @ notification = complete\n')
        file.write('# @ queue\n\n')
        file.write('date\n')
        file.write('mpirun -np $NSLOTS {executable} -output ascot_freia_$JOB_ID\n'.format(executable=executable))
        file.write('date\n')

    print("finished dumping ASCOT run file")

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
 
    filepath=support.dir_input_files / filename
    with open(filepath,'w') as file:

        file.write("# Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation\n")
        file.write("# range must cover [0,1] of normalised poloidal rho. It can exceed 1. \n")
        file.write("# 18Jan08 for testing (first 3 lines are comment lines)\n")
        file.write(str(temperature_e['flux_pol_norm'].size)+"   "+"1"+" # Nrad,Nion\n")
        file.write("1           # ion Znum\n")
        file.write("2           # ion Amass\n")
        file.write("1 1         # collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons\n")
        file.write("    RHO (pol)       Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)\n") #rho_pol is defined as np.sqrt((psi - psi_axis)/(psi_sep - psi_axis)), where psi is the poloidal magnetic flux, and psi_axis/psi_sep are its value evaluated on-axis/at the separatrix
 
        flux_pol_norm_sqrt=np.sqrt(np.abs(temperature_e['flux_pol_norm'])) #calculate profiles vs np.sqrt(flux_pol)
        flux_pol_norm_sqrt,te,ne,rot,ti,ni=processing.utils.sort_arrays(flux_pol_norm_sqrt,temperature_e['T'],density_e['n'],rotation_toroidal,temperature_i['T'],density_i['n']) #check order

        for RHO,Te,Ne,Vtor_I,Ti1,Ni1 in zip(flux_pol_norm_sqrt,te,ne,rot,ti,ni): 
            line=processing.utils.fortran_string(RHO,16,7)+processing.utils.fortran_string(Te,16,7)+processing.utils.fortran_string(Ne,16,7)+processing.utils.fortran_string(Vtor_I,15,7)+processing.utils.fortran_string(Ti1,17,7)+processing.utils.fortran_string(Ni1,15,7)+"\n"
            file.write(line)

    print("finished dumping profiles to ASCOT format")

def dump_input_options_ASCOT(filename='input.options'):
    """
    notes:
    """

    print("dumping ASCOT input options")

    filepath=support.dir_input_files / filename
    with open(filepath,'w') as file:
        file_contents="""
         !  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

         !!! This is an input.options file for ASCOT4.     !!!
         !!! You may freely alter the values listed here.  !!!
         
         !!! Options related to orbit trajectory integration:
         !   real(kind=params_wp) :: ORB%RUNGEKUTTATOL:
         !    - a relative tolerance for adaptive runge-kutta
         !      integration methods.
         !   real(kind=params_wp) :: ORB%WRITINGINTERVAL:
         !    - a frequency for writing the particle info into the hard drive.
         !    - Negative value means that only the starting values and endstate
         !      are written into the file.
         !    - zero, means that every value, the orbit integrator produces
         !      is written into the file, and logically positive value
         !      means the frequency for writing. The default value is -1.0.
         !   real(kind=params_wp) :: ORB%MAXDELTAPHIPER2PI:
         !    - an estimate for the maximum allowed change in the phi/2pi,
         !      i.e, the fraction of the toroidal propagation.
         !    - restricts the length of the time step. The default value is 1.0
         !   real(kind=params_wp) :: ORB%MAXDELTARHO:
         !    - an estimate for the maximum allowed change in rho,
         !      i.e, radial propagation.
         !    - restricts the length of the time step. The default value is 1.0
         !   integer :: ORB%LEAPFROGSTEPDIVIDOR:
         !    - the leap-frog (and runge-kutta) time step is defined as
         !      dt = m /( B * q * leapFrogStepDividor )
         !      when fixedTimeStep is negative
         !   integer :: ORB%ORBITMETHOD:
         !    - 1 == pure guiding-center runge-kutta integration
         !    - 2 == full-orbit leap-frog integration (conserves energy)
         !    - 3 == full-orbit runge-kutta integration (accurate on place)
         !    - 4 == guiding-center runge-kutta with full-orbit
         !           runge-kutta wall collision checks
         !   integer :: ORB%ERAD:
         !    - 0 == no Erad profile from file/cpo
         !    - 1 == read Erad profile from file/cpo
         !   real(kind=params_wp) :: ORB%ALD_FORCE:
         !    - radiation reaction force. Relevant only for (runaway) electrons.
         !    - defines the strength of the force: <0 = off (default), 1 = physical.
         !    - the force is proportional to B^2 so adjusting this parameter shows
         !      what happens if the field strength is altered.
         !   real(kind=params_wp) :: ORB%FIXEDTIMESTEP:
         !    - Fixed orbit integration time step in seconds.
         !    - Negative value for adaptive time step
         !      and to use leapFrogStepDividor.
         !   real(kind=params_wp) :: ORB%GLOBALTIMESTEP:
         !    - Time step for ensemble time stepping mode (all particles stepwise)
         !    - Zero value disables ensemble time stepping mode.
        &ORBIT_TRACKING_OPTIONS
         ORB%RUNGEKUTTATOL=  2.0000000000000000E-006,
         ORB%WRITINGINTERVAL=  -1.0000000000000000,
         ORB%MAXDELTAPHIPER2PI=  2.0000000000000000E-001,
         ORB%MAXDELTARHO=  5.0000000000000000E-002,
         ORB%LEAPFROGSTEPDIVIDOR=         10,
         ORB%ORBITMETHOD=          1,
         ORB%ERAD=          0,
         ORB%ALD_FORCE= -1.0000000000000000     ,
         ORB%FIXEDTIMESTEP= -1,
         ORB%GLOBALTIMESTEP= 0.0,
         /
         
         !!! Options related to interactions with background:
         !   integer :: INTERACT%coulombCollMode:
         !    - tells what coulomb collision model is in use
         !    - 0 == no coulomb collisions
         !    - 1 == weakly relativistic coulomb collision kernel
         !   integer :: INTERACT%mhd:
         !    - tells whether mhd modes are being used
         !    - 0 == no mhd modes
         !    - 1 == mhd modes according to april 2012, input.mhd
         !   integer :: INTERACT%atomic:
         !    - tells whether effective ionization and recombination are in use
         !    - 0 == not in use
         !    - 1 == in use
         !   integer :: INTERACT%flow:
         !    - tells whether plasma flow/rotation is in use
         !    - 0 == not in use
         !    - 1 == toroidal rotation in [rad/s]
         !    - 2 == parallel flow in [m/s]
         !   integer :: INTERACT%ACC:
         !    - Interaction timescale acceleration toggle, off by default
         !    - 0 == no acceleration
         !    - 1 == adaptive acceleration, NOTE: set ENERGYCOLLSTEPLIM to a
         !           reasonable value, each time step will try to reach it.
         !           Reasonable = probably much lower than the default value.
         !    ->1 == max acceleration factor is user defined (this parameter)
         !   integer :: INTERACT%ICRH:
         !    - ICRH heating (In development stages)
         !    - 0 == no icrh heating (default)
         !    - 1 == Self-consistent ICRH heating with RFOF according to
         !           parameters in rfof_codeparam.xml and rfof_codeparam.xsd
         !   real(kind=params_wp) :: INTERACT%ICRHT0:
         !    - ICRH switch-on time (s)
         !   real(kind=params_wp) :: INTERACT%PROF_EXTRAP:
         !    - sets the characteristic width of extrapolated plasma 1D-profiles
         !      in units of rho when this data is not available all the way to 
         !      the wall and the outermost data points are non-zero
         !    - Zero means that 1D-profile extrapolation is disabled
         !    - Negative values means that 1D-profiles are extrapolated by an 
         !      exponential fit matching the derivative at the last data point
         !   real(kind=params_wp) :: INTERACT%PITCHCOLLTIMESTEPLIM:
         !    - sets limit for collision timestep in terms of pitch collisions.
         !    - timestep is proportional to this parameter.
         !    - squareroot of this parameter is the maximum stochastic
         !    - variation of the pitch. The default value is 0.003.
         !   real(kind=params_wp) :: INTERACT%FIXEDTIMESTEP:
         !    - Fixed collision time step in seconds.
         !    - Negative value for adaptive time step.
        &INTERACTION_OPTIONS
         INTERACT%COULOMBCOLLMODE=          1,
         INTERACT%MHD=          0,
         INTERACT%ATOMIC=          0,
         INTERACT%ACC=          20,
         INTERACT%FLOW=          0,
         INTERACT%ICRH=          0,
         INTERACT%ICRHT0=        0.0E0,
         INTERACT%PROF_EXTRAP=   0.0,
         INTERACT%PITCHCOLLTIMESTEPLIM=  1.0000000000000000E-003,
         INTERACT%FIXEDTIMESTEP= -1,
         /
         
         !!! Options defining the end conditions for the particle:
         !   real(kind=params_wp) :: ENDCOND%TMAX:
         !    - Maximum time for the single particle simulation.
         !    - If this variable is smaller than the possible field
         !      tmax in "input.particles", then this variable dominates.
         !    - the particle%tmax is defined to be:
         !      min(particle%tmax,ENDCOND%TMAX)
         !   real(kind=params_wp) :: ENDCOND%MINKINETICENERGY:
         !    - If the particles kinetic energy gets below this value,
         !      then the simulation ends.
         !   real(kind=params_wp) :: ENDCOND%:TIMESLOCALTHERMALENERGY
         !    - If the particle energy gets below the local thermal
         !      energy times this variable, then the simulation ends.
         !   real(kind=params_wp) :: ENDCOND%:CPUTMAX
         !    - Maximum CPU time per particle.
         !   real(kind=params_wp) :: ENDCOND%:RHOMAX
         !    - Maximum value of rho (1.0 for separatrix).
        &END_CONDITIONS
         ENDCOND%TMAX=  0.100000000000000     ,
         ENDCOND%MINKINETICENERGY=  50.000000000000     ,
         ENDCOND%TIMESLOCALTHERMALENERGY=  1.5000000000000000     ,
         ENDCOND%CPUTMAX=  10000000.0000000000000     ,
         ENDCOND%RHOMAX=  10000000.0000000     ,
         /
         
         !!! Options to decide what particle features to write into output:
         !   integer :: anum:
         !    - the mass number.
         !   integer :: mass:
         !    - the particle mass [amu].
         !   integer :: znum:
         !    - the particle nucleus charge [e]
         !   integer :: charge:
         !    - the particle charge [e]
         !   integer :: energy:
         !    - particle energy [eV].
         !   integer :: pitch:
         !    - particle pitch angle cosine (vpar/v) [].
         !   integer :: R:
         !    - cylindrical coordinate R [m] for the particle.
         !   integer :: phi:
         !    - cylindrical coordinate phi [rad] for the particle.
         !   integer :: Z:
         !    - cylindrical coordinate z [m] for the particle.
         !   integer :: Id:
         !    - particle identification number.
         !   integer :: rho:
         !    - sqrt of normalized poloidal flux []
         !   integer :: Time:
         !    - the time how long the particle has been followed [s].
         !   integer :: Bfield:
         !    - magnetic field components [T].
         !   integer :: cputime:
         !    - the elapsed cputime [s].
         !   integer :: dtLimiters:
         !    - which time step limiters (for orbit following and collisions) were strictest.
        &PARTICLE_DATA_FIELDS
         FIELDS%ANUM=          0,
         FIELDS%MASS=          0,
         FIELDS%ZNUM=          0,
         FIELDS%CHARGE=          0,
         FIELDS%ENERGY=          1,
         FIELDS%PITCH=          1,
         FIELDS%PHI=          1,
         FIELDS%R=          1,
         FIELDS%Z=          1,
         FIELDS%ID=          1,
         FIELDS%RHO=          1,
         FIELDS%TIME=          1,
         FIELDS%ENDCOND=          0,
         FIELDS%FBOUNCE=          0,
         FIELDS%FORBIT=          0,
         FIELDS%FTORPREC=          0,
         FIELDS%ORBITTYPE=          0,
         FIELDS%PPHI=          0,
         FIELDS%BFIELD=          0,
         FIELDS%CPUTIME=          0,
         FIELDS%DTLIMITERS=          0,
         /
         
         !!! Options for distributions
         !   integer: simpleDists
         !    - add whole time step (0) or only the end point (1) to distributions
         !   integer: rhodists
         !    - [rho,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzdists
         !    - [r,z,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzvdist
         !    - [r,z,vpar,vperp,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzpdists
         !    - [r,z,ppar,pperp,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzpitchedist
         !    - [r,z,pitch,energy,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzmuedist
         !    - [r,z,magnetic moment,energy,time,test particle species] -distributions off(0) or on(1) or with variance(2)
         !   integer: rhophipitchEdist
         !    - [rho,phi,pitch,energy,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: rzvvdist
         !    - [r,z,vtor/vtot,vz/vtot,time,test particle species] - distributions off (0) or on (1) or with variance (2)
         !   integer: nrho
         !    - number of rho slots in rhodists and rho-phi-E-pithc dists
         !   integer: nphi
         !    - number of phi slots in rho and rho-phi-E-pithc dists
         !   integer: ntime
         !    - number of time slots in all distributions
         !   real: rho1
         !    - lower bound of rho in rho and rho-phi-E-pithc dists
         !   real: rho2
         !    - upper bound of rho in rho and rho-phi-E-pithc dists
         !   real: phi1
         !    - lower bound of phi in rho-phi-E-pithc dist (deg)
         !   real: phi2
         !    - upper bound of phi in rho-phi-E-pithc dist (deg)
         !   real: time1
         !    - lower bound of time in all distributions
         !   real: time2
         !    - upper bound of time in all distributions, endcond%tmax is used if not given
         !   integer: nr
         !    - number of r slots in all rz - distributions
         !   integer: R_auto
         !    - automatically set R bounds from wall input
         !   integer: nz
         !    - number of z slots in all rz - distributions
         !   integer: z_auto
         !    - automatically set z bounds from wall input
         !   integer: nvpa
         !    - number of parallel velocity slots in rzvdists
         !   integer: vpa_auto
         !    - automatically set vpa bounds from particle input
         !   integer: nvpe
         !    - number of perpendicular velocity slots in rzvdists
         !   integer: vpe_auto
         !    - automatically set vpe bounds from particle input
         !   integer: nppa
         !    - number of parallel momentum slots in rzpdists
         !   integer: nppe
         !    - number of perpendicular momentum slots in rzpdists
         !   real: r1
         !    - lower bound of r in all rz - distributions
         !   real: r2
         !    - upper bound of r in all rz - distributions
         !   real: z1
         !    - lower bound of z in all rz - distributions
         !   real: z2
         !    - upper bound of z in all rz - distributions
         !   real: vpa1
         !    - lower bound of parallel velocity in rzvdists
         !   real: vpa2
         !    - upper bound of parallel velocity in rzvdists
         !   real: vpe1
         !    - lower bound of perpendicular velocity in rzvdists
         !   real: vpe2
         !    - upper bound of perpendicular velocity in rzvdists
         !   real: ppa1
         !    - lower bound of parallel momentum in rzpdists
         !   real: ppa2
         !    - upper bound of parallel momentum in rzpdists
         !   real: ppe1
         !    - lower bound of perpendicular momentum in rzpdists
         !   real: ppe2
         !    - upper bound of perpendicular momentum in rzpdists
         !   integer: npitch
         !    - number of pitch slots in rzpitchedist (pitch bounds will always be [-1,1])
         !   integer: nenergy
         !    - number of energy slots in rzpitchedist
         !   integer: energy_auto
         !    - automatically set energy bounds from particle input
         !   real: energy1
         !    - lower bound of energy in rzpitchedist [eV]
         !   real: energy2
         !    - upper bound of energy in rzpitchedist [eV]
         !   integer: nmu
         !    - number of mu slots in rzmuedist
         !   real: mu1
         !    - lower bound of mu in rzmuedist [eV/T]
         !   real: mu2
         !    - upper bound of mu in rzmuedist [eV/T]
         !   integer: timedepmode
         !    - Time dependence mode: (0) Steady-state, (1) time-dependent using CPOs
        &DIST_OPTIONS
         DIST%SIMPLEDISTS=        0,
         DIST%RHODISTS=           2,
         DIST%RZDISTS=            0,
         DIST%RZVDIST=            1,
         DIST%RZPDIST=            0,
         DIST%RZPITCHEDIST=       1,
         DIST%RZVVDIST=           0,
         DIST%RZMUEDIST=          0,
         DIST%RHOPHIPITCHEDIST=   0,
         DIST%NTIME=              1,
         DIST%TIME1=  0.0000000000000000     ,
         DIST%TIME2=  10000000.00000000000000     ,
         DIST%NRHO=             200,
         DIST%RHO1=  0.0000000000000000     ,
         DIST%RHO2=  1.0000000000000000     ,
         DIST%NR=                100,
         DIST%R_AUTO=             0,
         DIST%R1=  1.000     ,
         DIST%R2=  2.500     ,
         DIST%NPHI=               1,
         DIST%PHI1=  0.0000000000000000     ,
         DIST%PHI2=  360.00000000000000     ,
         DIST%NZ=                100,
         DIST%Z_AUTO=             0,
         DIST%Z1= -1.45     ,
         DIST%Z2=  1.45     ,
         DIST%NVPA=              30,
         DIST%VPA_AUTO=           0,
         DIST%VPA1= -1.5e6     ,
         DIST%VPA2=  1.5e6     ,
         DIST%NVPE=              30,
         DIST%VPE_AUTO=           0,
         DIST%VPE1=  0.0     ,
         DIST%VPE2=  1.5e6     ,
         DIST%NPPA=              10,
         DIST%PPA1= -5.9999999999999998E-021,
         DIST%PPA2=  5.9999999999999998E-021,
         DIST%NPPE=              10,
         DIST%PPE1=  1.0000000000000000E-022,
         DIST%PPE2=  5.9999999999999998E-021,
         DIST%NENERGY=           100,
         DIST%ENERGY_AUTO=        0,
         DIST%ENERGY1=  0.000000000000     ,
         DIST%ENERGY2=  100000.0000000000     ,
         DIST%NMU=               30,
         DIST%MU1=  10000.000000000000     ,
         DIST%MU2=  120000.0000000000     ,
         DIST%NVTOR=             40,
         DIST%NVZ=               20,
         DIST%NPITCH=            100,
         DIST%TIMEDEPMODE=          1,
         /

         !!! Options related to particle source:
         !   integer :: PARTICLESOURCE%READ:
         !    - read test particles from somewhere yes/no (1/0)
         !   integer :: PARTICLESOURCE%DISTRIBUTION:
         !    - sample test particles from given 4D distrib. yes/1st step/no (1/2/0)
         !   integer :: PARTICLESOURCE%DYNAMNBI:
         !    - dynamic NBI particle generation on/off (1/0)
         !   integer :: PARTICLESOURCE%REWEIGHTING:
         !    - particle reweighting on/off (1/0)
         !   integer :: PARTICLESOURCE%REMLOWWGT:
         !    - reweighting: remove particles with low weight yes/no (1/0)
         !   integer :: PARTICLESOURCE%NRHOBINS:
         !    - number of rho bins for reweighting
         !   integer :: PARTICLESOURCE%NMAX:
         !    - number of test particles to aim at when using dynamic source
        &PARTICLESOURCE_OPTIONS
         PARTICLESOURCE%READ = 1,
         PARTICLESOURCE%DISTRIBUTION = 0,
         PARTICLESOURCE%DYNAMNBI = 1,
         PARTICLESOURCE%REWEIGHTING = 1,
         PARTICLESOURCE%REMLOWWGT = 1,
         PARTICLESOURCE%NRHOBINS = 10,
         PARTICLESOURCE%NMAX = 1000,
         /"""

        file.write(file_contents)

    print("finished dumping ASCOT input options")

def read_inputs_ASCOT(input_path=pathlib.Path(''),beam_depo_GC=True,species_numbers=[1],wall_type='2D'):
    """
    reads full run input data from ASCOT inputs

    notes:
        XXX currently no way to asssimilate ASCOT equilibrium so ignore for now
        currently assumes 2D wall
        assumes all kinetic profiles for all desired species stored in same location
        uses default ASCOT filenames
        takes wall from ASCOT input, not output
    args:
        input_path - path to target in input_files dir (input_files/path/)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        species_numbers - species number labels from input particle list to read in
        wall_type - set wall type to '2D' or '3D'
    """

    print("read_inputs_ASCOT()")

    filepath_temperature_i=input_path / 'input.plasma_1d' 
    filepath_density_i=input_path / 'input.plasma_1d' 
    filepath_temperature_e=input_path / 'input.plasma_1d' 
    filepath_density_e=input_path / 'input.plasma_1d' 
    filepath_beam_deposition=input_path / 'input.particles'
    filepath_wall=input_path / 'input.wall_2d'

    temperature_array=[] #arrays to hold kinetic profile data for each species
    density_array=[]
    for species_number in species_numbers: #cycle through species and read corresponding kinetic profiles

        try:
            temperature=classes.input_classes.temperature.Temperature(ID='made using read_inputs_ASCOT()',data_format='ASCOT',filename=filepath_temperature_i,species_number=species_number,species='ions')
            temperature_array.append(temperature)
        except:
            print("WARNING: read_inputs_ASCOT() could not read ion temperature for species {species_number} from LOCUST_IO/input_files/{input_path}".format(species_number=species_number,input_path=filepath_temperature_i))
        try:
            density=classes.input_classes.number_density.Number_Density(ID='made using read_inputs_ASCOT()',data_format='ASCOT',filename=filepath_density_i,species_number=species_number,species='ions')
            density_array.append(density)
        except:
            print("WARNING: read_inputs_ASCOT() could not read ion number density for species {species_number} from LOCUST_IO/input_files/{input_path}".format(species_number=species_number,input_path=filepath_density_i))

    try:
        temperature_e=classes.input_classes.temperature.Temperature(ID='made using read_inputs_ASCOT()',data_format='ASCOT',filename=filepath_temperature_e,species_number=species_number,species='electrons')
    except:
        print("WARNING: read_inputs_ASCOT() could not read electron temperature from LOCUST_IO/input_files/{}".format(filepath_temperature_e))

    try:
        density_e=classes.input_classes.number_density.Number_Density(ID='made using read_inputs_ASCOT()',data_format='ASCOT',filename=filepath_density_e,species_number=species_number,species='electrons')
    except:
        print("WARNING: read_inputs_ASCOT() could not read electron number density from LOCUST_IO/input_files/{}".format(filepath_density_i))

    if beam_depo_GC:
        data_format_beam_depo='ASCOT_GC'
    else: 
        data_format_beam_depo='ASCOT_FO'
    try:
        beam_deposition=classes.input_classes.beam_deposition.Beam_Deposition(ID='made using read_inputs_ASCOT()',data_format=data_format_beam_depo,filename=filepath_beam_deposition)
    except:
        print("WARNING: read_inputs_ASCOT() could not read beam deposition from LOCUST_IO/input_files/{} - returning None".format(filepath_beam_deposition))
        beam_deposition=None

    if wall_type=='2D' or wall_type=='3D':
        wall_type+='_input'
    else:
        wall_type='2D_input'
        print("read_inputs_ASCOT() assuming 2D wall")
    data_format_wall='ASCOT_'+wall_type

    try:
        wall=classes.input_classes.wall.Wall(ID='made using read_inputs_ASCOT()',data_format=data_format_wall,filename=filepath_wall)
    except:
        print("WARNING: read_inputs_ASCOT() could not read wall from LOCUST_IO/input_files/{} - returning None".format(filepath_wall))
        wall=None

    print("finished read_inputs_ASCOT()")

    return temperature_array,density_array,temperature_e,density_e,beam_deposition,wall

def dump_inputs_ASCOT(temperature_i,temperature_e,density_i,density_e,rotation_toroidal,equilibrium,beam_deposition,wall,beam_depo_GC=True,input_path=pathlib.Path(''),tag=''):
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
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        input_path - path to target in input_files dir (input_files/path/)
        tag - optional identifier tag for each set of run files produced
    """

    print("dump_inputs_ASCOT creating ASCOT inputs")

    try:
        dump_run_file_ASCOT(input_path=input_path,tag=tag) #generate run file
    except:
        print("WARNING: dump_inputs_ASCOT() could not dump_run_file_ASCOT to LOCUST_IO/input_files/{}".format(input_path))

    try:
        dump_profiles_ASCOT(filename=input_path / str('input.plasma_1d'+tag),temperature_i=temperature_i,temperature_e=temperature_e,density_i=density_i,density_e=density_e,rotation_toroidal=rotation_toroidal)
    except:
        print("WARNING: dump_inputs_ASCOT() could not dump_profiles_ASCOT to LOCUST_IO/input_files/{}".format(input_path / str('input.plasma_1d'+tag)))

    try:
        dump_input_options_ASCOT(filename=input_path / str('input.options'+tag))
    except:
        print("WARNING: dump_inputs_ASCOT() could not dump_input_options_ASCOT to LOCUST_IO/input_files/{}".format(input_path / str('input.options'+tag)))

    try:
        wall.dump_data(data_format='ASCOT_2D_input',filename=input_path / str('input.wall_2d'+tag))
    except:
        print("WARNING: dump_inputs_ASCOT() could not dump ion temperature to LOCUST_IO/input_files/{}".format(input_path / str('input.wall_2d'+tag)))

    try:
        if beam_depo_GC:
            beam_deposition.dump_data(data_format='ASCOT_GC',filename=input_path / str('input.particles'+tag),equilibrium=equilibrium)
        else:
            beam_deposition.dump_data(data_format='ASCOT_FO',filename=input_path / str('input.particles'+tag),equilibrium=equilibrium)
    except:
        print("WARNING: dump_inputs_ASCOT() could not dump beam_deposition to LOCUST_IO/input_files/{}".format(input_path / str('input.particles'+tag)))

    print("dump_inputs_ASCOT finished creating ASCOT inputs")

################################################################################################### MARSF classes and functions

def dump_rotation_MARSF(filename,some_rotation):
    """
    writes rotation profile to MARSF Mogui ASCII format 

    notes
        assumes structure similar number_density or temperature, with normalised poloidal flux given against rotation
        re-usable if rotation class written for LOCUST_IO
        writes out a header line for number of points
        MARSF mogui written by David Ryan
    args:
        filename - output filename
        some_rotation - data structure holding 'rotation' against normalised poloidal flux
    """

    print("writing rotation to MARSF mogui")

    filepath=support.dir_input_files / filename

    with open(filepath,'w') as file: #open file

        flux_pol_norm_sqrt=np.sqrt(np.abs(some_rotation['flux_pol_norm'])) #calculate profiles vs np.sqrt(flux_pol)
        flux_pol_norm_sqrt,rotation=processing.utils.sort_arrays(flux_pol_norm_sqrt,some_rotation['rotation']) #check order
 
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

    filepath=support.dir_input_files / filename

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
        self.filepath=support.dir_output_files / filename
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

        print("reading LOCUST FINT distribution")

        with open(filepath,'r') as file:

            lines=file.readlines() #reading entire file anyway so grab all at once
            if not lines: #check to see if the file opened
                raise IOError("ERROR: read_beam_depo_LOCUST() cannot read from "+str(filepath))


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
            self['PFC_power']=np.asarray(self['PFC_power'])*1.e6 #convert from [MW] to [W]
            self['dfn']=np.asarray(self['dfn'])
            self['dfn']=self['dfn'].reshape(int(len(self['time'])),int(self['nE'])) #XXX check order of this

        print("finished reading LOCUST FINT distribution")

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
        ax.set_title(self.ID)

        E,time=np.meshgrid(self['E'],self['time'])

        if fill:
            ax.set_facecolor(colmap(np.amin(self['dfn'])))
            mesh=ax.pcolormesh(time,E,self['dfn'],cmap=colmap,vmin=np.amin(self['dfn']),vmax=np.amax(self['dfn']))
        else:
            mesh=ax.contour(time,E,self['dfn'],levels=np.linspace(np.amin(self['dfn']),np.amax(self['dfn']),num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(self['dfn']),vmax=np.amax(self['dfn']))
            #ax.clabel(mesh,inline=1,fontsize=10)

        ax.set_xlabel('time [s]')
        ax.set_ylabel('energy [eV]')                
        
        if fig_flag is False:    
            fig.colorbar(mesh,ax=ax,orientation='horizontal')

        if ax_flag is True or fig_flag is True: #return the plot object
            return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()


def dump_perturbation_point_data_input(BCHECK=1,**kwargs):
    """
    generates the point_data.inp file for checking magnetic perturbations using LOCUST -DBCHECK

    args:
        BCHECK - coordinate format setting for LOCUST field checking (1=RPhiZ,2=XYZ)  
    notes:
    usage:
       dump_perturbation_point_data_input(R=[1],phi=[2],Z=[3],time=[0]) 
       dump_perturbation_point_data_input(BCHECK=2,X=[1],Y=[2],Z=[3],time=[0])
    """

    print("writing point_inp.dat test points")

    filepath=support.dir_input_files / 'point_data.inp'
 
    with open(filepath,'w') as file: #open file

        if BCHECK==1:
            for R,Phi,Z,time in zip(kwargs['R'],kwargs['phi'],kwargs['Z'],kwargs['time']):
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
            for X,Y,Z,time in zip(kwargs['X'],kwargs['Y'],kwargs['Z'],kwargs['time']):
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

def dump_inputs_LOCUST(temperature_i=None,temperature_e=None,density_e=None,equilibrium=None,beam_deposition=None,wall=None,perturbation=None,beam_depo_GC=False,beam_depo_weighted=True,BCHECK=False,wall_type='2D',input_path=pathlib.Path(''),tag=''):
    """
    generates full run input data for LOCUST

    notes:
        dumps everything into input_files folder
        uses default LOCUST filenames
        some inputs from list below are missing and must be generated by hand for now
        LOCUST run requires:
            - collisions.dat = cross-section data (needs to be retrieved externally)
            - ptcles.dat = particle birth list
            - profile_ne.dat = electron density profile
            - profile_Te.dat = electron temperature profile
            - profile_Ti.dat = ion temperature profile
            - GEQDSK = B field
            - wall = 2D (supported) or 3D wall (not supported) - optional
            - point_data.inp = B field sample points for use with -DBCHECK mode - optional
            - dB_map.dat = perturbation input - optional
            - BPLASMA_n = perturbation input - optional
    args:
        temperature_i - ion temperature object (eV)
        temperature_e - electron temperature object (eV)
        density_e - electron density object (#/m^3)
        equilibrium - equilibrium object 
        beam_deposition - beam deposition object
        wall - wall object
        perturbation - perturbation object
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        beam_depo_weighted - toggle dumping weighted birth list 
        BCHECK - toggle dumping point_data.inp file for use with LOCUST -DBCHECK mode
        wall_type - set wall type to '2D' or '3D'
        input_path - path to target in input_files dir (input_files/path/)
        tag - optional identifier tag for each set of run files produced
    """

    print("dump_inputs_LOCUST creating LOCUST inputs")

    filepath_temperature_i=input_path / str('profile_Ti.dat'+tag)
    filepath_temperature_e=input_path / str('profile_Te.dat'+tag)
    filepath_number_density_e=input_path / str('profile_ne.dat'+tag)
    filepath_equilibrium=input_path / str('LOCUST_GEQDSK'+tag)
    filepath_beam_deposition=input_path / str('ptcles.dat'+tag)
    filepath_wall=input_path / str('LOCUST_wall'+tag)
    filepath_perturbation=input_path / str('pert'+tag)
    filepath_point_data=input_path / str('point_data.inp'+tag)

    if temperature_i:
        try:
            temperature_i.dump_data(data_format='LOCUST',filename=filepath_temperature_i)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump ion temperature to LOCUST_IO/input_files/{}".format(filepath_temperature_i))

    if temperature_e:
        try:
            temperature_e.dump_data(data_format='LOCUST',filename=filepath_temperature_e)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump electron temperature to LOCUST_IO/input_files/{}".format(filepath_temperature_e))

    if density_e:
        try:
            density_e.dump_data(data_format='LOCUST',filename=filepath_number_density_e)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump number density to LOCUST_IO/input_files/{}".format(filepath_number_density_e))

    if equilibrium:
        try:
            equilibrium.dump_data(data_format='GEQDSK',filename=filepath_equilibrium)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump equilibrium to LOCUST_IO/input_files/{}".format(filepath_equilibrium))

    if beam_deposition:
        try:
            if beam_depo_GC:
                if beam_depo_weighted:
                    beam_deposition.dump_data(data_format='LOCUST_GC_weighted',filename=filepath_beam_deposition,equilibrium=equilibrium)
                else:
                    beam_deposition.dump_data(data_format='LOCUST_GC',filename=filepath_beam_deposition,equilibrium=equilibrium) #XXX not implemented yet            
            else:
                if beam_depo_weighted:
                    beam_deposition.dump_data(data_format='LOCUST_FO_weighted',filename=filepath_beam_deposition)
                else:
                    beam_deposition.dump_data(data_format='LOCUST_FO',filename=filepath_beam_deposition)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump beam deposition to LOCUST_IO/input_files/{}".format(filepath_beam_deposition))

    if wall:
        try:
            data_format_wall='LOCUST_'+wall_type
            wall.dump_data(data_format=data_format_wall,filename=filepath_wall)
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump wall to LOCUST_IO/input_files/{}".format(filepath_wall))

    if perturbation:
        try:
            perturbation.dump_data(data_format='LOCUST',filename=filepath_perturbation) #XXX not implemented yet
            if BCHECK:
                try:
                    perturbation.dump_data(data_format='point_data',filename=filepath_point_data)
                except:
                    print("WARNING: dump_inputs_LOCUST() could not dump point_data.inp to LOCUST_IO/input_files/{}".format(filepath_point_data))                    
        except:
            print("WARNING: dump_inputs_LOCUST() could not dump perturbation to LOCUST_IO/input_files/{}".format(filepath_perturbation))

    print("dump_inputs_LOCUST finished")


def calc_coulomb_logarithm(Et,n,T,Ai,At,Z,Zt,Bmod):
    """
    calculates the Coulomb logarithm at single test particle energy Et for arbitrary background species

    notes: 
        must pass all species information to this routine, since Rmax is sum over species
        returns lnL, an array with a Coulomb logarithm corresponding to each species
    args:
        Et - test particle energy [eV]
        n - array of densities for background species [m**-3]
        T - array of temperatures for background species [eV]
        Ai - array of atomic masses for background species [amu]
        At - atomic mass for test particle species [amu]
        Z - array of atomic charges for background species [int]
        Zt - atomic charge for test particle species [int]
        Bmod - magnetic field strength [T]
    """
    Et/=1000. #convert to keV
    T/=1000. #convert to KeV

    n=np.asarray(n)
    T=np.asarray(T)
    Ai=np.asarray(Ai)
    Z=np.asarray(Z)

    omega2 = 1.74*Z**2/Ai*n + 9.18e15*Z**2/Ai**2*Bmod**2
    vrel2  = 9.58e10*(T/Ai + 2.0*Et/At)
    rmx    = np.sqrt(1.0/np.sum(omega2/vrel2))

    vrel2  = 9.58e10*(3.0*T/Ai+2.0*Et/At)
    rmincl = 0.13793*np.abs(Z*Zt)*(Ai+At)/Ai/At/vrel2
    rminqu = 1.9121e-08*(Ai+At)/Ai/At/np.sqrt(vrel2)

    rmn=[]
    for rmincl_,rminqu_ in zip(rmincl,rminqu): 
        rmn.extend([np.max([rmincl_,rminqu_])])
    rmn=np.asarray(rmn)

    lnL=np.log(rmx/rmn)

    return lnL

def calc_chandrasekhar_function(x,mi,mt):
    """
    calculate Chandrasekhar function G(x)=(ERF(x) - x.d{ERF(x)}/dx)/2x**2

    notes:
        N.B. ERF(x) using function erf_7_1_26 breaks down at small x.
        Routine resorts to small x analytic approximation for x<0.05.
        Routine is only valid for x>=0.0
        Handbook of Mathematical Functions. Edited by Milton Abramowitz and 
        Irene A. Stegun, Dover Publications, Inc., New York, 1965.
        Error function and Fresnel Integrals, EQN. 7.1.26.
        Valid to |E(x)| <= 1.5e-7. Calculation in gpu precision
        f1 - = ERF(x) - (1+mi/Mb).x.d{ERF(x)}/dx
        f2 - = ERF(x) - G(x)
    args:
        x - X
        mi - background particle mass [kg]
        mt - test particle mass [kg]
    """

    a1 = 0.254829592
    a2 =-0.284496736
    a3 = 1.421413741
    a4 =-1.453152027
    a5 = 1.061405429
    p  = 0.327591100
          
    i  = np.where(x<0.0)[0]
    if len(i)>0:
        print("ERROR: calc_chandrasekhar_function invalid for x<0!\nreturning!\n")
        return

    t     = 1.0/(1.0+p*x)          
    exp_  = np.exp(-x**2)
    erf_s = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp_
    G     = ( erf_s - x*2.0*exp_/np.sqrt(constants.pi) )/( 2.0*x**2)

    f1    = erf_s - (1.0+(mi/mt))*x*2.0*exp_/np.sqrt(constants.pi)
    f2    = erf_s - G

    i     = np.where( x<0.005)[0]
    if len(i)!=0:     
        #eqn. for G breaks down at low x - revert to analytic approximation:
        x     = x[i]
        G[i]  = 2.0*x/(3.0*np.sqrt(constants.pi))
        f1[i] = (2.0/np.sqrt(constants.pi))*( ( (2.0/3.0) + (mi/mt) )*x**3 - x*(mi/mt) ) 
        f2[i] = (2.0/(3.0*np.sqrt(constants.pi)))*(2.0-x**2)*x

    return G,f1,f2

################################################################################################### misc functions



#################################
 
##################################################################
 
###################################################################################################
