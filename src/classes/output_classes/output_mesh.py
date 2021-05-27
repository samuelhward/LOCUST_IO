#output_mesh.py

'''
Samuel Ward
23/08/2020
----
class to handle LOCUST output PFC mesh 
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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_output 
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/base_output.py could not be imported!\nreturning\n")
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

################################################################## output_mesh functions

def read_output_mesh_LOCUST(filepath,**properties):
    """
    notes:
        assumes niter=1
    """

    try:
        import vtk 
    except:
        raise ImportError("ERROR: read_output_mesh_LOCUST could not import vtk module!\nreturning\n")
        return

    input_data={}

    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(str(filepath))
    reader.Update()
    input_data['X'],input_data['Y'],input_data['Z']=np.array(reader.GetOutput().GetPoints().GetData()).T
    input_data['R'],input_data['phi']=processing.utils.XYZ_to_RphiZ(input_data['X'],input_data['Y'])

    return input_data

################################################################## output_mesh class

class Output_Mesh(classes.base_output.LOCUST_output):
    """
    class describing a output_mesh map output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'output_mesh'
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
    """

    LOCUST_output_type='output_mesh'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read output_mesh from file 

        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_output_mesh_LOCUST(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))            

    def plot(self,style='line',axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,label='',ax=False,fig=False):
        """
        plots beam deposition

        notes:
            style - choose from scatter or histogram
            axes - list of strings specifying which axes should be plotted
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - sets r,z scale to real tokamak cross section
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            line_style - set 1D line style
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,colmap_val=run_scripts.utils.literal_eval(axes,colmap_val)

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
            
        ax.set_title(self.ID)

        ndim=len(axes)
        if ndim==2: #plot 2D histograms

            if style=='line':
                ax.plot(self[axes[0]],self[axes[1]],color=colmap(colmap_val),label=self.ID)

            elif style=='scatter':
                ax.scatter(self[axes[0]],self[axes[1]],color=colmap(colmap_val),marker='x',s=1,label=self.ID)

            if axes==['R','Z']:
                if LCFS: #plot plasma boundary
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')

            elif axes==['X','Y']:

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
            
            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])

            if real_scale is True: #set x and y plot limits to real scales
                ax.set_aspect('equal')
           
        elif ndim==3: #plot 3D scatter - assume X,Y,Z

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
 
            if style=='line':
                ax.plot(self[axes[0]],self[axes[1]],self[axes[2]],color=colmap(colmap_val),label=self.ID)

            elif style=='scatter':
                ax.scatter(self[axes[0]],self[axes[1]],zs=self[axes[2]],color=colmap(colmap_val),s=0.1,label=self.ID)
        
        if ax_flag is False and fig_flag is False:
            plt.show() 

#################################

##################################################################

###################################################################################################
