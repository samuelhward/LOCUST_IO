#orbits.py

"""
Samuel Ward
15/01/2018
----
class to handle LOCUST orbit output data
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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_output 
except:
    raise ImportError("ERROR: LOCUST_IO/classes/base_output.py could not be imported!\nreturning\n")
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

################################################################## Orbit read functions

def read_orbits_LOCUST(filepath):
    """
    reads orbits stored in LOCUST format - r phi z

    notes:
        reads in a headerline for number of particles
        reads in a footerline for number of time steps
    """

    print("reading orbits from LOCUST")

    with open(filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_orbits_LOCUST() cannot read from "+filepath)
    
        number_particles=int(lines[0]) #extract number of particles
        number_timesteps=int(lines[-1])-1 #extract number of time steps of each trajectory
        number_coords=int(3)

        del(lines[0])
        del(lines[-1])

        input_data = {} #initialise the dictionary to hold the dat
        input_data['orbits']=np.array([[float(value) for value in line.split()] for line in lines]).reshape((number_timesteps+1,number_coords,number_particles)) #read all data and reshape accordingly
        input_data['number_particles']=np.asarray(number_particles)
        input_data['number_timesteps']=np.asarray(number_timesteps)
   
    print("finished reading orbits from LOCUST")

    return input_data

################################################################## Orbit write functions

'''
def dump_orbits_LOCUST(output_data,filepath): 
    """
    writes orbits to LOCUST format - r phi z
    
    notes:
        writes out a headerline for number of particles
        writes out a footerline for number of time steps
    """

    print("writing orbits to LOCUST")

    with open(filepath,'w') as file: #open file

        file.write("{}\n".format(output_data['number_particles'].size)) #re-insert line containing number of particles


        #some stuff here


        file.write("{}".format(output_data['number_timesteps'].size)) #re-insert line containing number of time steps

    print("finished writing orbits to LOCUST")
'''
'''
def dump_orbits_vtk(output_data,filepath)
    """
    notes:
    
    here just need to print x,y,z \n locations at all times for particle 1 THEN all times for particle 2 etc
    """

    printF, 1, '# vtk DataFile Version 2.0'
    printF, 1, 'Unstructured Grid ORBIT data'
    printF, 1, 'ASCII'
    printF, 1, 'DATASET UNSTRUCTURED_GRID'
    printF, 1, 'POINTS '+strcompress(string(npoints),/re)+' float'
'''

################################################################## Orbits class

class Orbits(classes.base_output.LOCUST_output):
    """
    class describing orbits output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'orbits'
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
        data is stored such that coordinate i at time t for particle p is my_orbit['orbits'][t,i,p]
        in this way, a single particle's trajectory is my_orbit['orbits'][:,i,p] where i=0,1,2=r,phi,z
    """

    LOCUST_output_type='orbits'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read orbits from file 

        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties={**properties}
                self.data=read_orbits_LOCUST(self.filepath) #read the file
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write orbits to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_orbits_LOCUST(self.data,filepath)
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST)\n")

    def plot(self,some_equilibrium=False,particles=[0],axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,start_mark=False,colmap='blue',ax=False,fig=False):
        """
        simple orbits plot in the R,Z/X,Y planes
         
        args:
            some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
            particles - iterable list of particle numbers
            axes - define plot axes
            LCFS - show plasma boundary outline (requires equilibrium arguement)
            limiters - toggles limiters on/off in 2D plots
            real_scale - plot to Tokamak scale (requires equilibrium arguement)
            start_mark - include marker to show birth point
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
        
        ndim=len(axes)
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)

        if ndim==2: #2D plotting
            
            if axes==['R','Z']: #if we're plotting along poloidal projection, then give options to include full cross-section and plasma boundary
               
                if real_scale is True:
                    if some_equilibrium:
                        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                    ax.set_aspect('equal')

                if LCFS is True: #plot plasma boundary
                    ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],plot_style_LCFS) 
                if limiters is True: #add boundaries if desired
                    ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)

                for particle in particles: #plot all the particle tracks one by one
                    ax.plot(self['orbits'][1::2,0,particle],self['orbits'][1::2,2,particle],color=colmap,linewidth=0.5) #plot every other position along trajectory
                    if start_mark: #show birth point
                        ax.plot(self['orbits'][0,0,particle],self['orbits'][0,2,particle],color=colour_start_mark,marker='o',markersize=1)

            elif axes==['X','Y']: #plotting top-down
                
                if real_scale is True:
                    if some_equilibrium:
                        ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                        ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_aspect('equal')

                if LCFS is True: #plot concentric rings to show inboard/outboard plasma boundaries
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

                for particle in particles: #plot all the particle tracks one by one
                    x_points=self['orbits'][1::2,0,particle]*np.cos(self['orbits'][1::2,1,particle]) #calculate using every other position along trajectory
                    y_points=self['orbits'][1::2,0,particle]*np.sin(self['orbits'][1::2,1,particle])   
                    ax.plot(x_points,y_points,color=colmap,linewidth=1) 
                    
                    if start_mark: #show birth point
                        ax.plot(x_points[0],y_points[0],color=colour_start_mark,marker='o',markersize=1)

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
            ax.set_title(self.ID)


        elif ndim==3: #3D plotting

            if ax_flag is False:
                ax = fig.gca(projection='3d')

            if real_scale:
                if some_equilibrium:
                    ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')

            if LCFS is True: #plot periodic poloidal cross-sections in 3D 
                for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                    x_points=some_equilibrium['lcfs_r']*np.cos(angle)
                    y_points=some_equilibrium['lcfs_r']*np.sin(angle)
                    z_points=some_equilibrium['lcfs_z']
                    ax.plot(x_points,y_points,zs=z_points,color='m')
                if limiters is True: #add boundaries if desired
                    ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)

            for particle in particles: #plot all the particle tracks one by one
                x_points=self['orbits'][1::2,0,particle]*np.cos(self['orbits'][1::2,1,particle]) #calculate using every other position along trajectory
                y_points=self['orbits'][1::2,0,particle]*np.sin(self['orbits'][1::2,1,particle])   
                z_points=self['orbits'][1::2,2,particle]
                ax.plot(x_points,y_points,zs=z_points,color=colmap) 
                
                if start_mark: #show birth point
                    mesh=ax.scatter(x_points[0],y_points[0],z_points[0],color=colour_start_mark,marker='o',s=10,label=self.ID)

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
            ax.set_zlabel(axes[2])
            ax.set_title(self.ID)

        if ax_flag is False and fig_flag is False:
            plt.show()

#################################

##################################################################

###################################################################################################