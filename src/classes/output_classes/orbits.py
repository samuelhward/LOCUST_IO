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

################################################################## Orbit read functions

def read_orbits_LOCUST(filepath,**properties):
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
            raise IOError("ERROR: read_orbits_LOCUST() cannot read from "+str(filepath))
    
        number_particles=int(lines[0]) #extract number of particles
        number_timesteps=int(lines[-1])-1 #extract number of time steps of each trajectory
        number_coords=int(3)

        del(lines[0])
        del(lines[-1])

        input_data = {} #initialise the dictionary to hold the dat
        orbit_data=np.array([[float(value) for value in line.split()] for line in lines]).reshape((number_timesteps+1,number_coords,number_particles)) #read all data and reshape accordingly
        input_data['R']=orbit_data[:,0,:]
        input_data['phi']=orbit_data[:,1,:]
        input_data['Z']=orbit_data[:,2,:]
        input_data['number_particles']=np.asarray(number_particles)
        input_data['number_timesteps']=np.asarray(number_timesteps)
   
    print("finished reading orbits from LOCUST")

    return input_data

def read_orbits_ASCOT(filepath,**properties):
    """
    reads orbits stored in ASCOT

    notes:
        reads only one particle
        optionally read the particle with ID stored in properties['ID'] 
    """

    print("reading orbits from ASCOT")

    with open(filepath) as file:

        if 'ID' in properties:
            ID_desired=properties['ID']
        else:
            ID_desired=0 #by default read particle 0

        undesired_variables=['rho']

        input_data={} #to store the desired data that is output
        data_key={} #to tell us which column corresponds to which data field

        for line in file:
            if '# Number of different fields for each particle [10 first letters are significant]' in line:
                number_of_data_fields=int(line.split()[0])
                break #now at the point where the data fields are defined

        for data_field in range(number_of_data_fields):
            variable=file.readline().split()[0]
            if variable=='pitch':
                variable='V_pitch'
            if variable=='z':
                variable='Z'
            if variable=='energy':
                variable='E'
            if variable=='id':
                variable='ID'
            input_data[variable]=[]
            data_key[variable]=int(data_field)

        for variable in undesired_variables: #get rid of data we do not want to read
            del(data_key[variable])

        line=file.readline() #read a blank line 

        for line in file:
            split_line=line.split()
            if len(split_line)<number_of_data_fields: #reached end of file
                break
            elif int(split_line[data_key['ID']])!=ID_desired: #check if this is the particle we care about
                pass
            else:
               for variable_name,field_number in data_key.items():
                    input_data[variable_name].append(float(split_line[field_number]))

        for variable in input_data:
            input_data[variable]=np.asarray(input_data[variable])[:,np.newaxis] #add fake axis to match LOCUST format

    print("finished reading orbits from ASCOT")

    return input_data

################################################################## Orbit write functions

def dump_orbits_LOCUST(output_data,filepath,**properties): 
    """
    writes orbits to LOCUST format - r phi z
    
    notes:
        writes out a headerline for number of particles
        writes out a footerline for number of time steps
    """

    print("ERROR: dump_orbits_LOCUST not yet implemented!\nreturning\n")
    
    '''
    print("writing orbits to LOCUST")
    with open(filepath,'w') as file: #open file

        file.write("{}\n".format(output_data['number_particles'].size)) #re-insert line containing number of particles


        #some stuff here


        file.write("{}".format(output_data['number_timesteps'].size)) #re-insert line containing number of time steps

    print("finished writing orbits to LOCUST")
    '''

'''
def dump_orbits_vtk(output_data,filepath,**properties)
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
        data is stored such that coordinate coord at time t for particle p is my_orbit[coord][t,p]
    """

    LOCUST_output_type='orbits'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read orbits from file 

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
                self.data=read_orbits_LOCUST(self.filepath,**properties)

        elif data_format=='ASCOT': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from ASCOT - filename required\n".format(self.ID),filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_orbits_ASCOT(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/ASCOT)\n".format(self.ID))            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write orbits to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_output_files / filename
                dump_orbits_LOCUST(self.data,filepath,**properties)
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))

    def plot(self,particles=[0],axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,start_mark=False,colmap=settings.cmap_default,colmap_val=np.random.uniform(),ax=False,fig=False):
        """
        simple orbits plot in the R,Z/X,Y planes
         
        args:
            particles - iterable list of particle numbers
            axes - define plot axes
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - plot to Tokamak scale (requires equilibrium argument)
            start_mark - include marker to show birth point
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,particles,colmap_val=run_scripts.utils.literal_eval(axes,particles,colmap_val)

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
        
        ndim=len(axes)
        if fig_flag is False:
            fig = plt.figure() #if user has not externally supplied figure, generate
        
        if ax_flag is False: #if user has not externally supplied axes, generate them
            ax = fig.add_subplot(111)
        ax.set_title(self.ID)

        if ndim==2: #2D plotting

            for particle in particles:
                ax.plot(self[axes[0]][:,particle],self[axes[1]][:,particle],color=colmap(colmap_val),linewidth=settings.plot_linewidth)
                if start_mark:
                    ax.plot(self[axes[0]][0,particle],self[axes[1]][0,particle],color=settings.colour_start_mark,marker=settings.marker_start_mark,markersize=settings.markersize_start_mark)

            if axes==['R','Z']: #if we're plotting along poloidal projection, then give options to include full cross-section and plasma boundary
               
                if real_scale is True:
                    ax.set_aspect('equal')

                if LCFS: #plot plasma boundary
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],settings.plot_style_LCFS) 
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],settings.plot_style_limiters)

            elif axes==['X','Y']: #plotting top-down
                
                if real_scale is True:
                    ax.set_aspect('equal')

                if LCFS: #plot concentric rings to show inboard/outboard plasma boundaries
                    plasma_max_R=np.max(LCFS['lcfs_r'])
                    plasma_min_R=np.min(LCFS['lcfs_r'])
                    ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),settings.plot_style_LCFS)
                    ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),settings.plot_style_LCFS)
                if limiters: #add boundaries if desired
                    ax.set_xlim(-1.0*np.max(limiters['rlim']),np.max(limiters['rlim']))
                    ax.set_ylim(-1.0*np.max(limiters['rlim']),np.max(limiters['rlim']))
                    limiters_max_R=np.max(limiters['rlim'])
                    limiters_min_R=np.min(limiters['rlim'])
                    ax.plot(limiters_max_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_max_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),settings.plot_style_limiters)
                    ax.plot(limiters_min_R*np.cos(np.linspace(0,2.0*constants.pi,100)),limiters_min_R*np.sin(np.linspace(0.0,2.0*constants.pi,100)),settings.plot_style_limiters)           

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])

        elif ndim==3: #3D plotting

            if ax_flag is False:
                ax = fig.gca(projection='3d')

            if real_scale:

                ax.set_aspect('equal')

            if LCFS: #plot periodic poloidal cross-sections in 3D 
                for angle in np.linspace(0.0,2.0*constants.pi,4,endpoint=False):
                    x_points=LCFS['lcfs_r']*np.cos(angle)
                    y_points=LCFS['lcfs_r']*np.sin(angle)
                    z_points=LCFS['lcfs_z']
                    ax.plot(x_points,y_points,zs=z_points,color=settings.plot_style_LCFS)
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],settings.plot_style_limiters)

            for particle in particles:
                ax.plot(self[axes[0]][:,particle],self[axes[1]][:,particle],zs=self[axes[2]][:,particle],color=colmap(colmap_val),linewidth=settings.plot_linewidth)
                if start_mark:
                    ax.plot(self[axes[0]][0,particle],self[axes[1]][0,particle],zs=self[axes[2]][0,particle],color=settings.colour_start_mark,marker=settings.marker_start_mark,s=settings.markersize_start_mark)

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
            ax.set_zlabel(axes[2])

        if ax_flag is False and fig_flag is False:
            plt.show()

#################################

##################################################################

###################################################################################################
