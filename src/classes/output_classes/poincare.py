#poincare.py

'''
Samuel Ward
23/08/2020
----
class to handle LOCUST output poincare maps
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

################################################################## Poincare functions

def read_poincare_LOCUST(filepath,**properties):
    """
    notes:
        assumes niter=1
    """

    try:
        from scipy.io import FortranFile 
    except:
        raise ImportError("ERROR: read_poincare_LOCUST could not import scipy.io.FortranFile!\nreturning\n")
        return

    file=FortranFile(filepath,'r')

    input_data={}

    input_data['EQBM_md5']=file.read_reals(dtype=np.float32)

    input_data['nRP'], input_data['nZP'], input_data['nTP'], input_data['npln'] = file.read_ints()

    for _ in range(6):
        file.read_record(dtype=np.byte)
    input_data['phi']=np.linspace(0,2.*np.pi,input_data['nTP'])/input_data['npln']

    input_data['phi']=np.linspace(0,2.*np.pi,input_data['nTP'])/input_data['npln']

    input_data['R0P'],input_data['R1P'],input_data['Z0P'],input_data['Z1P']=file.read_reals(dtype=np.float32)
    
    input_data['map']=np.empty(shape=(input_data['nTP'],input_data['nZP'],input_data['nRP']))
    for i in range(input_data['nTP']):
        input_data['map'][i,:,:]=file.read_ints().reshape(input_data['nZP'],input_data['nRP'])
    input_data['map']=np.array(input_data['map'],dtype=float).T

    input_data['nEQ_R']=file.read_ints()
    input_data['nEQ_Z']=file.read_ints()
    input_data['R_1D']=file.read_reals(dtype=np.float32)
    input_data['Z_1D']=file.read_reals(dtype=np.float32)
    input_data['psi_equil_h']=file.read_reals(dtype=np.float32)
    
    #calculate extra info
    input_data['R']=np.linspace(input_data['R0P'],input_data['R1P'],input_data['nRP'])
    input_data['Z']=np.linspace(input_data['Z0P'],input_data['Z1P'],input_data['nZP'])

    return input_data

################################################################## Poincare class

class Poincare(classes.base_output.LOCUST_output):
    """
    class describing a poincare map output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'poincare'
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

    LOCUST_output_type='poincare'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read poincare from file 

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
                self.data=read_poincare_LOCUST(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))            

    def plot(self,phi=0,style='histogram',LCFS=False,limiters=False,real_scale=True,colmap=settings.cmap_default,label='',ax=False,fig=False):
        """
        plot the poincare map
         
        args:
            phi - toroidal angle to slice at (nearest neighbour) [rad]
            style - choose from scatter or histogram
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - plot to Tokamak scale (requires equilibrium argument)
            colmap - set the colour map (use get_cmap names)
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        phi=run_scripts.utils.literal_eval(phi)

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
            ax = fig.add_subplot(111)

        ax.set_title(self.ID)
        X=self['R'] #make a mesh
        Y=self['Z']
        ax.set_xticks(X) #set axes ticks
        ax.set_yticks(Y)
        dx,dy=X[1]-X[0],Y[1]-Y[0]
        Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres
        phi_slice=np.argmin(np.abs(self['phi']-phi))

        mesh=None

        if style is 'histogram':
            mesh=ax.pcolormesh(X,Y,self['map'][:,:,phi_slice],cmap=colmap,edgecolor='none',antialiased=True)
        elif style is 'scatter':
            inds=np.where(self['map'][:,:,phi_slice].flatten()==1)[0]
            ax.scatter(X.flatten()[inds],Y.flatten()[inds],cmap=colmap,edgecolor='none',antialiased=True,s=0.03)

        if LCFS:
            ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
        if limiters: #add boundaries if desired
            ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall') 

        ax.set_xlim(np.min(self['R']),np.max(self['R']))
        ax.set_ylim(np.min(self['Z']),np.max(self['Z']))
        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')

        if ax_flag is True or fig_flag is True: #return the plot object
            return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()

    def animate(self,save=False,*args):
        """
        plot the poincare map
         
        args:
            save - toggle whether to save gif
            *args - positional args to pass to .plot()
        """

        fig = plt.figure() #if user has not externally supplied figure, generate
        phi=np.linspace(0.,2.*np.pi)
        animation=FuncAnimation(fig,self.plot,frames=phi,fargs=[*args],repeat=True,interval=1)
        plt.show()
        if save: animation.save('poincare_animation.gif',writer='pillow')

#################################

##################################################################

###################################################################################################
