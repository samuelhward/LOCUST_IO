#particle_list.py

'''
Samuel Ward
01/04/2018
----
class to handle LOCUST particle list output data
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
    import processing.process_output
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/process_output.py could not be imported!\nreturning\n") 
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

################################################################## Final_Particle_List read functions

def read_final_particle_list_LOCUST(filepath,**properties):
    """
    reads final particle list stored in LOCUST format

    notes:
        filename typically ptcl_cache.dat
        contains lots of references to process_ptcles.pro, written by Rob Akers
        status_flag describes each particle's final status (guide stored in status_flags, verbose guide in LOCUST/ctrk_mod.f90/ctrk_kernel) 
    args:
        compression - toggle True for large files for efficient processing
        coordinates - the particle coordinates to read in if compression enabled        
    """

    print("reading final particle list from LOCUST")

    if 'compression' in properties and properties['compression'] is True: #for now just revert to the efficient compression version
        if 'coordinates' in properties: 
            input_data=processing.process_output.particle_list_compression(filepath,coordinates=properties['coordinates'])
        else:
            input_data=processing.process_output.particle_list_compression(filepath)
        print("finished reading compressed final particle list from LOCUST with compression")
        return input_data

    else: #otherwise we can read the whole file - WARNING this can be huge

        with open(filepath) as file:
            
            try:
                lines=file.readlines() #return lines as list            
            except:
                raise IOError("ERROR: read_final_particle_list_LOCUST() cannot read from "+str(filepath))

            #read in headerlines
            header=lines[0].split()
            n=int(header[0]) 
            ngpu=int(header[1]) #n*ngpu=number of particles
            niter=int(header[2]) #time iterations
            npt_=int(header[3]) #info slots
            nphc=int(header[4]) #levels in -DSPLIT split cache (always use the first)
            ntri=int(header[5]) #triangle grid "dimension" 

            #initialise particle_list and data dictionary
            input_data={}
            input_data['R']=np.array([])
            input_data['phi']=np.array([])
            input_data['Z']=np.array([])
            input_data['V_R']=np.array([])
            input_data['V_phi']=np.array([])
            input_data['V_Z']=np.array([])
            input_data['time']=np.array([])
            input_data['status_flag']=np.array([])
            input_data['additional_flag1']=np.array([])
            input_data['PFC_intercept']=np.array([])
            input_data['psi']=np.array([])
            input_data['V_R_final']=np.array([])
            input_data['V_phi_final']=np.array([])
            input_data['V_Z_final']=np.array([])
            input_data['additional_flag7']=np.array([])
            input_data['additional_flag8']=np.array([])
            input_data['additional_flag9']=np.array([])
            input_data['n']=np.array(n)
            input_data['ngpu']=np.array(ngpu)
            input_data['niter']=np.array(niter)
            input_data['npt_']=np.array(npt_)
            input_data['nphc']=np.array(nphc)
            input_data['ntri']=np.array(ntri)
            input_data['number_particles']=np.array(n*ngpu)

            #get rid of white space and completely flatten IDL/FORTRAN-style
            lines=[[float(number) for number in line.split()] for line in lines]
            lines=[number for line in lines for number in line]
            del(lines[0:6])
            input_data['f']=np.array(lines[-1])
            del(lines[-1])

            for i in range(niter):
                
                #transfer chunk from lines to file_buffer and assimilate into dictionary
                file_buffer=np.array([lines[0:npt_*nphc*(n*ngpu)]]).reshape(npt_,nphc,(n*ngpu),order='F')
                del(lines[0:npt_*nphc*(n*ngpu)])
                
                input_data['R']=np.append(input_data['R'],file_buffer[0,0,:])
                input_data['phi']=np.append(input_data['phi'],file_buffer[1,0,:])
                input_data['Z']=np.append(input_data['Z'],file_buffer[2,0,:])
                input_data['V_R']=np.append(input_data['V_R'],file_buffer[3,0,:])
                input_data['V_phi']=np.append(input_data['V_phi'],file_buffer[4,0,:])
                input_data['V_Z']=np.append(input_data['V_Z'],file_buffer[5,0,:])
                input_data['time']=np.append(input_data['time'],file_buffer[6,0,:])
                input_data['status_flag']=np.append(input_data['status_flag'],file_buffer[7,0,:])
                input_data['additional_flag1']=np.append(input_data['additional_flag1'],file_buffer[8,0,:])
                input_data['PFC_intercept']=np.append(input_data['PFC_intercept'],file_buffer[9,0,:])
                input_data['psi']=np.append(input_data['psi'],file_buffer[10,0,:])
                input_data['V_R_final']=np.append(input_data['V_R_final'],file_buffer[11,0,:])
                input_data['V_phi_final']=np.append(input_data['V_phi_final'],file_buffer[12,0,:])
                input_data['V_Z_final']=np.append(input_data['V_Z_final'],file_buffer[13,0,:])
                input_data['additional_flag7']=np.append(input_data['additional_flag7'],file_buffer[14,0,:])
                input_data['additional_flag8']=np.append(input_data['additional_flag8'],file_buffer[15,0,:])
                input_data['additional_flag9']=np.append(input_data['additional_flag9'],file_buffer[16,0,:])

            input_data['status_flags']={} #nested dictionary to hold possible status flags for the particle list
            input_data['status_flags']['ok_if_greater']=0.0 
            input_data['status_flags']['undefined']=0.0
            input_data['status_flags']['left_space_grid']=-1.0
            input_data['status_flags']['not_poss_on_1st_call']=-1000.0 
            input_data['status_flags']['track_failure']=-2000.0
            input_data['status_flags']['unresolved_hit']=-3.0
            input_data['status_flags']['left_mesh']=-3000.0
            input_data['status_flags']['track_problem']=-4000.0
            input_data['status_flags']['PFC_intercept_2D']=-4.0
            input_data['status_flags']['ptcl_disconnect']=-5000.0
            input_data['status_flags']['PFC_intercept_3D']=-5.0
            input_data['status_flags']['left_field_grid']=-6.0
            input_data['status_flags']['goose_fail']=-7000.0
            input_data['status_flags']['left_plasma']=-8.0
            input_data['status_flags']['thermalised']=-9.0
            input_data['status_flags']['coll_op_fail']=-10000.0
            input_data['status_flags']['GC_calc_fail']=-10.0 
            input_data['status_flags']['CX_loss']=-11.0
            input_data['status_flags']['gc_init_fail']=-11000.0
            input_data['status_flags']['bin_fail_soft']=-12.0
            input_data['status_flags']['bin_fail_hard_1']=-13000.0
            input_data['status_flags']['time_limit_reached']=-14.0
            input_data['status_flags']['cross_open_face']=-15.0
            input_data['status_flags']['bin_fail_hard_2']=-16000.0
            input_data['status_flags']['generic_fail_hard']=-99999.0

            #calculate some additional things
            input_data['E']=.5*constants.species_mass*(input_data['V_R']**2+input_data['V_phi']**2+input_data['V_Z']**2)/constants.species_charge
            input_data['weight']=np.full(len(input_data['R']),1.)

        print("finished reading final particle list from LOCUST")

        return input_data

'''
def read_final_particle_list_TRANSP_FLOST(filepath,**properties):
    """
    reads final particle list stored in TRANSP format

    notes:
    """

    print("reading final particle list from TRANSP")

    with open(filepath,'r') as file:
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_final_particle_list_TRANSP() cannot read from "+str(filepath))

        for counter,line in enumerate(lines): #look for the start of the data, marked by a certain string
            if 'R(cm)' in line: #check for line which has unit headers
                del(lines[0:counter+1])
                break

        input_data = {} #initialise the dictionary to hold the data
        input_data['R']=[] #initialise the arrays
        input_data['Z']=[]
        input_data['V_pitch']=[]
        input_data['energy']=[]    
        input_data['number_particles']=0.0

        for line in lines:
            split_line=line.split()
            split_line[0]=float(split_line[0])
            split_line[1]=float(split_line[1])
            split_line[2]=float(split_line[2])
            split_line[3]=float(split_line[3])
            if split_line: #try ignore whitespace at the end
                input_data['number_particles']+=1

            input_data['R'].append(split_line[0]) 
            input_data['Z'].append(split_line[1])
            input_data['V_pitch'].append(split_line[2])
            input_data['energy'].append(split_line[3])

        input_data['R']=0.01*np.asarray(input_data['R'])
        input_data['Z']=0.01*np.asarray(input_data['Z'])
        input_data['V_pitch']=np.asarray(input_data['V_pitch'])
        input_data['energy']=np.asarray(input_data['energy'])

        input_data['number_particles']=np.array(input_data['number_particles'])
        input_data['status_flags']['any_loss']=0.0
        input_data['status_flag']=np.zeros(len(input_data['R']))
        
    print("finished reading final particle list from TRANSP")

    return input_data
'''

def read_final_particle_list_TRANSP(filepath,**properties):
    """
    read final particle list from TRANSP fi CDF output file

    notes:
        
    """

    print("reading final particle list from TRANSP")

    try:
        from scipy.io import netcdf as ncdf 
    except:
        raise ImportError("ERROR: read_final_particle_list_TRANSP could not import scipy.io.netcdf module!\nreturning\n")
        return

    file=ncdf.netcdf_file(filepath,'r')
    input_data={}

    input_data['R']=file.variables['RMJION'].data*.01
    input_data['Z']=file.variables['ZYION'].data*.01
    input_data['V_pitch']=file.variables['XKSIDLST'].data
    input_data['E']=file.variables['ZELST'].data
    input_data['weight']=file.variables['WGHTLST'].data
    input_data['time']=file.variables['TIMELST'].data

    input_data['status_flags']={} #initialise generic status_flags - TRANSP does store these in ['GOOSELST'] but all are 1.00 or 1.0001 anyway
    input_data['status_flags']['all_losses']=0.0
    input_data['status_flag']=np.zeros(len(input_data['R']))

    file.close()
    file.close()
    
    print("finished reading final particle list from TRANSP")
    
    return input_data

def read_final_particle_list_ASCOT(filepath,**properties):
    """
    read final particle list from ASCOT HDF5 output file

    notes:
        endstates
             1 - maximum tracking time reached
             2 - minimum kinetic energy reached
             3 - collision with wall  
             4 - thermalisation
            -1 - particle was rejected at initialisation, possibly due to starting outside wall or having negative energy etc
            -2 - particle aborted during tracking, possible numerical error etc
    """

    print("reading final particle list from ASCOT")

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_final_particle_list_ASCOT could not import h5py module!\n") 
        return

    with h5py.File(filepath,'r') as file:

        input_data={}

        input_data['R']=file['endstate/R'].value
        input_data['phi']=file['endstate/phi'].value
        input_data['Z']=file['endstate/z'].value
        input_data['V_R']=file['endstate/vR'].value
        input_data['V_phi']=file['endstate/vphi'].value
        input_data['V_Z']=file['endstate/vz'].value
        input_data['E']=file['endstate/energy'].value
        input_data['V_pitch']=file['endstate/pitch'].value
        input_data['time']=file['endstate/time'].value
        input_data['weight']=file['endstate/weight'].value

        input_data['status_flag']=file['endstate/endcond'].value
        input_data['status_flags']={}
        input_data['status_flags']['maximum_time']=1.
        input_data['status_flags']['minimum_kinetic_energy']=2.
        input_data['status_flags']['wall_collision']=3.
        input_data['status_flags']['thermalisation']=4.
        input_data['status_flags']['particle_initial_reject']=-1.
        input_data['status_flags']['particle_abort']=-2.
        input_data['status_flags']['outside_rho_max']=10.
        input_data['status_flags']['outside_wall']=16.
        input_data['status_flags']['cpu_tmax']=17.

    print("finished reading final particle list from ASCOT")

    return input_data

################################################################## Final_Particle_List write functions

def dump_final_particle_list_LOCUST(output_data,filepath,**properties): 
    """
    writes final particle list to LOCUST format
    
    notes:

    """

    pass

################################################################## Final_Particle_List class

class Final_Particle_List(classes.base_output.LOCUST_output):
    """
    class describing final particle list output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'final particle list'
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
        my_final_particle_list['status_flags'] contains a guide to the values a particle's status flag may contain 
    """

    LOCUST_output_type='final particle list'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read final particle list from file 

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
                self.data=read_final_particle_list_LOCUST(self.filepath,**properties)

        elif data_format=='TRANSP': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from TRANSP - filename required\n".format(self.ID),filename):

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_final_particle_list_TRANSP(self.filepath,**properties)

        elif data_format=='ASCOT': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from ASCOT - filename required\n".format(self.ID),filename):

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_final_particle_list_ASCOT(self.filepath,**properties) 

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/TRANSP/ASCOT)\n".format(self.ID))            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write final particle list to file

        notes: 
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_output_files / filename
                dump_final_particle_list_LOCUST(self.data,filepath,**properties)

        elif data_format=='TRANSP':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() to TRANSP - filename required\n".format(self.ID),filename):
                filepath=support.dir_output_files / filename
                dump_final_particle_list_TRANSP(self.data,filepath,**properties)
        
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST/TRANSP)\n".format(self.ID))

    def plot(self,grid=False,style='histogram',number_bins=20,fill=True,vminmax=None,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,status_flags=['PFC_intercept_3D'],weight=False,colmap=settings.cmap_default,colmap_val=np.random.uniform(),colfield='status_flag',line_style=settings.plot_line_style,label='',ax=False,fig=False):
        """
        plot the final particle list
         
        args:
            grid - corresponding distribution_function for matching binning axes etc. 
            style - choose from scatter or histogram
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            vminmax - set mesh Vmin/Vmax values
            axes - define plot axes
            LCFS - object which contains LCFS data lcfs_r and lcfs_z
            limiters - object which contains limiter data rlim and zlim
            real_scale - plot to Tokamak scale (requires equilibrium argument)
            status_flags - plot particles with these statuses
            weight - toggle whether to include marker weights in histograms
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            colfield - set the quantity which is associated with colmap e.g. time (defaults to status_flag, where the numerical value of the status_flag will dictate the colour)
            line_style - set 1D line style
            label - plot label for legends
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        #do some preliminary parsing of variables in case supplied as strings from command line etc.
        axes,status_flags,number_bins,vminmax,colmap_val=run_scripts.utils.literal_eval(axes,status_flags,number_bins,vminmax,colmap_val)

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
        if ndim==1: #plot 1D histograms

            for status in status_flags:
                p=np.where(self['status_flag']==self['status_flags'][status]) #find the particle indices which have the desired status_flag
                if weight:
                    self_binned,self_binned_edges=np.histogram(self[axes[0]][p],bins=number_bins,weights=self['weight'][p])                
                else:
                    self_binned,self_binned_edges=np.histogram(self[axes[0]][p],bins=number_bins)
                self_binned_centres=(self_binned_edges[:-1]+self_binned_edges[1:])*0.5
                ax.plot(self_binned_centres,self_binned,color=colmap(colmap_val),label=label,linestyle=line_style)
                ax.set_xlabel(axes[0])

        elif ndim==2: #plot 2D histograms

            for status in status_flags: #XXX THIS MIGHT BE CAUSING THE BUG FOR PLOTTING MULTIPLE STATUS FLAGS, AS AXES COULD BE RESET BETWEEN EACH PLOT
                p=np.where(self['status_flag']==self['status_flags'][status]) #find the particle indices which have the desired status_flag
                
                if style=='histogram':

                    if grid is not False: #bin according to pre-defined grid
                        if weight:
                            self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]][p],self[axes[1]][p],bins=[grid[axes[0]],grid[axes[1]]],weights=self['weight'][p])
                        else:
                            self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]][p],self[axes[1]][p],bins=[grid[axes[0]],grid[axes[1]]])
                    else:
                        if weight:
                            self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]][p],self[axes[1]][p],bins=number_bins,weights=self['weight'][p])
                        else:
                            self_binned,self_binned_x,self_binned_y=np.histogram2d(self[axes[0]][p],self[axes[1]][p],bins=number_bins)

                    #self_binned_x and self_binned_x are first edges then converted to centres
                    self_binned_x=(self_binned_x[:-1]+self_binned_x[1:])*0.5
                    self_binned_y=(self_binned_y[:-1]+self_binned_y[1:])*0.5

                    dx,dy=self_binned_x[1]-self_binned_x[0],self_binned_y[1]-self_binned_y[0]
                    ax.set_xticks(self_binned_x) #set axes ticks
                    ax.set_yticks(self_binned_y)
                    for index,xlabel,ylabel,xtick,ytick in enumerate(zip(ax.xaxis.get_ticklabels(),ax.yaxis.get_ticklabels(),ax.xaxis.get_ticklines(),ax.yaxis.get_ticklines())):
                        for label in [xlabel,ylabel,xtick,ytick]: label.set_visible(True) if (index % settings.tick_frequency==0) else label.set_visible(False)

                    self_binned_y,self_binned_x=np.meshgrid(self_binned_y-dy/2.,self_binned_x-dx/2.) #offset ticks onto bin centres

                    if vminmax:
                        vmin=vminmax[0]
                        vmax=vminmax[1]
                    else:
                        vmin=np.amin(self_binned)
                        vmax=np.amax(self_binned)

                    if fill:
                        ax.set_facecolor(colmap(vmin))
                        mesh=ax.pcolormesh(self_binned_x,self_binned_y,self_binned,cmap=colmap,vmin=vmin,vmax=vmax)
                        #ax.contourf(self_binned_x,self_binned_y,self_binned,levels=np.linspace(np.amin(self_binned),np.amax(self_binned),num=20),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=np.amin(self_binned),vmax=np.amax(self_binned))
                    else:
                        mesh=ax.contour(self_binned_x,self_binned_y,self_binned,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',linestyles=line_style,antialiased=True,vmin=vmin,vmax=vmax)
                        if settings.plot_contour_labels:
                            ax.clabel(mesh,inline=1,fontsize=10)
                        
                    if fig_flag is False:    
                        fig.colorbar(mesh,ax=ax,orientation='horizontal')
                        
                elif style=='scatter':
                    mesh=ax.scatter(self[axes[0]][p],self[axes[1]][p],c=self[colfield][p],cmap=colmap,marker='x',s=1,label=self.ID)

            if axes==['R','Z']:
                if LCFS: #plot plasma boundarys
                    ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                if limiters: #add boundaries if desired
                    ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall') 
                if real_scale is True: #set x and y plot limits to real scales
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')

            elif axes==['X','Y']:
                if 'X' not in self.data or 'Y' not in self.data:
                    try:
                        self['X'],self['Y']=processing.utils.RphiZ_to_XYZ(self['R'],self['phi'])
                    except:
                        pass
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
                if real_scale is True: #set x and y plot limits to real scales
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')

            ax.set_xlabel(axes[0])
            ax.set_ylabel(axes[1])
            
            if ax_flag is True or fig_flag is True: #return the plot object
                return mesh

        if ax_flag is False and fig_flag is False:
            plt.show()

#################################

##################################################################

###################################################################################################
