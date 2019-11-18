#phase_shift_parameter_scan_plot.py
 
"""
Samuel Ward
16/11/2019
----
script to plot phase shift scan
---
usage:
 
notes:         

to do:
---
"""

import context
from classes.output_classes.distribution_function import Distribution_Function
from classes.output_classes.particle_list import Final_Particle_List
import glob
import support
import processing.utils
import numpy as np
import pathlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import copy
import ast

def parse_scan_directories(parent_path='RMP_phase_scan',return_filetype='distribution_function'):
    """
    for a parameter scan grabs all filepaths and determines parameters varied over too 
    
    notes:
        assumes the <parameter_name>_<parameter_value>... structure
    returns:
        filepaths, parameters - list of desired filepaths corresponding to LOCUST output, list of dicts holding runtime parameters
    """
    #now each code directory structure begins to diverge
    #find all sub directories for this run, for example directories containing time slices in the case for LOCUST
    parameters=[]
    filepaths=[]
    filetype_dispatch={}
    filetype_dispatch['distribution_function']='F_*.dfn'
    filetype_dispatch['particle_list']='ptcl_cache.dat'
    dirs=list((support.dir_output_files / parent_path ).glob('*')) #dirs representing each value of parameter that is varied
    for dir_ in dirs: #cycling through dirs         
        split_dir=str(dir_.parts[-1]).split('_')       
        _={}
        for parameter_names,parameter_values in zip(split_dir[::2],split_dir[1::2]):
            _[parameter_names]=ast.literal_eval(parameter_values)
        parameters.append(_)
        filepath=list(dir_.glob(filetype_dispatch[return_filetype]))[0]
        filepaths.append(filepath.relative_to(support.dir_output_files))
    return filepaths,parameters

def filter_simulations(filepaths,parameters,**criteria):
    """
    extract simulations from results of parse_scan_directories based on some criteria
    notes:
    args:
        filepaths, parameters - list of desired filepaths corresponding to LOCUST output, list of dicts holding runtime parameters
        criteria - criteria to filter according, if list type then list holds allowed values to e.g. resistive=True, phase=[0,10,20]
    returns:
        filepaths, parameters - list of desired filepaths corresponding to LOCUST output, list of dicts holding runtime parameters
    """

    filepaths_=copy.deepcopy(filepaths)
    parameters_=copy.deepcopy(parameters)

    for criteria_key,criteria_value in criteria.items():
        for counter,(filepath,parameter) in reversed(list(enumerate(zip(filepaths,parameters)))): #traverse backwards since we are deleting items as we go
            if type(criteria_value)==type([]): #if criteria is given as a list then multiple options allowed as defined in list
                if any([parameter[criteria_key]==criteria for criteria in criteria_value]):
                    pass
                else:             
                    del(filepaths_[counter])
                    del(parameters_[counter])   
            else:
                if parameter[criteria_key]!=criteria_value:
                    del(filepaths_[counter])
                    del(parameters_[counter])
    return filepaths_,parameters_

def draw_frame(output,filenames,fig,ax1,ax2,axes):
    """
    notes:
    args:
        output - 
        filenames - list of all filenames in this parameter sweep, so we know where we are at in the cycle 
        fig - fig object
        ax1 - axis object
        axes - .plot() axes arg
    """

    for counter,filename in enumerate(filenames): #find where we are in the cycle
        if filename==output.filename:
            colmap_value=np.linspace(0.,1.,len(filenames))[counter] #set colourmap according to where we are in cycle

    ax1.collections=[]
    ax2.cla()

    mesh=output.plot(fig=fig,ax=ax1,axes=axes,number_bins=100,colmap_val=colmap_value,real_scale=True)

    cbar=fig.colorbar(mesh,cax=ax2,orientation='vertical')
    ax2.set_xlim([0,20])
    ax2.tick_params(axis="y",direction="in", pad=-150,length=0)
    
    ax1.set_title('')
    ax1.set_xlabel(axes[0])
    #ax1.set_ylabel(axes[1])
    #ax1.set_xlim(-1.1*np.pi,1.1*np.pi)
    ax1.set_xlim(1.,2.3)
    ax1.set_ylim(0,35)

def yield_outputs(filenames,output_type=Distribution_Function,**init):
    """
    generator to yield desired output object
    notes:
        assumes all read settings are the same for each object i.e. only filename varies
    args:
        filenames - list of target filenames
        output_type - corresponding LOCUST_IO class
        **init - kwargs to be passed to output_type.__init__(**init) constructor 
    """
    for filename in filenames:
        yield output_type(filename=filename,**init)


fig,(ax1,ax2)=plt.subplots(1,2)
filenames,parameters=parse_scan_directories(return_filetype='particle_list')
#filenames,parameters=filter_simulations(filenames,parameters,phase=0) #optionally filter out some simulations
phases=[float(parameter['phase']) for parameter in parameters]
phases,filenames,parameters=processing.utils.sort_arrays(np.array(phases),np.array(filenames),np.array(parameters))
axes=['R','Z']
outputs=yield_outputs(filenames,output_type=Final_Particle_List,ID='',data_format='LOCUST',compression=True,coordinates=axes+['status_flag'])
animation=FuncAnimation(fig,draw_frame,frames=outputs,fargs=[filenames,fig,ax1,ax2,axes],repeat=True) #cycle through phases and make animation
#plt.show()
animation.save('RMP_phase_scan_RZ.gif',writer='pillow')


#################################
 
##################################################################
 
###################################################################################################
