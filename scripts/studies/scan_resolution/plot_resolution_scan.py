#plot_resolution_scan.py
 
"""
Samuel Ward
28/05/20
----
script resolution scan
---
 
notes:         
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    import os
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
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

try:
    cwd=pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    sys.path.append(str(cwd.parents[1]))
    import templates.plot_mod
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])

import scan_resolution_launch as batch_data

""" 
outputs=templates.plot_mod.get_divergence_files(batch_data)
for run_number,output in enumerate(outputs):
    if run_number==1: continue
    if output is not None:
        field_data_RZ,field_data_XY=output
        p=np.where((np.abs(field_data_RZ['dB_field_R'])>1.0) | (np.log(np.abs(field_data_RZ['divB']))<-10.))
        field_data_RZ['divB'][p]=1.e-3
        field_data_RZ['divB']=np.log10(np.abs(field_data_RZ['divB']))
        field_data_RZ.plot(key='divB')
"""

outputs=templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'fpl',batch_data,processes=16,chunksize=1)
fig,ax1=plt.subplots(1,1,constrained_layout=True)
Pinj=33.e6
Einj=1.e6
for run_number,output in enumerate(outputs):
    if output is not None:
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        PFC_power=100.*1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron/Pinj
        #weight_factor=Pinj/np.sum(output['weight']*Einj*constants.charge_e)
        #print(1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron)
        #print(np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron)
        #print(batch_data.parameters__perturbation_resolutions_R[run_number])
        if 'B3D_EX' not in batch_data.args_batch['LOCUST_run__flags'][run_number]: #this is the 2D case
            ax1.axhline(PFC_power,color='black',label='2D field losses',linewidth=2)
        else:
            ax1.scatter(np.log10(batch_data.parameters__perturbation_resolutions_R[run_number]),PFC_power,color='black',marker='x',linestyle='solid',alpha=1,linewidth=2,s=80)
ax1.set_xlabel("")
ax1.set_xlabel(r"Perturbation grid spacing $\Delta x_{\widetilde{\mathbf{B}}}$ (log$_{10}$) [m]")
ax1.set_ylabel("NBI power loss [%]")
ax1.legend()
#ax1.tick_params(
#    axis='x',          
#    which='both',      
#    bottom=False,      
#    top=False,         
#    labelbottom=False) 


# plot divergence here

#static length scales
equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][0],GEQDSKFIX1=True,GEQDSKFIX2=True)
grid_spacing_equilibrium=np.abs(equilibrium['R_1D'][1]-equilibrium['R_1D'][0])
boris_step_length=1.e-9*np.sqrt(2.*1.e6*constants.charge_e/constants.mass_deuteron)
print(boris_step_length)

#static B scales
B_0=np.abs(equilibrium['bcentr'])

ax2=ax1.twinx()
outputs=templates.plot_mod.get_divergence_files(batch_data)
for run_number,output in enumerate(outputs):
    if output is not None:
        field_data_RZ,field_data_XY=output
        #remove values outside LCFS since they can dominate errors
        #p=np.where(np.abs(field_data_RZ['dB_field_R'])<1.0)
        p=processing.utils.within_LCFS(
            field_data_RZ['R'].flatten(),field_data_RZ['Z'].flatten(),equilibrium).reshape(
            len(field_data_RZ['R_1D']),len(field_data_RZ['Z_1D']))
        
        for scaling_factor,scaling_label,marker in zip(
            [batch_data.parameters__perturbation_resolutions_R[run_number]/B_0,boris_step_length/B_0],
        [r'$A=\Delta x_{\widetilde{\mathbf{B}}}/B_{\mathrm{0}}$','$A=\Delta x_{\mathrm{Boris}}/B_{\mathrm{0}}$'],
        ['x','.'],
        ):
            divB=np.log10(scaling_factor*np.mean(np.abs(field_data_RZ['divB'][p])))
            if 'B3D_EX' not in batch_data.args_batch['LOCUST_run__flags'][run_number]: #this is the 2D case
                pass#ax2.axhline(divB,color=cmap_r(0.),label='2D case')
            else:
                label=scaling_label if run_number == 5 else ''
                ax2.scatter(np.log10(batch_data.parameters__perturbation_resolutions_R[run_number]),divB,color=cmap_b(0.),marker=marker,linestyle='solid',alpha=1,linewidth=2,s=80,label=label)
                #ax2.scatter(np.log10(batch_data.parameters__perturbation_resolutions_R[run_number]),divB,color='black',marker='x',linestyle='-',alpha=1,linewidth=2)
ax2.set_ylabel(r"Mean $\left|\nabla\cdot A\mathbf{B}\right|$ (log$_{10}$)",color=cmap_b(0.))
ax2.tick_params('y',colors=cmap_b(0.))
ax2.legend(bbox_to_anchor=(.0,.9),ncol=1,loc='upper left')
plt.show()





#cycle through poincare maps
outputs=templates.plot_mod.get_output_files(batch_data,'poinc')
for run_number,output in enumerate(outputs):
    equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
    fig,ax=plt.subplots(1)
    if 'B3D_EX' not in batch_data.args_batch['LOCUST_run__flags'][run_number] and output is not None: #this is the 2D case    
        equilibrium.plot()
    elif output is not None:
        output.plot(ax=ax,fig=fig,LCFS=equilibrium,limiters=equilibrium,real_scale=True)
        ax.set_xlim([4.5,6.0])
        ax.set_ylim([-4.0,-1.3])
    plt.show()

'''
outputs=templates.plot_mod.get_divergence_files(batch_data)
fig,ax=plt.subplots(1)
def plot_divergence(output,fig,ax): 
    try:
        if output:   
            ax.cla()
            field_data_RZ,field_data_XY=output
            field_data_RZ.plot(key='divB',fill=False,ax=ax,fig=fig)
        else:
            ax.cla()
    except:
        pass
animation=FuncAnimation(fig,plot_divergence,frames=outputs,fargs=[fig,ax],repeat=True,interval=1)
plt.show()
animation.save('divergence_animation.gif',writer='pillow')
'''

#################################
 
##################################################################
 
###################################################################################################