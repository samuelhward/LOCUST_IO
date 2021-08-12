#plot_loss_grid.py
 
"""
Samuel Ward
28/05/20
----
script for plotting 2D leak rate grid
---
 
notes:
    assumes fundamental-only spectrum is first in its subgroup i.e. [3],[3,6] not [3,6],[3]          
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import scan_upper_lower_all_scenarios_launch as batch_data

Pinj=33.e6

def calc_PFC_power(filename,classtype='fpl'):
    fpl=templates.plot_mod.read_locust_io_obj(filename=filename,classtype=classtype)
    if fpl:
        i=np.where(fpl['status_flag']=='PFC_intercept_3D')[0]
        PFC_power=1.e6*fpl['f']*np.sum((fpl['V_R'][i]**2+fpl['V_phi'][i]**2+fpl['V_Z'][i]**2)*fpl['FG'][i])*0.5*constants.mass_deuteron
    else:
        PFC_power=-1
    return PFC_power  

def calc_mean_loss_time(filename,classtype='fpl'):
    fpl=templates.plot_mod.read_locust_io_obj(filename=filename,classtype=classtype)
    if fpl:
        i=np.where(fpl['status_flag']=='PFC_intercept_3D')[0]
        mean_loss_time=np.average(fpl['time'][i],weights=fpl['FG'][i])
    else:
        mean_loss_time=-1
    return mean_loss_time  

calc_function=calc_PFC_power
title='Mean loss time (s)'
title='NBI deposited power loss (%)'






quantity=templates.plot_mod.apply_func_parallel(calc_function,'fpl',batch_data,processes=24,chunksize=4)
quantity=np.array(quantity).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper),len(batch_data.parameters__phases_lower))

fig,ax=plt.subplots(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers)+int(len(batch_data.parameters__toroidal_mode_numbers)/2.),constrained_layout=True) #add pair of axes for ploting XPD alongside
ax=np.array(ax,ndmin=2)
counter=0 #counter over different scenarios
spectra_per_fundamental=2 #how many spectra are we looping over per fundamental mode number? e.g. [[3],[3,6]] is 2 spectra
for plasma_state_counter,database in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):

        if len(mode_number)==1:
            n=int(np.abs(mode_number[0]))
            case=int(batch_data.parameters__databases[plasma_state_counter][-1])
            plot_X_point_displacement(case=case,n=n,LVV=90,fig=fig,ax=ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)],colourbar=False,coord_system='ITER', #use MARS system because contours might be generated differently if using ITER system, so points might not match up
            plot_points=False
            ) 
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlabel('')
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylabel('')
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xticks(batch_data.parameters__phases_lowers[mode_number_counter][::2])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_yticks(batch_data.parameters__phases_uppers[mode_number_counter][::2])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),np.max(batch_data.parameters__phases_uppers[mode_number_counter])])

            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].tick_params(
                axis='both',          
                which='both',      
                labelbottom=False,
                labelleft=True,
                )

        dLPHASE,dUPHASE=batch_data.parameters__phases_lowers[mode_number_counter][1]-batch_data.parameters__phases_lowers[mode_number_counter][0],batch_data.parameters__phases_uppers[mode_number_counter][1]-batch_data.parameters__phases_upper[0]
        LPHASE,UPHASE=np.append(batch_data.parameters__phases_lowers[mode_number_counter],batch_data.parameters__phases_lowers[mode_number_counter][-1]+dLPHASE),np.append(batch_data.parameters__phases_uppers[mode_number_counter],batch_data.parameters__phases_uppers[mode_number_counter][-1]+dUPHASE)
        #settings.discrete_colmap(colmap_name='Greys',face_colour='red',number_bins=11)
        mesh=ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].pcolormesh(LPHASE-dLPHASE/2.,UPHASE-dUPHASE/2.,quantity[plasma_state_counter,mode_number_counter],cmap='Greys',edgecolor='none',antialiased=True)
        #mesh=ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].contourf(LPHASE[0:-1]-dLPHASE/2.,UPHASE[0:-1]-dUPHASE/2.,quantity[plasma_state_counter,mode_number_counter],cmap=settings.discrete_colmap(colmap_name='Greys',face_colour='black',number_bins=11),edgecolor='none',antialiased=True)#,vmin=0,vmax=10)
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xticks(batch_data.parameters__phases_lowers[mode_number_counter][::2])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_yticks(batch_data.parameters__phases_uppers[mode_number_counter][::2])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
        #ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].text(x=0.05,y=0.1,s=r'Pr={}, $\tau_{{\Phi}}/\tau_{{\mathrm{{E}}}}={}, \mathrm{{NBI}}_{{1}},\mathrm{{NBI}}_{{2}}$=[{},{}], $n$={}'.format(batch_data.parameters__kinetic_profs_Pr[0],batch_data.parameters__kinetic_profs_tF_tE[0],batch_data.config_beam_1,batch_data.config_beam_2,mode_number),horizontalalignment='left',transform=ax[plasma_state_counter,mode_number_counter].transAxes,color=settings.cmap_w(0.))
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xlabel('')
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_ylabel('')
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_title('EP losses: '+'$n$ = '+','.join(str(abs(mode)) for mode in mode_number))
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].tick_params(
            axis='both',          
            which='both',      
            labelbottom=False,
            labelleft=False,
            )


for axis in ax[:,0]:
    axis.set_ylabel('$\Delta\Phi_{\mathrm{U}}$ [deg]')

for axis in ax[-1,:]:
    axis.set_xlabel('$\Delta\Phi_{\mathrm{L}}$ [deg]')
    axis.tick_params(
        axis='both',          
        which='both',      
        labelbottom=True,
        )

fig.suptitle(title)
fig.colorbar(mesh,ax=ax.ravel().tolist())
#fig.set_tight_layout(True)  
plt.show()

#################################

##################################################################

###################################################################################################