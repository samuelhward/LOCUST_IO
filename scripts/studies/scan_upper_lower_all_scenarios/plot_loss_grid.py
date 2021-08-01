#plot_loss_grid.py
 
"""
Samuel Ward
28/05/20
----
script for plotting 2D leak rate grid
---
 
notes:
    phase shift origin moved to centre of coils as per Todd Evans convention         
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

PFC_power=templates.plot_mod.apply_func_parallel(calc_PFC_power,'fpl',batch_data,processes=16,chunksize=8)

PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper),len(batch_data.parameters__phases_lower))

fig,ax=plt.subplots(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers)*2) #double axes for plotting XPD alongside
ax=np.array(ax,ndmin=2)
counter=0 #counter over different scenarios
for plasma_state_counter,database in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        dLPHASE,dUPHASE=batch_data.parameters__phases_lowers[mode_number_counter][1]-batch_data.parameters__phases_lowers[mode_number_counter][0],batch_data.parameters__phases_uppers[mode_number_counter][1]-batch_data.parameters__phases_upper[0]
        LPHASE,UPHASE=np.append(batch_data.parameters__phases_lowers[mode_number_counter],batch_data.parameters__phases_lowers[mode_number_counter][-1]+dLPHASE),np.append(batch_data.parameters__phases_uppers[mode_number_counter],batch_data.parameters__phases_uppers[mode_number_counter][-1]+dUPHASE)
        #mesh=ax[plasma_state_counter,2*mode_number_counter].pcolormesh(LPHASE-dLPHASE/2.,UPHASE-dUPHASE/2.,100*PFC_power[plasma_state_counter,mode_number_counter]/Pinj,cmap=settings.discrete_colmap(colmap_name='Greys',face_colour='black',number_bins=11),edgecolor='none',antialiased=True,vmin=0,vmax=10)
        mesh=ax[plasma_state_counter,2*mode_number_counter].contourf(LPHASE-dLPHASE/2.,UPHASE-dUPHASE/2.,100*PFC_power[plasma_state_counter,mode_number_counter]/Pinj,cmap=settings.discrete_colmap(colmap_name='Greys',face_colour='black',number_bins=11),edgecolor='none',antialiased=True)#,vmin=0,vmax=10)
        ax[plasma_state_counter,2*mode_number_counter].set_xticks(batch_data.parameters__phases_lowers[mode_number_counter])
        ax[plasma_state_counter,2*mode_number_counter].set_yticks(batch_data.parameters__phases_uppers[mode_number_counter])
        ax[plasma_state_counter,2*mode_number_counter].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
        ax[plasma_state_counter,2*mode_number_counter].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
        ax[plasma_state_counter,2*mode_number_counter].set_xlabel('$\Delta\Phi_{\mathrm{L}}$ [deg]')
        ax[plasma_state_counter,2*mode_number_counter].set_ylabel('$\Delta\Phi_{\mathrm{U}}$ [deg]')
        ax[plasma_state_counter,2*mode_number_counter].text(x=0.05,y=0.1,s=r'Pr={}, $\tau_{{\Phi}}/\tau_{{\mathrm{{E}}}}={}, \mathrm{{NBI}}_{{1}},\mathrm{{NBI}}_{{2}}$=[{},{}], $n$={}'.format(batch_data.parameters__kinetic_profs_Pr[0],batch_data.parameters__kinetic_profs_tF_tE[0],batch_data.config_beam_1,batch_data.config_beam_2,mode_number),horizontalalignment='left',transform=ax[plasma_state_counter,mode_number_counter].transAxes,color=settings.cmap_w(0.))
        ax[plasma_state_counter,2*mode_number_counter].set_title(f'{database}')
        n=int(np.abs(mode_number[0]))
        case=int(batch_data.parameters__databases[plasma_state_counter][-1])
        plot_X_point_displacement(case=case,n=n,LVV=90,fig=fig,ax=ax[plasma_state_counter,2*mode_number_counter+1],colourbar=False,coord_system='ITER') #use MARS system because contours might be generated differently if using ITER system, so points might not match up
        ax[plasma_state_counter,2*mode_number_counter+1].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
        ax[plasma_state_counter,2*mode_number_counter+1].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
        counter+=1

fig.suptitle('Total NBI power loss (%)')
fig.colorbar(mesh,ax=ax.ravel().tolist())
#fig.set_tight_layout(True)  
plt.show()

#################################

##################################################################

###################################################################################################