#plot_loss_grid_paper.py
 
"""
Samuel Ward
28/05/20
----
script for plotting 2D leak rate grid
for paper 3
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

title='X-point displacement and NBI deposited power loss (%)'

rundata=templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'rund',batch_data,processes=32,chunksize=1)
quantity=[]
components={
            '1': 'first wall',
            '2': 'dome',
            '3': 'inner vertical supports',
            '4': 'outer vertical supports',
            '5': 'horizontal supports',
            '6': 'inner plate',
            '7': 'outer plate',
            '8': 'outer pipes',
            '9': 'inner pipes',
            '10': 'outer baffle',
            '11': 'inner baffle',
            '12': 'base',
            } # component IDs
component_number='10'
for run in rundata:
    try:
        PFC_power=1.e6*float(dict(run['data']['write_vtk']['Integrated power to component'])[component_number].strip('MW'))
        quantity.append(PFC_power)
        #quantity.append(run['PFC_power']['total'])
    except:
        quantity.append(0.)

#quantity=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,'fpl',batch_data,processes=4,chunksize=1)
quantity=np.array(quantity).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_upper),
    len(batch_data.parameters__phases_lower))

fig=plt.figure()
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
ax_total.set_xlabel('Lower row phase shift $\Phi_{\mathrm{L}}$ [deg]')
ax_total.set_ylabel('Upper row phase shift $\Phi_{\mathrm{U}}$ [deg]')

ax=[]
for ax_counter in range(len(batch_data.parameters__databases)*(len(batch_data.parameters__toroidal_mode_numbers)+int(len(batch_data.parameters__toroidal_mode_numbers)/2.))):
    ax.append(fig.add_subplot(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers)+int(len(batch_data.parameters__toroidal_mode_numbers)/2.),ax_counter+1))

ax=np.array(ax).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers)+int(len(batch_data.parameters__toroidal_mode_numbers)/2.))

counter=0 #counter over different scenarios
spectra_per_fundamental=2 #how many spectra are we looping over per fundamental mode number? e.g. [[3],[3,6]] is 2 spectra
for plasma_state_counter,database in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        dLPHASE,dUPHASE=batch_data.parameters__phases_lowers[mode_number_counter][1]-batch_data.parameters__phases_lowers[mode_number_counter][0],batch_data.parameters__phases_uppers[mode_number_counter][1]-batch_data.parameters__phases_upper[0]
        if len(mode_number)==1:
            n=int(np.abs(mode_number[0]))
            case=int(batch_data.parameters__databases[plasma_state_counter][-1])
            plot_X_point_displacement(case=case,n=n,LVV=90,fig=fig,ax=ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)],colourbar=False,coord_system='ITER', #use MARS system because contours might be generated differently if using ITER system, so points might not match up
            plot_points=False
            )
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),dLPHASE/2.+np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),dUPHASE/2.+np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlabel('')
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylabel('')
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xticks(batch_data.parameters__phases_lowers[mode_number_counter][::4])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_yticks(batch_data.parameters__phases_uppers[mode_number_counter][::4])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),dLPHASE/2.+np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),dUPHASE/2.+np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].set_title(f'case {case}, $n={n}$',fontsize=13)
            ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)].tick_params(
                axis='both',          
                which='both',      
                labelbottom=False,
                labelleft=True,
                )
        LPHASE,UPHASE=np.append(batch_data.parameters__phases_lowers[mode_number_counter],batch_data.parameters__phases_lowers[mode_number_counter][-1]+dLPHASE),np.append(batch_data.parameters__phases_uppers[mode_number_counter],batch_data.parameters__phases_uppers[mode_number_counter][-1]+dUPHASE)
        #settings.discrete_colmap(colmap_name='Greys',face_colour='red',number_bins=11)
        mesh=ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].pcolormesh(LPHASE-dLPHASE/2.,UPHASE-dUPHASE/2.,quantity[plasma_state_counter,mode_number_counter],cmap='Greys',edgecolor='none',antialiased=True,
        vmin=np.min(quantity[plasma_state_counter,2*int(mode_number_counter/2):2*int(mode_number_counter/2)+1]),
        vmax=np.max(quantity[plasma_state_counter,2*int(mode_number_counter/2):2*int(mode_number_counter/2)+1])) #have same min/max for all modes within scenario 
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xticks(batch_data.parameters__phases_lowers[mode_number_counter][::4])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_yticks(batch_data.parameters__phases_uppers[mode_number_counter][::4])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xlim([np.min(batch_data.parameters__phases_lowers[mode_number_counter]),dLPHASE/2.+np.max(batch_data.parameters__phases_lowers[mode_number_counter])])
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_ylim([np.min(batch_data.parameters__phases_uppers[mode_number_counter]),dUPHASE/2.+np.max(batch_data.parameters__phases_uppers[mode_number_counter])])
        #ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].text(x=0.05,y=0.1,s=r'Pr={}, $\tau_{{\Phi}}/\tau_{{\mathrm{{E}}}}={}, \mathrm{{NBI}}_{{1}},\mathrm{{NBI}}_{{2}}$=[{},{}], $n$={}'.format(batch_data.parameters__kinetic_profs_Pr[0],batch_data.parameters__kinetic_profs_tF_tE[0],batch_data.config_beam_1,batch_data.config_beam_2,mode_number),horizontalalignment='left',transform=ax[plasma_state_counter,mode_number_counter].transAxes,color=settings.cmap_w(0.))
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_xlabel('')
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_ylabel('')
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].set_title(f'$n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}',fontsize=10)
        ax[plasma_state_counter,mode_number_counter+int(mode_number_counter/2)+1].tick_params(
            axis='both',          
            which='both',      
            labelbottom=False,
            labelleft=False,
            )

for axis in ax[-1,:]:
    axis.tick_params(
        axis='both',          
        which='both',      
        labelbottom=True,
        )

#fig.suptitle(title)
#fig.colorbar(mesh,ax=ax.ravel().tolist())
#fig.set_tight_layout(True)  
fig.subplots_adjust(left=0.1,bottom=0.1)
plt.show()

#################################

##################################################################

###################################################################################################