#plot_stage_1_10.py
 
"""
Samuel Ward
05/05/21
----
script for plotting stage 1.10 of Simon's studies 
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
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])

import stage_1_10_launch as batch_data

#outputs=templates.plot_mod.get_output_files(batch_data,'fpl')

Pinj=33.e6
def calc_PFC_power(filename,classtype='fpl'):
    fpl=templates.plot_mod.read_locust_io_obj(filename=filename,classtype=classtype)
    if fpl:
        i=np.where(fpl['status_flag']=='PFC_intercept_3D')[0]
        PFC_power=1.e6*fpl['f']*np.sum((fpl['V_R'][i]**2+fpl['V_phi'][i]**2+fpl['V_Z'][i]**2)*fpl['FG'][i])*0.5*constants.mass_deuteron
    else:
        PFC_power=-1
    return PFC_power  

PFC_power=templates.plot_mod.apply_func_parallel(calc_PFC_power,'fpl',batch_data,processes=24,chunksize=1)
PFC_power=np.array(PFC_power).reshape(
                    len(batch_data.parameters__kinetic_profs_Pr),
                    len(batch_data.parameters__currents_upper),
                    len(batch_data.parameters__toroidal_mode_numbers),
                    len(batch_data.parameters__phases_upper)
                    )

# one figure per plasma
fig,axs=plt.subplots(len(batch_data.parameters__kinetic_profs_Pr),constrained_layout=True)
colours=[cmap_r(0.),cmap_r(0.),cmap_b(0.),cmap_b(0.)] #one colour per mode number
colours=['r','r',cmap_b(0.),cmap_b(0.)] #one colour per mode number

for plasma_state_counter,(ax,Pr,tftE) in enumerate(zip([axs],batch_data.parameters__kinetic_profs_Pr,batch_data.parameters__kinetic_profs_tF_tE)):
    for current_counter,(current,linestyle) in enumerate(zip(batch_data.parameters__currents_upper,['solid','solid'])):
        for mode_number_counter,(mode_number,colour) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,colours)):
            marker='.' if len(mode_number)>1 else 'x' #XXX if only one value of current, here override what linestyle denotes

            label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
            #label+=', '+r'$\Delta\Phi_{\mathrm{u,m}}$'
            #label+=f'={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}'
            #label+=', '+r'$\Delta\Phi_{\mathrm{l,m}}$='
            #label+=f'{int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}'
            
            #find relative phase between coils for this toroidal mode number
            relative_phases_upper_middle=batch_data.parameters__phases_uppers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
            relative_phases_lower_middle=batch_data.parameters__phases_lowers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
            ax.scatter(batch_data.parameters__phases_middles[mode_number_counter],100*PFC_power[plasma_state_counter,current_counter,mode_number_counter]/Pinj,label=label,color=colour,linestyle=linestyle,marker=marker,s=80)

ax.legend()
ax.set_xlabel('$\Phi_{\mathrm{m}}$ [deg]') #\Phi for absolute
ax.set_ylabel('NBI power loss [%]')
# remove ticks from total ax
#ax_total = fig.add_subplot(111,frameon=False)
#ax_total.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False)

plt.show()

#################################
 
##################################################################
 
###################################################################################################