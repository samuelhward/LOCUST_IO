#plot_scan_DNB_phase.py
 
"""
Samuel Ward
07/11/21
----
script for plotting scan_DNB_phase 
---
 
notes:         
fpl is preferred since some output rundata files cannot be read by python readlines for some reason
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
    import scipy.constants
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

import scan_DNB_phase_launch as batch_data

Pinj=.13e6 #RMS is 0.13MW, since fires for 0.1s every 1.4s
PFC_power=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,'fpl',batch_data,processes=32,chunksize=1)/2. #if using fpl remember to re-scale beam species!
PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

fig,ax=plt.subplots(1,constrained_layout=True)
colours=plt.rcParams['axes.prop_cycle'].by_key()['color']
for scenario_counter,scenario in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        label=f'{scenario.split("_")[-1]}'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
        ax.plot(batch_data.parameters__phases_middles_cases_all[scenario_counter][mode_number_counter],100.*PFC_power[scenario_counter,mode_number_counter]/Pinj,
        label=label,
        color=colours[int(2*scenario_counter)+2+mode_number_counter],
        )

ax.set_xlabel('$\Delta\Phi$ [deg]') #\Phi for absolute
ax.set_xlabel('Absolute phase shift of RMP ($\Phi_{\mathrm{m}}$) [deg]') #\Phi for absolute
ax.set_ylabel("NBI power loss [%]")
fig.subplots_adjust(right=0.85)
ax.legend(loc='center',bbox_to_anchor=(0.5,1.1),ncol=4)
plt.show()
""" 
 
eqs=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'eq',batch_data,processes=32,chunksize=1,GEQDSKFIX1=True,GEQDSKFIX2=True))
fpls=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'fpl',batch_data,processes=32,chunksize=1))
beam_depos=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'beamdepo',batch_data,processes=32,chunksize=1))

for fpl,eq,beam_depo in zip(fpls,eqs,beam_depos):
    fpl.set(V_pitch_initial_2D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in fpl.data.items() if 'initial' in key},equilibria=[eq],perturbations=None,i3dr=-1,phase=0.))
    beam_depo.set(V_pitch=processing.utils.pitch_calc(beam_depo,equilibria=[eq],perturbations=None,i3dr=-1,phase=0.))

eqs=np.array(eqs).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))
fpls=np.array(fpls).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))
beam_depos=np.array(beam_depos).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

fig,ax=plt.subplots(1)
for beam_depo in beam_depos[:,0,1]:
    beam_depo.plot(axes=['V_pitch'],number_bins=50,fig=fig,ax=ax)

plt.show()
    

fig,ax=plt.subplots(1)
for fpl in fpls[:,0,1]:
    fpl.plot(axes=['V_pitch_initial_2D'],number_bins=50,fig=fig,ax=ax)

plt.show()

""" 

#################################
 
##################################################################
 
###################################################################################################