#plot_scan_current.py
 
"""
Samuel Ward
07/11/21
----
script for plotting scan_current 
---
 
notes:         
---
"""

###################################################################################################
#Preamble
 
import sys

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

import scan_current_launch as batch_data

Pinj=33.e6 #RMS is 0.13MW, since fires for 0.1s every 1.4s
PFC_power=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,'rund',batch_data,processes=32,chunksize=1) #if using fpl remember to re-scale beam species!
PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__currents_upper))

fig,ax=plt.subplots(1,constrained_layout=False)
colours=plt.rcParams['axes.prop_cycle'].by_key()['color']
for scenario_counter,scenario in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        linestyle='solid' if np.abs(mode_number[0])==3 else 'dashed'
        label=f'{scenario.split("_")[-1]}'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
        ax.plot(
            batch_data.parameters__currents_upper/1.e3,
            100.*PFC_power[scenario_counter,mode_number_counter]/Pinj,
            label=label,
            linestyle=linestyle,
            color=colours[int(2*scenario_counter)+2+mode_number_counter],
        )

ax.set_xlabel('ECC current [kAt]')
ax.set_ylabel('Deposited power lost [%]')
fig.subplots_adjust(right=0.85)
ax.legend(loc='center',bbox_to_anchor=(1.1,0.5),ncol=1)
plt.show()


#################################
 
##################################################################
 
###################################################################################################