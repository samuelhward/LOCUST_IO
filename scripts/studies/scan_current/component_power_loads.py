#component_power_loads.py
 
"""
Samuel Ward
30/12/21
----
script for reading component-level power loads
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
    import matplotlib
    import matplotlib.pyplot as plt
    import os
    import datetime
    import itertools
    import random
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
    from classes.output_classes.rundata import Rundata
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/rundata.py could not be imported!\nreturning\n") 
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

import scan_current_launch as batch_data


output_rundata=templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'rund',batch_data,processes=32,chunksize=1)
output_rundata=np.array(output_rundata).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__currents_upper),
    )


#component ID dict
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

marker = itertools.cycle((
    'X',
    '+',
    '.',
    'o',
    '*',
    'x',
    'P',
    'p',
    '1',
    'v',
    'D',
    's',
    )) 

number_colours=len(list(components.keys()))
colours=matplotlib.cm.get_cmap('hsv',number_colours+1)
colours=reversed(list(colours(np.arange(number_colours)))[0:number_colours])
colours=itertools.cycle(colours)
#next(colours)
#colours=matplotlib.colors.ListedColormap(colours[:-1], "")

#output_rundata[0,0,4]['data']['write_vtk']['Integrated power to component'][0]

fig=plt.figure(constrained_layout=False)
ylim=[-5.,3.5]
ylim=[0.,3.5]
axs=[]
ax_total=fig.add_subplot(1,1,1)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
num_axes=len(batch_data.parameters__databases)*len(batch_data.parameters__toroidal_mode_numbers)
counter=0
for database_counter,database in enumerate(batch_data.parameters__databases):
    for mode_counter,modes in enumerate(batch_data.parameters__toroidal_mode_numbers):
        counter+=1
        axs.append(fig.add_subplot(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),counter))

axs=np.array(axs).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
)

for database_counter,(ax_col,database) in enumerate(zip(axs,batch_data.parameters__databases)):
    for mode_counter,(ax,modes) in enumerate(zip(ax_col,batch_data.parameters__toroidal_mode_numbers)):
        for component_number in [str(number) for number in np.arange(12)+1]:

            try:
                PFC_power=np.array([float(dict(rd['data']['write_vtk']['Integrated power to component'])[component_number].strip('MW'))  if rd else 0 for rd in output_rundata[database_counter,mode_counter]])
                print(f'case={database}: n={modes}: {components[component_number]}: {np.min(PFC_power)}-{np.max(PFC_power)}MW')
                PFC_power=np.nan_to_num(np.log10(PFC_power*1.e3),nan=ylim[0])
                PFC_power[PFC_power>ylim[-1]]=ylim[-1]
                PFC_power[PFC_power<ylim[0]]=ylim[0]
                ax.plot(batch_data.parameters__currents_upper/1.e3,PFC_power,label=components[component_number],color=next(colours),marker=next(marker),markersize=10)
                title=f'{database.split("_")[-1]}'# ({batch_data.configs_beam_species[database_counter][0].capitalize()})'
                title+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in modes)}'
                ax.set_title(title,loc='left',y=1.0)
                ax.set_ylim(ylim)
                ax.set_yticks(np.arange(ylim[0],ylim[1]))
                ax.grid(True,axis='y',which='both')
            except:
                pass

for ax in axs[:,-1]:
    ax.tick_params(
        axis='y',
        which='both',
        labelleft=False
        )
for row in axs[:-1,:]: 
    for ax in row:
        ax.tick_params(
            axis='x',
            which='both',
            labelbottom=False
            )

ax_total.set_xlabel('ECC current [kAt]')
ax_total.set_ylabel(r'Component power flux [kW] ($\mathrm{log}_{10}$)')

axs[0,0].legend(loc='center',ncol=6,bbox_to_anchor=(0.5,1.1),bbox_transform=ax_total.transAxes)
plt.show()