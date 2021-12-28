#plot_XPD_correlation_paper.py
 
"""
Samuel Ward
28/05/20
----
script for correlating particle loss vs XPD
paper_3
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
    import scipy.stats
    import matplotlib,itertools
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import scan_upper_lower_all_scenarios_launch as batch_data

Pinj=33.e6

def calculate_quantity(rundata,component_numbers=None):
    quantity=[]
    for run in rundata:
        try:
            if component_numbers:
                PFC_power=0
                for component_number in component_numbers:
                    PFC_power+=1.e6*float(dict(run['data']['write_vtk']['Integrated power to component'])[str(component_number)].strip('MW'))
            else:
                PFC_power=run['PFC_power']['total']
            PFC_power*=100./Pinj
            quantity.append(PFC_power)
        except:
            quantity.append(0.)
    quantity=np.array(quantity).reshape(
        len(batch_data.parameters__databases),
        len(batch_data.parameters__toroidal_mode_numbers),
        len(batch_data.parameters__phases_upper),
        len(batch_data.parameters__phases_lower))
    return quantity

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



markers = itertools.cycle((
    '+',
    '.',
    '*',
    'x',
    ))
    
#number_colours=len(batch_data.parameters__databases)
#colours=matplotlib.cm.get_cmap('tab10',10)
#colours=itertools.cycle(list(colours(np.arange(number_colours)))[0:number_colours])
#next(colours)

rundata=templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'rund',batch_data,processes=32,chunksize=1)

#correlation plot
quantity=calculate_quantity(rundata,component_numbers=[2])
quantity=calculate_quantity(rundata)
fig=plt.figure(figsize=(20, 10))
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
top_ax=fig.add_subplot(2,1,1)
bottom_axs=[]
bottom_axs.append(fig.add_subplot(2,2,3))
bottom_axs.append(fig.add_subplot(2,2,4))
colours=plt.rcParams['axes.prop_cycle'].by_key()['color']
for plasma_state_counter,plasma_state in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        marker=next(markers)
        n=int(np.abs(mode_number[0]))
        #linestyle='dashed' if len(mode_number)>1 else 'solid' 
        bottom_ax=1 if len(mode_number)>1 else 0
        case=int(batch_data.parameters__databases[plasma_state_counter][-1])
        PHI_L,PHI_U,XPD_mars=plot_X_point_displacement(case=case,n=n,LVV=90,return_grid=True,coord_system='ITER')
        XPD_interpolator=processing.utils.interpolate_2D(
                        np.sort(PHI_L[0,:]),
                        np.sort(PHI_U[:,0]),
                        XPD_mars
                        )
        XPD=XPD_interpolator(
                batch_data.parameters__phases_lowers[mode_number_counter], 
                batch_data.parameters__phases_uppers[mode_number_counter],
            )
        label=f'{plasma_state.split("_")[-1]} ({batch_data.configs_beam_species[plasma_state_counter][0].capitalize()})'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'        
        top_ax.scatter(
        XPD.T.flatten(),
        quantity[plasma_state_counter,mode_number_counter].flatten(),
        label=label,
        color=colours[int(2*plasma_state_counter)+2+int(mode_number_counter/2)],
        marker=marker,
        )
        m,c,r,p,err=scipy.stats.linregress(
        XPD.T.flatten(),
        quantity[plasma_state_counter,mode_number_counter].flatten(),
        )
        bottom_axs[bottom_ax].plot(
            np.linspace(np.min(XPD),np.max(XPD),10),
            m*np.linspace(np.min(XPD),np.max(XPD),10)+c,
            marker=marker,
            color=colours[int(2*plasma_state_counter)+2+int(mode_number_counter/2)],
            label=label,
            linewidth=1
        #    linestyle=linestyle,
        )
        m,c,r,p,err=scipy.stats.linregress(
        XPD.T.flatten(),
        quantity[plasma_state_counter,mode_number_counter].flatten(),
        )
        Pearson=scipy.stats.pearsonr(
            XPD.T.flatten(),
            quantity[plasma_state_counter,mode_number_counter].flatten(),)
        print(f'Pearson correlation for {label}={Pearson}')
        print(f'Gradient for {label}={m}')

for ax in bottom_axs:
    ax.set_ylim([0,15.5])

bottom_axs[-1].tick_params(top=False, bottom=True, left=False, right=False)
ax_total.set_ylabel('NBI power loss [%]')
ax_total.set_xlabel('X-point displacement [mm]')
top_ax.legend(bbox_to_anchor=(1.15,.5),ncol=1,loc='center',bbox_transform=ax_total.transAxes)
fig.subplots_adjust(right=0.8)

#plt.savefig(
#    fname=f'/home/ITER/wards2/pics/paper_3_XPD_correlation.pdf',
#    dpi=800,
#)

plt.show()




# component-resolved Pearson correlation vs component
fig,ax=plt.subplots(1)
for plasma_state_counter,plasma_state in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        marker=next(markers)
        for component_number,component in components.items():
            quantity=calculate_quantity(rundata,component_numbers=[component_number])
            n=int(np.abs(mode_number[0]))
            case=int(batch_data.parameters__databases[plasma_state_counter][-1])
            PHI_L,PHI_U,XPD_mars=plot_X_point_displacement(case=case,n=n,LVV=90,return_grid=True,coord_system='ITER')
            XPD_interpolator=processing.utils.interpolate_2D(
                            np.sort(PHI_L[0,:]),
                            np.sort(PHI_U[:,0]),
                            XPD_mars
                            )
            XPD=XPD_interpolator(
                    batch_data.parameters__phases_lowers[mode_number_counter], 
                    batch_data.parameters__phases_uppers[mode_number_counter],
                )
            label=f'{plasma_state.split("_")[-1]}'# ({batch_data.configs_beam_species[plasma_state_counter][0].capitalize()})'
            label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
            label+=f', {component}'
            Pearson=scipy.stats.pearsonr(
                XPD.T.flatten(),
                quantity[plasma_state_counter,mode_number_counter].flatten(),
                )
            m,c,r,p,err=scipy.stats.linregress(
            XPD.T.flatten(),
            quantity[plasma_state_counter,mode_number_counter].flatten(),
            )
            label+=f': {Pearson}'
            label+=f' : {m}'
            print(label)
            ax.scatter(
                int(component_number),
                #Pearson[0],
                m,
                marker=marker, 
                color=colours[int(2*plasma_state_counter)+2+int(mode_number_counter/2)],
            )
            #try:
            #except:
            #    pass

plt.show()

#################################

##################################################################

###################################################################################################