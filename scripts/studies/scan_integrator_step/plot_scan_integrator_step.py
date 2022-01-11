#plot_integrator_step_scan.py
 
"""
Samuel Ward
17/08/20
----
check convergence of integrator step size
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

import scan_integrator_step_launch as batch_data

outputs=templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'fpl',batch_data,processes=16,chunksize=1)

fig=plt.figure(constrained_layout=False)
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
ax1=fig.add_subplot(2,1,1)
ax2=fig.add_subplot(2,1,2)

Pinj=33.e6
for counter,(output,col_val) in enumerate(zip(outputs,np.linspace(0,0.9,len(batch_data.args_batch['LOCUST_run__dir_output'])))):
    if output: 
        if batch_data.args_batch['LOCUST_run__flags'][counter]['LEIID']==6:
            colmap=settings.cmap_inferno
        elif batch_data.args_batch['LOCUST_run__flags'][counter]['LEIID']==7:
            colmap=settings.cmap_inferno
        number_markers=len(output['weight'])
        print(f'number markers={number_markers}')
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        UNBOR=batch_data.args_batch['LOCUST_run__flags'][counter]['UNBOR']
        print(f'UNBOR={UNBOR}, #lost={len(i)}')
        PFC_power=1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron
        print(f'PFC power = {100*PFC_power/Pinj}%')
        print(f"f={output['f']}")
        
        
        label=''
        #label+=r'$\mathrm{log}_{10}(\Delta t)=$'
        #label+=f'{format(np.log10(1.e-7/UNBOR),".2f")}, '
        label+=r'$\mathrm{n}_{\mathrm{coll}}=$'
        label+=f'{UNBOR}'
        label+=r'$, \mathrm{P}_{\mathrm{loss}}=$'
        label+=f'{format(100*PFC_power/Pinj,".2f")}%'
        output['weight']*=1.e6*output['f']*output['E']*constants.charge_e/1.e3 #plot loss power of markers in kW
        output.plot(fig=fig,ax=ax1,axes=['time'],fill=False,label=label,colmap=colmap,colmap_val=col_val,number_bins=200,weight=True)
        output['E']/=1000. #convert to keV
        output.plot(fig=fig,ax=ax2,axes=['E'],fill=False,label=label,colmap=colmap,colmap_val=col_val,number_bins=200,weight=True)
ax1.legend()
ax1.set_xlabel('Time [s]')
ax2.set_xlabel('Energy [keV]')
ax_total.set_ylabel('Loss power [a.u.]')
ax1.set_title('')
ax2.set_title('')
plt.show()

#################################
 
##################################################################
 
###################################################################################################