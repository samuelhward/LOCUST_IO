#plot_stage_2.py
 
"""
Samuel Ward
15/05/21
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
    import os,itertools
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

import stage_2_launch as batch_data

outputs=templates.plot_mod.get_output_files(batch_data,'fpl')


Pinj=33.e6
PFC_power=[]
fig1,ax1=plt.subplots(1)
fig2,ax2=plt.subplots(1)
currents=itertools.cycle(batch_data.parameters__currents_upper)
for output,current,col_val in zip(outputs,currents,np.linspace(0,1,len(batch_data.args_batch['LOCUST_run__dir_output']))):
    if output:
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        PFC_power.append([1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron])
        output['E']/=1000. #convert to keV
        output['weight'][i]=1.e6*output['f']*(output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i]*0.5*constants.mass_deuteron/Pinj
        output.plot(fig=fig1,ax=ax1,axes=['time'],fill=False,label=f'{current/1000.}kAt',colmap=settings.cmap_default,colmap_val=col_val,number_bins=200,weight=True)
        output.plot(fig=fig2,ax=ax2,axes=['E'],fill=False,label=f'{current/1000.}kAt',colmap=settings.cmap_default,colmap_val=col_val,number_bins=200,weight=True)
    else:
        PFC_power.append([-10.])

ax1.legend()
ax2.legend()
ax1.set_xlabel('Time [s]')
ax2.set_xlabel('Energy [keV]')
ax1.set_ylabel('Deposited power lost [%]',fontsize=20)
ax2.set_ylabel('Deposited power lost [%]',fontsize=20)

plt.show()

PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__kinetic_profs_Pr),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__currents_upper))

colours=['r','b'] #one colour per mode number
fig,axs=plt.subplots(len(batch_data.parameters__kinetic_profs_Pr),constrained_layout=True)
for plasma_state_counter,(ax,power,Pr,tftE) in enumerate(zip([axs],PFC_power,batch_data.parameters__kinetic_profs_Pr,batch_data.parameters__kinetic_profs_tF_tE)):
    for mode_number_counter,(mode_number,colour) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,colours)):
        n=int(np.abs(mode_number[0]))
        relative_phases_upper_middle=batch_data.parameters__phases_uppers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        relative_phases_lower_middle=batch_data.parameters__phases_lowers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        ax.scatter(batch_data.parameters__currents_upper/1000.,100.*PFC_power[plasma_state_counter][mode_number_counter]/Pinj,color=colour,label=f'$n$ = {mode_number}, U:M={int(n*relative_phases_upper_middle[0])}, L:M={int(n*relative_phases_lower_middle[0])}, M={int(n*batch_data.parameters__phases_middles[mode_number_counter][0])}')

ax.set_xlabel('ECC current amplitude [kAt]',fontsize=20) #\Phi for absolute
ax.set_ylabel('Deposited power lost [%]',fontsize=20)
ax.legend()
plt.show()


"""
outputs=templates.plot_mod.get_output_files(batch_data,'dfn')

fig,ax=plt.subplots(1)
for output,current,col_val in zip(outputs,batch_data.parameters__currents_upper,np.linspace(0,1,len(batch_data.args_batch['LOCUST_run__dir_output']))):
    if output: output.plot(fig=fig,ax=ax,axes=['R'],label=current,colmap=settings.cmap_default,colmap_val=col_val,number_bins=200)
ax.legend()
plt.show()

#cycle through poincare maps
outputs=templates.plot_mod.get_output_files(batch_data,'poinc')
axisymm=None
for run_number,output in enumerate(outputs):
    equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
    if 'B3D_EX' in batch_data.args_batch['LOCUST_run__flags'][run_number] and output is not None: #this is the 2D case    
        fig,ax=plt.subplots(1)
        output.plot(LCFS=equilibrium,limiters=equilibrium,ax=ax,fig=fig)
        equilibrium.plot(fill=False,ax=ax,fig=fig)
        plt.show()
"""
#################################
 
##################################################################
 
###################################################################################################