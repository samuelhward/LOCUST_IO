#for use with studies_simon.stage1_10_1
import context,pathlib,sys,os
import stage_1_10_1_launch as batch_data
import numpy as np
import matplotlib.pyplot as plt 
import templates.plot_mod
from classes.input_classes.equilibrium import Equilibrium

equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][0],GEQDSKFIX1=True,GEQDSKFIX2=True)

poincares=np.array(list(templates.plot_mod.get_io_files(batch_data,'poinc3')),dtype='object').reshape(
    len(batch_data.parameters__currents_upper),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers[0]),
    )

fig=plt.figure(constrained_layout=False)

ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
axs=[]
axs.append(fig.add_subplot(2,2,1))
axs.append(fig.add_subplot(2,2,3))
axs.append(fig.add_subplot(2,2,2))
axs.append(fig.add_subplot(2,2,4))
axs=np.array(axs).T

phi=60 #phi slice to plot

for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers)):
    templates.plot_mod.plot_poincare_q_theta(poincares[0,mode_counter,3],equilibrium,phi=phi*np.pi/180.,ax=axs[mode_counter],fig=fig) #plot where the beam depo is max
    axs[mode_counter].set_ylim([1,4])
    axs[mode_counter].set_xlabel('')
    axs[mode_counter].set_ylabel('')
    title=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
    axs[mode_counter].set_title(title,fontsize=15)


axs[0].tick_params(axis='x',which='both',labelbottom=False)
axs[2].tick_params(axis='x',which='both',labelbottom=False)
axs[2].tick_params(axis='y',which='both',labelleft=False)
axs[3].tick_params(axis='y',which='both',labelleft=False)
ax_total.set_ylabel(r'Safety factor')#,fontsize=10)
ax_total.set_xlabel(r'Poloidal angle [deg]')#,fontsize=10)

plt.show()

plt.savefig(
    fname=f'/home/ITER/wards2/pics/paper_2_poincare_q_theta_phi_{phi}.png',
    dpi=800,
)


