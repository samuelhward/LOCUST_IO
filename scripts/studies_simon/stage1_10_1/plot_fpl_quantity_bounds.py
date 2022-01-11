def plot_fpl_q_lambda_bounds_1d(output_fpls,beam_deposition,batch_data,output_fpls_collisionless=[]):
    import matplotlib.pyplot as plt 
    import numpy as np 
    fig=plt.figure(constrained_layout=False)
    axs=[]
    ax_total=fig.add_subplot(1,1,1,frameon=False)
    ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_total.spines['top'].set_color('none')
    ax_total.spines['bottom'].set_color('none')
    ax_total.spines['left'].set_color('none')
    ax_total.spines['right'].set_color('none')
    ax_row1=fig.add_subplot(2,1,1,frameon=False)
    ax_row1.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_row1.spines['top'].set_color('none')
    ax_row1.spines['bottom'].set_color('none')
    ax_row1.spines['left'].set_color('none')
    ax_row1.spines['right'].set_color('none')
    ax_row2=fig.add_subplot(2,1,2,frameon=False)
    ax_row2.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_row2.spines['top'].set_color('none')
    ax_row2.spines['bottom'].set_color('none')
    ax_row2.spines['left'].set_color('none')
    ax_row2.spines['right'].set_color('none')
    axs.append(fig.add_subplot(2,2,1))
    axs.append(fig.add_subplot(2,2,3))
    axs.append(fig.add_subplot(2,2,2))
    axs.append(fig.add_subplot(2,2,4))
    axs=np.array(axs).reshape(2,2)
    for ax_col,output_fpl_group in zip(
        axs.T,
        [output_fpls_collisionless,output_fpls],
        ):
        for ax,quantity,grid in zip(
            ax_col,
            ['q_initial','V_pitch_initial_2D'],
            [np.linspace(2,4.05,200),np.linspace(0.34,.8,200)],
        ):
            for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower)):
                for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers)):
                    binned_fpls=[]
                    for phase_counter,(parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower)): 
                        p=np.where(output_fpl_group[current_counter,mode_counter,phase_counter]['status_flag']=='PFC_intercept_3D')[0] 
                        hist_fpl,bins_fpl=np.histogram(
                            output_fpl_group[current_counter,mode_counter,phase_counter][quantity][p],
                            bins=grid,
                            weights=output_fpl_group[current_counter,mode_counter,phase_counter]['weight'][p],
                            )
                        hist_bd,bins_bd=np.histogram(
                            beam_deposition[quantity.replace('_initial','')],
                            bins=grid,
                            weights=beam_deposition['weight'],
                            )
                        bins_fpl_centres=(bins_fpl[:-1]+bins_fpl[1:])*0.5
                        hist_fpl_bd=hist_fpl/hist_bd
                        #hist_fpl_bd[np.abs(hist_fpl_bd)>1e6]=0
                        binned_fpls.append(hist_fpl_bd)
                    binned_fpls=np.array(binned_fpls)
                    binned_fpls_min=np.min(binned_fpls,axis=0)
                    binned_fpls_max=np.max(binned_fpls,axis=0)
                    #ax.plot(bins_fpl_centres,binned_fpls[0])
                    label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
                    ax.fill_between(bins_fpl_centres,binned_fpls_min,binned_fpls_max,label=label,alpha=0.35)
                    #axs[1].plot(bins_fpl_centres,hist_bd)
                    ax.set_ylim([0,1])
    axs[-1,-1].legend()
    axs[0,-1].set_xlabel(r'$q_{t=0}$')
    axs[1,-1].set_xlabel(r'$\lambda_{t=0}$')
    axs[1,0].tick_params(axis='x',which='both',labelbottom=False)
    axs[0,0].tick_params(axis='x',which='both',labelbottom=False)
    axs[1,0].tick_params(axis='y',which='both',labelleft=False)
    axs[1,1].tick_params(axis='y',which='both',labelleft=False)
    axs[0,0].set_xlim([2,4.05])
    axs[0,1].set_xlim([2,4.05])
    axs[1,0].set_xlim([0.34,.8])
    axs[1,1].set_xlim([0.34,.8])
    ax_row1.set_title('Without plasma collisions')
    ax_row2.set_title('With plasma collisions')
    ax_total.set_ylabel('Loss probability')
    #plt.savefig(
    #    fname='paper_2_loss_probability_q_lambda_bounds.pdf',
    #    dpi=800,
    #)
    plt.show()

def plot_fpl_q_lambda_1d(output_fpls,beam_deposition,batch_data,run_number=0,output_fpls_collisionless=[]):
    import matplotlib.pyplot as plt 
    import numpy as np 
    fig=plt.figure(constrained_layout=False)
    axs=[]
    ax_total=fig.add_subplot(1,1,1,frameon=False)
    ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_total.spines['top'].set_color('none')
    ax_total.spines['bottom'].set_color('none')
    ax_total.spines['left'].set_color('none')
    ax_total.spines['right'].set_color('none')
    ax_row1=fig.add_subplot(2,1,1,frameon=False)
    ax_row1.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_row1.spines['top'].set_color('none')
    ax_row1.spines['bottom'].set_color('none')
    ax_row1.spines['left'].set_color('none')
    ax_row1.spines['right'].set_color('none')
    ax_row2=fig.add_subplot(2,1,2,frameon=False)
    ax_row2.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_row2.spines['top'].set_color('none')
    ax_row2.spines['bottom'].set_color('none')
    ax_row2.spines['left'].set_color('none')
    ax_row2.spines['right'].set_color('none')
    axs.append(fig.add_subplot(2,2,1))
    axs.append(fig.add_subplot(2,2,3))
    axs.append(fig.add_subplot(2,2,2))
    axs.append(fig.add_subplot(2,2,4))
    axs=np.array(axs).reshape(2,2)
    for ax_col,output_fpl_group in zip(
        axs.T,
        [output_fpls_collisionless,output_fpls],
        ):
        for ax,quantity,grid in zip(
            ax_col,
            ['q_initial','V_pitch_initial_2D'],
            [np.linspace(2,4.05,200),np.linspace(0.34,.8,200)],
        ):
            for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower)):
                for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers)):
                    binned_fpls=[]
                    p=np.where(output_fpl_group[current_counter,mode_counter,run_number]['status_flag']=='PFC_intercept_3D')[0] 
                    hist_fpl,bins_fpl=np.histogram(
                        output_fpl_group[current_counter,mode_counter,run_number][quantity][p],
                        bins=grid,
                        weights=output_fpl_group[current_counter,mode_counter,run_number]['weight'][p],
                        )
                    hist_bd,bins_bd=np.histogram(
                        beam_deposition[quantity.replace('_initial','')],
                        bins=grid,
                        weights=beam_deposition['weight'],
                        )
                    bins_fpl_centres=(bins_fpl[:-1]+bins_fpl[1:])*0.5
                    hist_fpl_bd=hist_fpl/hist_bd
                    label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
                    ax.plot(bins_fpl_centres,hist_fpl_bd,label=label,alpha=0.35)
                    ax.set_ylim([0,1])
    axs[-1,-1].legend()
    axs[0,-1].set_xlabel(r'$q_{t=0}$')
    axs[1,-1].set_xlabel(r'$\lambda_{t=0}$')
    axs[1,0].tick_params(axis='x',which='both',labelbottom=False)
    axs[0,0].tick_params(axis='x',which='both',labelbottom=False)
    axs[1,0].tick_params(axis='y',which='both',labelleft=False)
    axs[1,1].tick_params(axis='y',which='both',labelleft=False)
    axs[0,0].set_xlim([2,4.05])
    axs[0,1].set_xlim([2,4.05])
    axs[1,0].set_xlim([0.34,.8])
    axs[1,1].set_xlim([0.34,.8])
    ax_row1.set_title('Collisionless')
    ax_row2.set_title('Collisional')
    ax_total.set_ylabel('Loss probability')
    plt.savefig(
        fname=f'paper_2_loss_probability_q_lambda_run_{run_number}.pdf',
        dpi=800,
    )
    plt.show()

import context,time,matplotlib
import matplotlib.pyplot as plt
from plot_stage_1_10_1 import * 

for key,value in batch_data.args_batch.items():
    if value:
        for val_counter,val in enumerate(value):
            try:
                batch_data.args_batch[key][val_counter]=val.replace('1_10_1','1_10_2')
            except:
                pass

input_dBs_1_10_2=np.array(list(templates.plot_mod.get_io_files(batch_data=batch_data,classtype='pert',return_lists=True)),dtype='object').reshape(
    len(batch_data.parameters__currents_upper),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers[0]),
    )

output_fpls_1_10_2=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'fpl',batch_data=batch_data,processes=4,chunksize=1)).reshape(
    len(batch_data.parameters__currents_upper),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers[0]),
    )

#calculate some data
for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower)):
    for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers)):
        for phase_counter,(parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower)): #nest at same level == offset them together rigidly 
            if output_fpls_1_10_2[current_counter,mode_counter,phase_counter] and np.all(input_dBs_1_10_2[current_counter,mode_counter,phase_counter]):
                for mode in input_dBs_1_10_2[current_counter,mode_counter,phase_counter]:
                    mode.mode_number=-int(str(mode.filename)[-1])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(theta_initial=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['R_initial'],Z=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['Z_initial'],Z_major=0.597384943))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(theta_final=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['R'],Z=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['Z'],Z_major=0.597384943))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(dTheta=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['theta_final']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['theta_initial'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(dPhi=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['phi']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['phi_initial'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['power']=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['weight']*0.5*constants.mass_deuteron*(output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_R']**2+output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_phi']**2+output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_Z']**2)
                #output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['weight']=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['weight']*0.5*constants.mass_deuteron*(output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_R']**2+output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_phi']**2+output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_Z']**2)
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(V_pitch_initial_2D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpls_1_10_2[current_counter,mode_counter,phase_counter].data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(V_pitch_initial_3D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpls_1_10_2[current_counter,mode_counter,phase_counter].data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=[mode for mode in input_dBs_1_10_2[current_counter,mode_counter,phase_counter]],i3dr=-1,phase=0.))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(V_pitch_final_2D=processing.utils.pitch_calc(output_fpls_1_10_2[current_counter,mode_counter,phase_counter],equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(V_pitch_final_3D=processing.utils.pitch_calc(output_fpls_1_10_2[current_counter,mode_counter,phase_counter],equilibria=[equilibrium],perturbations=[mode for mode in input_dBs_1_10_2[current_counter,mode_counter,phase_counter]],i3dr=-1,phase=0.))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(dV_pitch_initial=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_initial_3D']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_initial_2D'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(dV_pitch_final=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_final_3D']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_final_2D'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(dV_pitch=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_final_3D']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['V_pitch_initial_3D'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(ddV_pitch=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['dV_pitch_final']-output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['dV_pitch_initial'])
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(psi_norm_sqrt_initial=np.sqrt(output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['psi_norm_initial']))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(q_initial=processing.utils.value_at_RZ(R=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['R_initial'],Z=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['Z_initial'],quantity=equilibrium['q_rz'],grid=equilibrium))
                output_fpls_1_10_2[current_counter,mode_counter,phase_counter].set(q_final=processing.utils.value_at_RZ(R=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['R'],Z=output_fpls_1_10_2[current_counter,mode_counter,phase_counter]['Z'],quantity=equilibrium['q_rz'],grid=equilibrium))
                #to avoid discontinuities in the colour map since lots of markers near to phi=0 and theta=0
                for key in output_fpls_1_10_2[current_counter,mode_counter,phase_counter].data:
                    if np.any([angle in key for angle in ['theta','phi']]):
                        output_fpls_1_10_2[current_counter,mode_counter,phase_counter][key][output_fpls_1_10_2[current_counter,mode_counter,phase_counter][key]>np.pi]-=2.*np.pi
                        output_fpls_1_10_2[current_counter,mode_counter,phase_counter][key]*=180./np.pi
            

plot_fpl_q_lambda_bounds_1d(output_fpls,beam_deposition,batch_data,output_fpls_collisionless=output_fpls_1_10_2)


