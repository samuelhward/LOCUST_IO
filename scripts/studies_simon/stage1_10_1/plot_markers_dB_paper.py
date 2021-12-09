import context,time,matplotlib
import matplotlib.pyplot as plt
from plot_stage_1_10_1 import * 

fig=plt.figure(constrained_layout=True)
axs=[]
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
axs.append(fig.add_subplot(2,1,1))
axs.append(fig.add_subplot(2,1,2))
axs=np.array(axs).flatten()

lims=[[[-180,0],[-120,140]],  
      [[-94.,-13],[-116,-89.]]]
for ax,lim in zip(axs,lims):
    mesh=templates.plot_mod.plot_perturbation_phi_theta(perturbations=input_dBs[0,1,2],equilibrium=equilibrium,psi=0.99,vminmax=None,number_bins=20,colmap=matplotlib.cm.get_cmap('Greys'),ax=ax,fig=fig)
    #fig.colorbar(mesh,ax=ax,orientation='vertical')
    mesh=output_fpls[0,0,2].plot(axes=['phi','theta_final'],style='scatter',colfield='dV_pitch_final',colmap=matplotlib.cm.get_cmap('plasma'),ax=ax,fig=fig)
    ax.set_xlim(lim[0])
    ax.set_ylim(lim[1])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')

ax_total.set_xlabel(r'Toroidal angle $\phi$ [deg]')
ax_total.set_ylabel(r'Poloidal angle $\theta$ [deg]',labelpad=15)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar=fig.colorbar(mesh,cax=cbar_ax,orientation='vertical')
#cbar.set_label(r'Pitch $v_{\parallel}/v$') #if instead plotting V_pitch_final_2D above instead
cbar_ax.set_title(r'$\widetilde{\lambda}-\lambda$')
plt.show()


'''
fig,ax=plt.subplots(1)
for output_fpl in output_fpls_1_10_2[0,0,:]:
    p=np.where(output_fpl['status_flag']=='PFC_intercept_3D')[0] 
    hist_fpl,bins_fpl=np.histogram(
        output_fpl['V_pitch_initial_2D'][p],
        bins=np.linspace(0,1,300),
        weights=output_fpl['weight'][p],
        )
    hist_bd,bins_bd=np.histogram(
        beam_deposition['V_pitch_initial_2D'.replace('_initial','')],
        bins=np.linspace(0,1,300),
        weights=beam_deposition['weight'],
        )
    bins_fpl_centres=(bins_fpl[:-1]+bins_fpl[1:])*0.5
    hist_fpl_bd=hist_fpl/hist_bd
    ax.plot(bins_fpl_centres,hist_fpl_bd)

plt.show()


inds=(0,1,0)
for output_fpl in [output_fpls_1_10_2[inds]]: 
    axes=['q_initial','V_pitch_initial_2D']
    grid={}
    grid['R_initial']=np.linspace(8,9.5,400)
    grid['q_initial']=np.linspace(2,4.5,400)
    grid['V_pitch_initial_2D']=np.linspace(0,1.,400)
    p=np.where((output_fpl['status_flag']=='PFC_intercept_3D') & (np.abs(output_fpl['theta_final'])<180))[0] 
    output_fpl_binned,output_fpl_binned_x,output_fpl_binned_y=np.histogram2d(output_fpl[axes[0]][p],output_fpl[axes[1]][p],bins=[grid[axes[0]],grid[axes[1]]],weights=output_fpl['weight'][p])
    bd_binned,bd_binned_x,bd_binned_y=np.histogram2d(beam_deposition[axes[0].replace('_initial','')],beam_deposition[axes[1].replace('_initial','')],bins=[grid[axes[0]],grid[axes[1]]],weights=beam_deposition['weight'])
    #output_fpl_binned_x and output_fpl_binned_x are first edges then converted to centres
    output_fpl_binned_x=(output_fpl_binned_x[:-1]+output_fpl_binned_x[1:])*0.5
    output_fpl_binned_y=(output_fpl_binned_y[:-1]+output_fpl_binned_y[1:])*0.5
    dx,dy=output_fpl_binned_x[1]-output_fpl_binned_x[0],output_fpl_binned_y[1]-output_fpl_binned_y[0]                        
    output_fpl_binned_y,output_fpl_binned_x=np.meshgrid(output_fpl_binned_y-dy/2.,output_fpl_binned_x-dx/2.) #offset ticks onto bin centres
    fig,ax=plt.subplots(1)
    mesh=ax.pcolormesh(output_fpl_binned_x,output_fpl_binned_y,output_fpl_binned/bd_binned,cmap=settings.cmap_inferno)
    plt.show()
'''