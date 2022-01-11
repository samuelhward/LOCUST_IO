#plot_scenario_data.py
 
"""
Samuel Ward
23/07/21
----
script for plotting stage 1.10.1 of Simon's studies 
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
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
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

import plot_scenario_data_launch as batch_data

beam_depositions=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'beamdepo',batch_data,processes=16,chunksize=1)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))
equilibria=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'eq',batch_data,processes=16,chunksize=1,GEQDSKFIX1=True,GEQDSKFIX2=True)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))
temperatures_e=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'temp_e',batch_data,processes=16,chunksize=1)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))
temperatures_i=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'temp_i',batch_data,processes=16,chunksize=1)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))
density_electrons=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'numden',batch_data,processes=16,chunksize=1)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))
rotations=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'rot',batch_data,processes=16,chunksize=1)).reshape(len(batch_data.configs_beam_1),len(batch_data.parameters__databases))


if __name__=='__main__':

    #'''
    #paper 2 figures

    #   KINETIC PROFILES FOR Q=10 SCENARIO

    fig1,ax1=plt.subplots(1)
    fig1.subplots_adjust(right=0.7)
   
    p1,=ax1.plot(
        equilibria[0,3]['flux_pol_norm'],
        equilibria[0,3]['qpsi'],   
        label='',
        color=settings.cmap_k(0.),
    )
    ax2 = ax1.twinx()
    ax2.spines.right.set_position(("axes", 1.))
    p2,=ax2.plot(
        density_electrons[0,3]['flux_pol_norm'],
        density_electrons[0,3]['n'],    
        label='',
        color=cmap_r(0.),
    )
    ax3 = ax1.twinx()
    ax3.spines.right.set_position(("axes", 1.1))
    p3,=ax3.plot(
        temperatures_e[0,3]['flux_pol_norm'],
        temperatures_e[0,3]['T']/1.e3,
        label='',
        color=cmap_g(0.),
    )

    p3_5,=ax3.plot(
        temperatures_i[0,3]['flux_pol_norm'],
        temperatures_i[0,3]['T']/1.e3,
        label='',
        color=cmap_g(0.),
        linestyle='dashed',
    )

    ax4 = ax1.twinx()
    ax4.spines.right.set_position(("axes", 1.2))
    p4,=ax4.plot(
        rotations[0,3]['flux_pol_norm'],
        np.abs(rotations[0,3]['rotation_ang'])/1000.,    
        label='',
        color=cmap_b(0.),
    )

    ax1.set_ylabel(r'Safety factor')
    ax2.set_ylabel(r'Electron number density $n_{\mathrm{e}}$ [$\mathrm{m}^{-3}$]')
    ax3.set_ylabel(r'Temperature $T_{\mathrm{e}}$, $T_{\mathrm{i}}$ (dashed) [keV]') 
    ax4.set_ylabel(r'Bulk rotation velocity $\omega_{\phi}$ [krad/s]')

    ax1.yaxis.label.set_color(p1.get_color())#settings.cmap_b(0.))
    ax2.yaxis.label.set_color(p2.get_color())#settings.cmap_r(0.))
    ax3.yaxis.label.set_color(p3.get_color())#settings.cmap_g(0.))
    ax4.yaxis.label.set_color(p4.get_color())#settings.cmap_k(0.))

    ax1.tick_params(axis='y', colors=p1.get_color(),)
    ax2.tick_params(axis='y', colors=p2.get_color(),)
    ax3.tick_params(axis='y', colors=p3.get_color(),)
    ax4.tick_params(axis='y', colors=p4.get_color(),)

    ax1.set_xlabel(r'Normalised poloidal flux $\psi_{\mathrm{n}}$')

    plt.show()

    #'''
    #paper 3 figures


    #   KINETIC PROFILES FOR ALL SCENARIOS

    fig=plt.figure(constrained_layout=True)

    ax_total=fig.add_subplot(1,1,1,frameon=False)
    ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    ax_total.spines['top'].set_color('none')
    ax_total.spines['bottom'].set_color('none')
    ax_total.spines['left'].set_color('none')
    ax_total.spines['right'].set_color('none')
    ax_total.set_xlabel(r'$\psi_{\mathrm{n}}$',labelpad=10)

    ax=[]
    ax.append([fig.add_subplot(2,2,counter+1) for counter in range(4)])
    ax=np.array(ax).flatten()

    colours=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for counter,scenario in enumerate(batch_data.parameters__databases):

        ax[0].plot(
            equilibria[0,counter]['flux_pol_norm'],
            equilibria[0,counter]['qpsi'],   
            label=f'{scenario}'.split('_')[-1],
            color=colours[2*counter+2],
        )
        ax[1].plot(
            temperatures_e[0,counter]['flux_pol_norm'],
            temperatures_e[0,counter]['T']/1.e3,
            label=f'{scenario}'.split('_')[-1],
            color=colours[2*counter+2],
        )
        ax[1].plot(
            temperatures_i[0,counter]['flux_pol_norm'],
            temperatures_i[0,counter]['T']/1.e3,
            label=f'{scenario}'.split('_')[-1],
            color=colours[2*counter+2],
            linestyle='dashed',
        )
        ax[2].plot(
            density_electrons[0,counter]['flux_pol_norm'],
            density_electrons[0,counter]['n'],    
            label=f'{scenario}'.split('_')[-1],
            color=colours[2*counter+2],
        )
        ax[3].plot(
            rotations[0,counter]['flux_pol_norm'],
            np.abs(rotations[0,counter]['rotation_ang'])/1000.,    
            label=f'{scenario}'.split('_')[-1],
            color=colours[2*counter+2],
        )

    ax[0].tick_params(
        axis='x',
        which='both',
        bottom=True,
        top=False,
        labelbottom=False
        )
    ax[1].tick_params(
        axis='x',
        which='both',
        bottom=True,
        top=False,
        labelbottom=True
        )
    ax[2].tick_params(
        axis='x',
        which='both',
        bottom=True,
        top=False,
        labelbottom=True
        )

    ax[0].set_ylabel(r'$q$')
    ax[1].set_ylabel(r'$T_{\mathrm{e}}$, $T_{\mathrm{i}}$ (dashed) [keV]') 
    ax[2].set_ylabel(r'$n_{\mathrm{e}}$ [$\mathrm{m}^{-3}$]')
    ax[3].set_ylabel(r'$\omega_{\phi}$ [krad/s]')

    ax[0].legend()

    plt.show()

    # beam deposition and NBI geometry for Q=10 scenario

    combined=Beam_Deposition('')

    for depo,mass in zip(beam_depositions[0:2,-1],[constants.mass_deuteron/2.,constants.mass_deuteron]):
        depo['weight']*=depo['E']*mass*constants.species_charge/constants.species_mass/1.e6 #deposition in MW

    combined.combine(*[depo for depo in beam_depositions[0:2,-1]])
    fig=plt.figure(constrained_layout=False)
    ax1 = fig.add_subplot(221)
    axes=['R','Z']
    run_scripts.utils.plot_NBI_geometry(axes=axes,beam_name='heating beam 1',fig=fig,ax=ax1,axis='off',colmap=settings.cmap_r)
    run_scripts.utils.plot_NBI_geometry(axes=axes,beam_name='heating beam 2',fig=fig,ax=ax1,axis='on',colmap=cmap_b)
    run_scripts.utils.plot_NBI_geometry(axes=axes,beam_name='diagnostic',fig=fig,ax=ax1,colmap=cmap_g)
    combined.plot(axes=axes,number_bins=500,fig=fig,ax=ax1,real_scale=True,LCFS=equilibria[0,-1],colmap=settings.discrete_colmap(colmap_name='plasma',face_colour='white',number_bins=50))
    ax1.set_xlabel('$R$ [m]')
    ax1.set_ylabel('$Z$ [m]')
    ax1.set_xlim((4.,8.5))
    ax1.set_xticks((4,6,8))
    ax2 = fig.add_subplot(222,polar=True)
    axes=['phi','R']
    run_scripts.utils.plot_NBI_geometry(axes=axes,label='heating beam 1',beam_name='heating beam 1',fig=fig,ax=ax2,axis='off',colmap=settings.cmap_r)
    run_scripts.utils.plot_NBI_geometry(axes=axes,label='heating beam 2',beam_name='heating beam 2',fig=fig,ax=ax2,axis='on',colmap=cmap_b)
    run_scripts.utils.plot_NBI_geometry(axes=axes,label='diagnostic',beam_name='diagnostic',fig=fig,ax=ax2,colmap=cmap_g)
    combined.plot(axes=axes,number_bins=500,fig=fig,ax=ax2,real_scale=True,LCFS=equilibria[0,-1],colmap=settings.discrete_colmap(colmap_name='plasma',face_colour='white',number_bins=50))
    ax2.set_xticks(np.pi/180. * np.linspace(0, 90, 3, endpoint=False))
    ax2.set_yticks((0.,3.,6.,9))
    ax2.set_ylim((0,9.))
    ax2.set_thetamin(0)
    ax2.set_thetamax(90)
    ax2.set_xlabel('$X$ [m]')
    ax2.set_ylabel('$Y$')

    ax3 = fig.add_subplot(234)
    ax4 = fig.add_subplot(235)
    ax5 = fig.add_subplot(236)
    ax4.tick_params(
        axis='y',
        which='both',
        left=False,
        bottom=True,
        top=False,
        right=False,
        labelleft=False,
        labelbottom=True
        )

    #XXX be careful with legend here, fudged/overplotted to add new lines to plot
    colours=[cmap_g,settings.cmap_r]
    for depo,colour,mass,label in zip(beam_depositions[0:2,-1],colours,[constants.mass_deuteron/2.,constants.mass_deuteron],['DNB','HNB1']):
        depo_binned_Z,depo_binned_Z_edges=np.histogram(depo['Z'],bins=np.linspace(-1.5,1.5,100),weights=depo['weight'])
        depo_binned_Z_centres=(depo_binned_Z_edges[:-1]+depo_binned_Z_edges[1:])*0.5        
        ax4.plot(depo_binned_Z_centres,depo_binned_Z/(depo_binned_Z_centres[1]-depo_binned_Z_centres[0]),color=colour(0.),label=label)
        depo_binned_R,depo_binned_R_edges=np.histogram(depo['R'],bins=np.linspace(4,8.5,100),weights=depo['weight'])
        depo_binned_R_centres=(depo_binned_R_edges[:-1]+depo_binned_R_edges[1:])*0.5
        ax3.plot(depo_binned_R_centres,depo_binned_R/(depo_binned_R_centres[1]-depo_binned_R_centres[0]),color=colour(0.),label=label)
        depo.set(V_pitch=processing.utils.pitch_calc(depo,equilibria=[equilibria[0,-1]],perturbations=None,i3dr=-1,phase=0.))
        depo_binned_pitch,depo_binned_pitch_edges=np.histogram(depo['V_pitch'],bins=np.linspace(-1.,1.,100),weights=depo['weight'])
        depo_binned_pitch_centres=(depo_binned_pitch_edges[:-1]+depo_binned_pitch_edges[1:])*0.5
        ax5.plot(depo_binned_pitch_centres,depo_binned_pitch/(depo_binned_pitch_centres[1]-depo_binned_pitch_centres[0]),color=colour(0.))

    ax3.plot(depo_binned_R_centres,depo_binned_R/(depo_binned_R_centres[1]-depo_binned_R_centres[0]),color=cmap_b(0.),label='HNB2')
    ax3.plot(depo_binned_R_centres,depo_binned_R/(depo_binned_R_centres[1]-depo_binned_R_centres[0]),color=settings.cmap_k(0.),label='HNB1 + HNB2')
    ax4.plot(depo_binned_Z_centres,depo_binned_Z/(depo_binned_Z_centres[1]-depo_binned_Z_centres[0]),color=settings.cmap_k(0.),label='HNB1 + HNB2')
    ax5.plot(depo_binned_pitch_centres,depo_binned_pitch/(depo_binned_pitch_centres[1]-depo_binned_pitch_centres[0]),color=settings.cmap_k(0.),label='HNB1 + HNB2')

    ax3.set_title('')
    ax4.set_title('')
    ax3.set_xlabel(r'$R$ [m]')
    ax4.set_xlabel(r'$Z$ [m]')
    ax5.set_xlabel(r'Pitch ($v_{\parallel}/v$)')
    ax4.set_ylabel(r'')
    ax3.set_ylabel(r'Power deposited [MWm$^{-1}$]')
    ax5.text(0.05, 0.95, r'[MW$\lambda^{-1}$]', ha='left', va='top', rotation='horizontal',transform=ax5.transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.)) 
    ax3.set_xlim((4.,8.5))
    ax4.set_xlim((-1.5,1.5))
    ax5.set_xlim((-1.,1.))
    ax3.set_ylim((0,35.))
    ax4.set_ylim((0,35.))
    ax3.axvline(np.min(equilibria[0,-1]['lcfs_r']),color='m',label='LCFS')
    ax3.axvline(np.max(equilibria[0,-1]['lcfs_r']),color='m')
    ax4.axvline(np.min(equilibria[0,-1]['lcfs_z']),color='m')
    ax4.axvline(np.max(equilibria[0,-1]['lcfs_z']),color='m')
    ax3.legend()

    plt.show()





