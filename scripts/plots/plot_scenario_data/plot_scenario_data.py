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

beam_depositions=np.array(list(templates.plot_mod.get_output_files(batch_data,'beamdepo'))).reshape(len(batch_data.parameters__databases),len(batch_data.configs_beam_1))

equilibria=np.array(list(templates.plot_mod.get_output_files(batch_data,'eq',GEQDSKFIX1=True,GEQDSKFIX2=True))).reshape(len(batch_data.parameters__databases),len(batch_data.configs_beam_1))

temperatures=np.array(list(templates.plot_mod.get_output_files(batch_data,'temp_e',species='electrons'))).reshape(len(batch_data.parameters__databases),len(batch_data.configs_beam_1))

density_electrons=np.array(list(templates.plot_mod.get_output_files(batch_data,'numden',species='electrons'))).reshape(len(batch_data.parameters__databases),len(batch_data.configs_beam_1))

rotations=np.array(list(templates.plot_mod.get_output_files(batch_data,'rot'))).reshape(len(batch_data.parameters__databases),len(batch_data.configs_beam_1))


if __name__=='__main__':


    #paper 2 figures


    #   KINETIC PROFILES FOR Q=10 SCENARIO

    fig1,ax1=plt.subplots(1)

    p1=ax1.plot(
        equilibria[2,0]['flux_pol_norm'],
        equilibria[2,0]['qpsi'],   
        label='',
        color=cmap_b(0.),
    )
    ax2 = ax1.twinx()
    p2=ax2.plot(
        temperatures[2,0]['flux_pol_norm'],
        temperatures[2,0]['T'],
        label='',
        color=cmap_r(0.),
    )
    ax3 = ax1.twinx()
    p3=ax3.plot(
        density_electrons[2,0]['flux_pol_norm'],
        density_electrons[2,0]['n'],    
        label='',
        color=cmap_g(0.),
    )
    ax4 = ax1.twinx()
    p4=ax4.plot(
        rotations[2,0]['flux_pol_norm'],
        rotations[2,0]['rotation_ang'],    
        label='',
        color=settings.cmap_k(0.),
    )

    ax1.set_ylabel(r'$q$')
    ax2.set_ylabel(r'$T_{\mathrm{e}}$ [keV]') 
    ax3.set_ylabel(r'$n_{\mathrm{e}}$ [$\mathrm{m}^{-3}$]')
    ax4.set_ylabel(r'$\omega_{\phi}$ [rad/s]')

    ax1.yaxis.label.set_color(p1.get_color())
    ax2.yaxis.label.set_color(p2.get_color())
    ax3.yaxis.label.set_color(p3.get_color())
    ax4.yaxis.label.set_color(p4.get_color())

    ax1.set_xlabel(r'$\psi_{\mathrm{n}}$')

    plt.show()


    #paper 3 figures


    #   KINETIC PROFILES FOR ALL SCENARIOS

    fig,ax=plt.subplots(4)

    for counter,scenario in enumerate(batch_data.parameters__databases):

        ax[0].plot(
            equilibria[2,0]['flux_pol_norm'],
            equilibria[2,0]['qpsi'],   
            label=f'{scenario}',
            color=settings.cmap_inferno(counter/len(batch_data.parameters__databases)),
        )
        ax[1].plot(
            temperatures[2,0]['flux_pol_norm'],
            temperatures[2,0]['T'],
            label=f'{scenario}',
            color=settings.cmap_inferno(counter/len(batch_data.parameters__databases)),
        )
        ax[2].plot(
            density_electrons[2,0]['flux_pol_norm'],
            density_electrons[2,0]['n'],    
            label=f'{scenario}',
            color=settings.cmap_inferno(counter/len(batch_data.parameters__databases)),
        )
        ax[3].plot(
            density_electrons[2,0]['flux_pol_norm'],
            density_electrons[2,0]['rotation_ang'],    
            label=f'{scenario}',
            color=settings.cmap_inferno(counter/len(batch_data.parameters__databases)),
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
    ax[1].set_ylabel(r'$T_{\mathrm{e}}$ [keV]') 
    ax[2].set_ylabel(r'$n_{\mathrm{e}}$ [$\mathrm{m}^{-3}$]')
    ax[3].set_ylabel(r'$\omega_{\phi}$ [rad/s]')

    ax[-1].set_xlabel(r'$\psi_{\mathrm{n}}$')

    plt.show()

    #   phi R for Q=10 scenario

    fig=plt.figure()
    ax = fig.add_subplot(111,polar=True)
    combined=beam_depositions[2,0].combine(*[depo for depo in beam_depositions[2,1:]])
    combined.plot(axes=['phi','R'],number_bins=500,fig=fig,ax=ax)