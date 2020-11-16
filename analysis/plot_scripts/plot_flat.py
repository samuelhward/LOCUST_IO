#plot_flat.py
 
"""
Samuel Ward
30/09/2020
----
plot flat re-runs of MAST data for ASCOT/LOCUST
---
usage:
 
notes:         
---
"""

import numpy as np
import context 
import copy
from classes.output_classes.distribution_function import Distribution_Function as dfn
from classes.input_classes.equilibrium import Equilibrium as equi
from classes.output_classes.particle_list import Final_Particle_List as fpl
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.wall import Wall as wall
from run_scripts.utils import  TRANSP_output as output_T
from run_scripts.utils import TRANSP_output_FI as dfn_t
import support
import constants
import settings
import pathlib
import matplotlib.pyplot as plt
import matplotlib
import processing.utils

#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_g_=settings.colour_custom([205,220,57,1])
cmap_b=settings.colour_custom([33,150,243,1])
cmap_grey=settings.colour_custom([97,97,97,1])

fig,ax=plt.subplots(1)

parent_dir='MAST_reruns_flat'
codes=['ASCOT','LOCUST']
tracking_types=['FO','GC']
axes=['E']

for code in codes:
    for tracking_type in tracking_types:
        glob='.h5' if code is 'ASCOT' else '.dfn'
        DFN=dfn(ID=f'{code} {tracking_type}',data_format=code,filename=list((support.dir_output_files / parent_dir / code / 'outputs' / tracking_type).glob(f'*{glob}'))[0])        
        colmap=cmap_b if code is 'ASCOT' else cmap_g

        line_style = 'dashed' if tracking_type is 'GC' else 'solid'

        if code is 'ASCOT':
            beam_power=1. #desired beam power
            ASCOT_beam_depo=bd(ID='ASCOT beam depo',data_format='ASCOT_FO',filename=pathlib.Path(parent_dir) / 'ASCOT' / 'inputs' / '29034W03' / 'input.particles')
            if 'E' not in ASCOT_beam_depo.data: 
                ASCOT_beam_depo['E']=0.5*constants.mass_deuteron*(ASCOT_beam_depo['V_R']**2+ASCOT_beam_depo['V_phi']**2+ASCOT_beam_depo['V_Z']**2)/constants.species_charge
            Pdep_ASCOT=np.sum(ASCOT_beam_depo['E']*constants.species_charge*ASCOT_beam_depo['weight'])
            DFN['dfn']*=beam_power/Pdep_ASCOT 

        DFN.plot(axes=axes,line_style=line_style,ax=ax,fig=fig,colmap=colmap)


plt.show()

#################################

##################################################################

###################################################################################################
