#blip test plot

#reads and plots the DFNs produced by blip test

import context
from classes.output_classes.distribution_function import Distribution_Function
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.output_classes.particle_list import Final_Particle_List
from classes.output_classes.moments import Moments as Moments
import processing.process_output
import processing.utils 
import run_scripts.utils
import constants
import numpy as np
import matplotlib.pyplot as plt
import copy
import matplotlib
from matplotlib import cm

cmap_ascot=matplotlib.cm.get_cmap('Reds') #ASCOT colourmap
cmap_locust=matplotlib.cm.get_cmap('Greens') #LOCUST colourmap

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq,GEQDSKFIX=0)

#ascot='ascot/GC/proper/ascot_freia_9860411.h5' #proper run
#locust='locust/proper/F_03-12-2018_22-42-42.478_TOTL.dfn' #proper run

#ascot='ascot/GC/flat_prof/ascot_freia_172004.h5' #flat profiles
#ascot='ascot/GC/flat_prof_fix_timestep/ascot_freia_199662.h5' #flat profiles with fixed timestep
#ascot='ascot/GC/flat_prof_no_acc/ascot_freia_691650.h5' #flat profiles with no acceleration
#ascot='ascot/GC/flat_prof_shifted/ascot_freia_1104091.h5' #shifted out the birth list by 0.1m
#ascot='ascot/GC/flat_prof_shifted/ascot_freia_1459286.h5'#shifted still but now colliding against a D background not H
#ascot='ascot/GC/flat_prof_shifted_more/ascot_freia_1104172.h5' #shifted out the birth list by 0.3m
#ascot='ascot/GC/flat_prof_shifted_more/ascot_freia_1453009.h5' #0.3m shift but now using D background instead of H 
#ascot='ascot/GC/flat_prof_shifted_more_high_stats/ascot_freia_1777904.h5' #0.3m shift but now using D background instead of H 
#ascot='ascot/GC/flat_prof_shifted_high_stats/ascot_freia_1700856.h5' #same as flat_prof_shifted but increased the number of particles
#ascot='ascot/FO/flat_prof_shifted/ascot_freia_1452884.h5' #full orbit with D background
ascot='ascot/GC/flat_prof_shifted_no_ac/ascot_freia_1906880.h5' #no acceleration, 10cm shift and D background
#ascot='ascot/GC/flat_prof_shifted_high_stats_no_ac/ascot_freia_1906881.h5' #high stats no acceleration, 10cm shift and D background
#locust='locust/GC/flat_prof/F_17-12-2018_16-03-01.537_TOTL.dfn' #flat profiles
#locust='locust/GC/flat_prof_shifted/F_08-02-2019_14-55-28.333_TOTL.dfn' #shifted out the birth list by 0.1m
#locust='locust/GC/flat_prof_shifted_more/F_21-03-2019_09-32-01.051_TOTL.dfn' #shifted the birth list out by 0.3m
#locust='locust/GC/flat_prof_shifted_high_stats/F_19-04-2019_12-17-32.730_TOTL.dfn' #same as flat_prof_shifted but increased the number of particles
locust='locust/GC/flat_prof_shifted_highunbor/F_08-02-2019_19-05-21.140_TOTL.dfn' #increased UNBOR to 100
#locust='locust/FO/flat_prof_shifted/F_12-02-2019_12-07-12.890_TOTL.dfn'  
#locust='locust/FO/flat_prof_no_ions/F_19-05-2019_18-41-36.857_TOTL.dfn'

ascot_dfn=Distribution_Function('ASCOT 10cm shift','ASCOT',ascot)
locust_dfn=Distribution_Function('LOCUST 10cm shift','LOCUST',locust)

ascot_list=Final_Particle_List('ascot','ASCOT',ascot)

#need corresponding beam depo to calculate the weight scaling - MAKE SURE THIS IS UPDATED
beam_depo=Beam_Deposition('ASCOT beam depo','ASCOT_GC','input.particles')
beam_power=1. #1 Watt beam power
if 'E' not in beam_depo.data.keys():
     beam_depo['E']=0.5*constants.species_mass*(beam_depo['V_tor']**2+beam_depo['V_R']**2+beam_depo['V_Z']**2)/constants.e_charge
BPCAP=np.sum(beam_depo['E']*constants.species_charge*beam_depo['weight'])
ascot_dfn['dfn']*=beam_power/BPCAP
ascot_dfn['V_pitch']*=-1.
 
print(BPCAP)

#locust_dfn['dfn']*=0.7935379718996407 #random factor I'm out by



fig,ax=plt.subplots(1)
quantity='density'
ascot_mom=Moments(ascot_dfn.ID,'ASCOT',ascot)
locust_mom=Moments(locust_dfn.ID,'LOCUST','locust/GC/flat_prof_shifted_high_stats/LOCUST_19-04-2019_12-17-32.730.h5')
ascot_mom[quantity]*=1e6*beam_power/BPCAP
ascot_mom.plot(key=quantity,fig=fig,ax=ax,colmap='r')
locust_mom.plot(key=quantity,fig=fig,ax=ax,colmap='g')
ax.legend((ascot_mom.ID,locust_mom.ID))
plt.show()

k=np.where(ascot_list['status_flag']!=4)[0] #check if any particles were rejected or hit a PFC
print(k)


fig,(ax1,ax2,ax3)=plt.subplots(1,3)

axes=['E']

vminmax=[0.,1.7e9]
#ascot_dfn=ascot_dfn.crop(E=[76000,83000])
#locust_dfn=locust_dfn.crop(E=[76000,83000])
ascot_dfn.plot(axes=axes,ax=ax3,fig=fig,colmap='r',vminmax=vminmax)
locust_dfn.plot(axes=axes,ax=ax3,fig=fig,colmap='g',vminmax=vminmax)
ax3.axvline(80000,color='m')
ax3.legend((ascot_dfn.ID,locust_dfn.ID))
ax3.set_title('')
ax3.set_ylabel('# per eV')

#'''
axes=['E','V_pitch']
mesh_ascot=ascot_dfn.plot(axes=axes,ax=ax1,fig=fig,fill=False,number_bins=30)
mesh_locust=locust_dfn.plot(axes=axes,ax=ax2,fig=fig,fill=False,number_bins=30)
for ax in [ax1,ax2]:
    ax.set_ylim([-1.0,-0.7])
    ax.set_xlim([0.,100000.])
#'''

'''
axes=['R','Z']
vminmax=[0,2.5e13]
mesh_ascot=ascot_dfn.plot(axes=axes,ax=ax1,fig=fig,LCFS=equi,fill=True,real_scale=True,vminmax=vminmax)
mesh_locust=locust_dfn.plot(axes=axes,ax=ax2,fig=fig,LCFS=equi,fill=True,real_scale=True,vminmax=vminmax)
for ax in [ax1,ax2]:
   ax.set_ylim([-0.9,0.9])
   ax.set_xlim([1.25,2.15])
'''

colourbars=True
colourbar_array=[]
if colourbars is True:
    for ax,mesh in zip([ax1,ax2],[mesh_ascot,mesh_locust]):
        cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
        colourbar_array.append(cbar)

plt.show()

axes=['V_pitch']
ascot_dfn_=processing.process_output.dfn_crop(ascot_dfn,V_pitch=[-1.,1])
locust_dfn_=processing.process_output.dfn_crop(locust_dfn,V_pitch=[-1.,1]) 
ascot_dfn_=ascot_dfn_.transform(axes=axes)
locust_dfn_=locust_dfn.transform(axes=axes)
plt.plot(ascot_dfn_[axes[0]],ascot_dfn_['dfn'],'r')
plt.plot(locust_dfn_[axes[0]],locust_dfn_['dfn'],'g')
plt.legend((ascot_dfn.ID,locust_dfn.ID))
plt.xlabel(axes[0])
plt.show()








'''


points={}
points['E']=[80000]
points['V_pitch']=[-0.961]
points['R']=[1.71]
points['Z']=[0.0251]
ascot_points=processing.utils.get_dfn_point(ascot_dfn,**points)
points['P']=[0] #ASCOT does not have this dimension
locust_points=processing.utils.get_dfn_point(locust_dfn,**points)


ascot_dfn__=copy.deepcopy(ascot_dfn)
locust_dfn__=copy.deepcopy(locust_dfn)

for counter,E in enumerate(ascot_dfn__['E']):
    ascot_dfn__['dfn'][counter,:,:,:]*=E

for counter,E in enumerate(locust_dfn__['E']):
    locust_dfn__['dfn'][:,counter,:,:,:]*=E

P_ascot=ascot_dfn__.dfn_transform(axes=['N'])
P_ascot=P_ascot['dfn']
P_locust=processing.process_output.dfn_transform(locust_dfn__,axes=['N'])
P_locust=P_locust['dfn']

'''
