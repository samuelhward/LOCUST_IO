#look at RZ, EP and EP sampled plots for TRANSP vs LOCUST vs ASCOT then animate through the limiter radii

import sys
import numpy as np
import context
import copy
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.wall import Wall
from classes.output_classes.distribution_function import Distribution_Function
from classes.output_classes.particle_list import Final_Particle_List
from processing import process_input
from processing import process_output
import run_scripts.utils
import constants
import settings


filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)



wall_files='input.wall_2d_'
radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
run_IDs=['U69','U70','U71','U72','U73','U74']
shot_number='157418'

colours=['r-','g-','b-','m-','k-','c-']

TRANSP_files_tail_FI='_fi_1_gc.cdf'
TRANSP_files_tail_CDF='.CDF'

ASCOT_files=['ascot_freia_1470025.h5','ascot_freia_1470029.h5','ascot_freia_1470032.h5','ascot_freia_1470036.h5','ascot_freia_1470040.h5','ascot_freia_1470044.h5']
ASCOT_run='ascot/run_1/' #this is with old kinetic profiles which are not extrapolated, ORBITMETHOD=1

ASCOT_files=['ascot_freia_1470026.h5','ascot_freia_1470030.h5','ascot_freia_1470033.h5','ascot_freia_1470037.h5','ascot_freia_1470041.h5','ascot_freia_1470045.h5']
ASCOT_run='ascot/run_2/' #changed ORBITMETHOD to 4, added extrapolated kinetic profiles

ASCOT_files=['ascot_freia_1470027.h5','ascot_freia_1470031.h5','ascot_freia_1470034.h5','ascot_freia_1470038.h5','ascot_freia_1470042.h5','ascot_freia_1470046.h5']
ASCOT_run='ascot/run_3/' #changed ORBITMETHOD back to 1, keep extrapolated kinetic profiles

ASCOT_files=['ascot_freia_1470028.h5','ascot_freia_1480719.h5','ascot_freia_1470035.h5','ascot_freia_1470039.h5','ascot_freia_1470043.h5','ascot_freia_1470047.h5']
ASCOT_run='ascot/run_4/' #ORBITMETHOD 4 and non-extrapolated kinetic profiles

LOCUST_beam_depo_tail='_ptcles.dat'
LOCUST_files=['F_04-12-2018_16-11-28.285_TOTL.dfn','F_03-12-2018_22-55-59.961_TOTL.dfn','F_03-12-2018_22-58-49.281_TOTL.dfn','F_03-12-2018_22-57-37.348_TOTL.dfn','F_04-12-2018_00-13-14.371_TOTL.dfn','F_04-12-2018_14-11-36.625_TOTL.dfn']
LOCUST_moments=['LOCUST_04-12-2018_16-11-28.285.h5','LOCUST_03-12-2018_22-55-59.961.h5','LOCUST_03-12-2018_22-58-49.281.h5','LOCUST_03-12-2018_22-57-37.348.h5','LOCUST_04-12-2018_00-13-14.371.h5','LOCUST_04-12-2018_14-11-36.625.h5']
LOCUST_run='locust/run_1/'
LOCUST_files=['F_04-12-2018_16-10-58.977_TOTL.dfn','F_04-12-2018_15-06-13.311_TOTL.dfn','F_04-12-2018_15-10-53.134_TOTL.dfn','F_04-12-2018_17-17-41.999_TOTL.dfn','F_04-12-2018_17-19-35.424_TOTL.dfn','F_05-12-2018_01-15-34.057_TOTL.dfn']
LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
LOCUST_run='locust/run_2/' #added -DSOLCOL

#plot just one radius (still need to comment out desired files where necessary)
rad=0
radii=[radii[rad]]
ASCOT_files=[ASCOT_files[rad]]
run_IDs=[run_IDs[rad]]
LOCUST_files=[LOCUST_files[rad]]

#start by creating some axis objects
fig,((ax1,ax2,ax3))=plt.subplots(1,3)

for radius,LOCUST_file,ASCOT_file,run_ID,colour in zip(radii,LOCUST_files,ASCOT_files,run_IDs,colours):

    #for ax in [ax1,ax2,ax3]:
    #    ax.cla()

    #toggle colourbars
    colourbars=True
    colourbar_array=[]

    LOCUST_dfn=Distribution_Function('LOCUST density (0.1s) - r = {}'.format(str(radius)),data_format='LOCUST',filename=LOCUST_run+LOCUST_file)
    TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP density (0.1s) - r = {}'.format(str(radius)),filename=shot_number+run_ID+TRANSP_files_tail_FI)
    ASCOT_dfn=Distribution_Function('ASCOT density (0.1s) - r = {}'.format(str(radius)),data_format='ASCOT',filename=ASCOT_run+ASCOT_file,datatype='distribution_function')

    #import wall, redefine pitch for ASCOT vs current in DIII-D and normalise DFNs according to beam powers
    wall=Wall('limiter - '+radius,data_format='ASCOT_2D_input',filename=wall_files+radius)
    ASCOT_dfn['V_pitch']*=-1.

    TRANSP_CDF=run_scripts.utils.TRANSP_output(ID='{}'.format(shot_number+run_ID+TRANSP_files_tail_CDF),filename=shot_number+run_ID+TRANSP_files_tail_CDF) #need to scale some things by captured beam power
    output_time=3.1 #time at which the distribution function was written out
    BPCAP_index=np.abs(TRANSP_CDF['TIME']-output_time).argmin()
    BPCAP=TRANSP_CDF['BPCAP'][BPCAP_index]
    TRANSP_dfn['dfn']/=BPCAP

    beam_depo=Final_Particle_List(ID=run_ID,data_format='ASCOT',filename=ASCOT_run+ASCOT_file)
    beam_power=1. #1 Watt beam power
    beam_energy=80000
    k=np.where(beam_depo['status_flag']>0)[0] #scale by all markers which actually contributed to the simulation
    BPCAP_ascot=np.sum(beam_energy*constants.species_charge*beam_depo['weight'][k])
    ASCOT_dfn['dfn']*=beam_power/(BPCAP_ascot)

    #RZ

    axes=['R','Z']

    vminmax=[0.,2.e10]

    #####optional crop of dfns in energy and pitch

    
    energy_min=0#80300
    energy_max=1000000#81001
    #pitch_min=-0.3601
    #pitch_max=0.3601
    
    TRANSP_dfn=process_output.dfn_crop(TRANSP_dfn,E=[energy_min,energy_max])
    ASCOT_dfn=ASCOT_dfn.crop(E=[energy_min,energy_max])
    LOCUST_dfn=LOCUST_dfn.crop(E=[energy_min,energy_max])
    #TRANSP_dfn=process_output.dfn_crop(TRANSP_dfn,V_pitch=[pitch_min,pitch_max],inside=False)
    #ASCOT_dfn=ASCOT_dfn.crop(V_pitch=[pitch_min,pitch_max],inside=False)
    #LOCUST_dfn=LOCUST_dfn.crop(V_pitch=[pitch_min,pitch_max],inside=False)    
    print('ASCOT='+str(ASCOT_dfn.transform(axes=['N'])['dfn']))
    print('LOCUST='+str(LOCUST_dfn.transform(axes=['N'])['dfn']))
    print('TRANSP='+str(TRANSP_dfn.dfn_integrate()['dfn']))

    print('ASCOT='+str(ASCOT_dfn['E']))
    print('LOCUST='+str(LOCUST_dfn['E']))
    print('TRANSP='+str(TRANSP_dfn['E']))
    #####

    #####optional crop of TRANSP/LOCUST/ASCOT dfns in R and Z
    
    #'''
    R_min=0#1.6
    R_max=100#1.84
    Z_min=-100#-0.1289
    Z_max=100#0.213

    i=np.where((R_min<TRANSP_dfn['R2D'])&(TRANSP_dfn['R2D']<R_max)&(TRANSP_dfn['Z2D']>Z_min)&(TRANSP_dfn['Z2D']<Z_max))[0]
    for quantity in ['dfn','dVOL','R2D','Z2D']: #remember space is first dimension in TRANSP array, so can do this
        TRANSP_dfn[quantity]=TRANSP_dfn[quantity][i]
    ASCOT_dfn=ASCOT_dfn.crop(R=[R_min,R_max],Z=[Z_min,Z_max])
    LOCUST_dfn=LOCUST_dfn.crop(R=[R_min,R_max],Z=[Z_min,Z_max])

    for ax in [ax1,ax2,ax3]: #show where we've cropped in R
        ax.axvline(R_min,color='m')
        ax.axvline(R_max,color='m')


    #ASCOT_dfn=ASCOT_dfn.crop(R=[1.4,2],Z=[-0.36,0.36],inside=False)
    #LOCUST_dfn=LOCUST_dfn.crop(R=[1.4,2],Z=[-0.36,0.36],inside=False)
    #'''
           

    TRANSP_mesh=TRANSP_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig,real_scale=True)
    ASCOT_mesh=ASCOT_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax2,fig=fig,real_scale=True)
    LOCUST_mesh=LOCUST_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax3,fig=fig,real_scale=True)

if colourbars is True:
    for ax,mesh in zip([ax1,ax2,ax3],[TRANSP_mesh,ASCOT_mesh,LOCUST_mesh]):
        cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
        colourbar_array.append(cbar)

#set labels and plot

for ax in [ax1,ax2,ax3]:
    ax.set_xlabel('R [m]')  
    ax.set_ylabel('Z [m]')  
    ax.set_xlim([np.min(equi['R_1D']),np.max(equi['R_1D'])])
    ax.set_ylim([np.min(equi['Z_1D']),np.max(equi['Z_1D'])])

plt.show()

axes=['R','Z']
fig,((ax1))=plt.subplots(1) #plot the difference
ASCOT_dfn_=ASCOT_dfn.transform(axes=axes)
LOCUST_dfn_=LOCUST_dfn.transform(axes=axes)
DFN_diff=copy.deepcopy(LOCUST_dfn)
DFN_diff.ID='LOCUST dfn - ASCOT dfn'
DFN_diff['dfn']=np.nan_to_num(np.log10(np.abs((LOCUST_dfn_['dfn']-ASCOT_dfn_['dfn'])/LOCUST_dfn_['dfn'])),nan=-5.)
DFN_diff['dfn'][DFN_diff['dfn']>1.e3]=-5.
DFN_diff_mesh=DFN_diff.plot(fig=fig,ax=ax1,axes=axes,transform=False,limiters=wall,LCFS=equi,real_scale=True)
cbar=fig.colorbar(DFN_diff_mesh,orientation='vertical')
ax1.set_xlabel('R [m]',fontsize=25)  
ax1.set_ylabel('Z [m]',fontsize=25)  
ax1.set_title('$log_{10}(f_{LOCUST}-f_{ASCOT})\slash f_{LOCUST}$',fontsize=25)
#ax1.set_xlim([np.min(equi['R_1D']),np.max(equi['R_1D'])])
#ax1.set_ylim([1.1*np.min(equi['lcfs_z']),1.1*np.max(equi['lcfs_z'])])
ax1.set_facecolor(settings.cmap_default(0.0))
plt.show()  