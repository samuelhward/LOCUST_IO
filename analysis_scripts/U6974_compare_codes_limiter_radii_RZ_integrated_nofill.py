#look at RZ, EP and EP sampled plots for TRANSP vs LOCUST vs ASCOT then animate through the limiter radii

import sys
import numpy as np
import context
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
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
#define some colourmaps
colors=[(1,0, 0), (1, 0, 0)]
reds=LinearSegmentedColormap.from_list('reds', colors, N=10)
colors=[(0,1, 0), (0, 1, 0)]
blues=LinearSegmentedColormap.from_list('blues', colors, N=10)
colors=[(0,0, 1), (0, 0, 1)]
greens=LinearSegmentedColormap.from_list('greens', colors, N=10)


filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)



wall_files='input.wall_2d_'
radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
run_IDs=['U69','U70','U71','U72','U73','U74']
shot_number='157418'

colours=['r-','g-','b-','m-','k-','c-']

TRANSP_files_tail_FI='_fi_1_gc.cdf'
TRANSP_files_tail_CDF='.CDF'

from U6974_files_ASCOT import * #import all the ASCOT filenames

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
fig,(ax1)=plt.subplots(1,1)

for radius,LOCUST_file,ASCOT_file,run_ID,colour in zip(radii,LOCUST_files,ASCOT_files,run_IDs,colours):

    #for ax in [ax1,ax2,ax3]:
    #    ax.cla()

    #toggle colourbars
    colourbars=False
    colourbar_array=[]

    LOCUST_dfn=Distribution_Function('LOCUST density (0.1s) - r = {}'.format(str(radius)),'LOCUST',filename=LOCUST_run+LOCUST_file)
    TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP density (0.1s) - r = {}'.format(str(radius)),filename=shot_number+run_ID+TRANSP_files_tail_FI)
    ASCOT_dfn=run_scripts.utils.ASCOT_output('ASCOT density (0.1s) - r = {}'.format(str(radius)),filename=ASCOT_run+ASCOT_file,datatype='distribution_function')

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

    vminmax=[0.1e11,7.e11]

    #####optional crop of dfns in energy

    
    energy_min=18000
    energy_max=100000

    #TRANSP_dfn=process_output.dfn_crop(TRANSP_dfn,E=[energy_min,energy_max])
    #ASCOT_dfn=process_output.dfn_crop(ASCOT_dfn,E=[energy_min,energy_max])
    #LOCUST_dfn=process_output.dfn_crop(LOCUST_dfn,E=[energy_min,energy_max])

    #for ax in [ax4,ax5,ax6,ax7,ax8,ax9]: #show where we've cropped in energy
        #ax.axvline(energy_min,color='m')
        #ax.axvline(energy_max,color='m')

    #####

    #####optional crop of LOCUST/ASCOT dfns in R and Z
    

    R_min=0.
    R_max=2.
    Z_min=-0.51
    Z_max=0.51

    i=np.where((R_min<TRANSP_dfn['R2D'])&(TRANSP_dfn['R2D']<R_max)&(TRANSP_dfn['Z2D']>Z_min)&(TRANSP_dfn['Z2D']<Z_max))[0]
    #for quantity in ['dfn','dVOL','R2D','Z2D']: #remember space is first dimension in TRANSP array, so can do this
        #TRANSP_dfn[quantity]=TRANSP_dfn[quantity][i]
    #ASCOT_dfn=process_output.dfn_crop(ASCOT_dfn,R=[R_min,R_max],Z=[Z_min,Z_max])
    #LOCUST_dfn=process_output.dfn_crop(LOCUST_dfn,R=[R_min,R_max],Z=[Z_min,Z_max])

    #for ax in [ax1,ax2,ax3]: #show where we've cropped in R
        #ax.axvline(R_min,color='m')
        #ax.axvline(R_max,color='m')

    #####       

    number_bins=6

    TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig,vminmax=vminmax,real_scale=True,fill=False,number_bins=number_bins,colmap=reds)
    ASCOT_mesh=ASCOT_dfn.dfn_plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig,vminmax=vminmax,real_scale=True,fill=False,number_bins=number_bins,colmap=blues)
    LOCUST_mesh=LOCUST_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig,vminmax=vminmax,real_scale=True,fill=False,number_bins=number_bins,colmap=greens)

    #add scatter point if you fancy

    r_sample_points=[1.7835,2.18]
    z_sample_points=[0.0143,0.37]

    for r_sample_point,z_sample_point in zip(r_sample_points,z_sample_points):


        #for point we have chosen, find closest dfn bin in real space
        diff_r=np.abs(TRANSP_dfn['R2D']-r_sample_point)**2 #irregular grid so need to minimise distance to every point
        diff_z=np.abs(TRANSP_dfn['Z2D']-z_sample_point)**2
        index_rz_transp=(diff_r+diff_z).argmin()
        r_transp=TRANSP_dfn['R2D'][index_rz_transp]
        z_transp=TRANSP_dfn['Z2D'][index_rz_transp]

        index_r_ascot=np.abs(ASCOT_dfn['R']-r_sample_point).argmin() #regular grids so can use this method
        index_z_ascot=np.abs(ASCOT_dfn['Z']-z_sample_point).argmin() 
        r_ascot=ASCOT_dfn['R'][index_r_ascot]
        z_ascot=ASCOT_dfn['Z'][index_z_ascot]

        index_r_locust=np.abs(LOCUST_dfn['R']-r_sample_point).argmin()
        index_z_locust=np.abs(LOCUST_dfn['Z']-z_sample_point).argmin()
        r_locust=LOCUST_dfn['R'][index_r_locust]
        z_locust=LOCUST_dfn['Z'][index_z_locust]

        #ax1.scatter(r_transp,z_transp,color='w',s=20)
        #ax1.scatter(r_ascot,z_ascot,color='w',s=20)
        #ax1.scatter(r_locust,z_locust,color='w',s=20)


    if colourbars is True:
        for ax,mesh in zip([ax1],[LOCUST_mesh]):
            cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
            colourbar_array.append(cbar)

    #set labels and plot

    for ax in [ax1]:
        ax.set_title('limiter radius = {}'.format(radii[0]))
        ax.set_xlabel('R [m]')  
        ax.set_ylabel('Z [m]')  
        ax.set_xlim([np.min(equi['R_1D']),np.max(equi['R_1D'])])
        ax.set_ylim([1.1*np.min(equi['lcfs_z']),1.1*np.max(equi['lcfs_z'])])

    plt.draw()
    plt.pause(0.0001)
    #fig.clear()

plt.show()
