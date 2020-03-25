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

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)



wall_files='input.wall_2d_'
radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
#radii=['1.50']
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

#LOCUST_files=['F_04-12-2018_16-10-58.977_TOTL.dfn','F_04-12-2018_15-06-13.311_TOTL.dfn','F_04-12-2018_15-10-53.134_TOTL.dfn','F_04-12-2018_17-17-41.999_TOTL.dfn','F_04-12-2018_17-19-35.424_TOTL.dfn','F_05-12-2018_01-15-34.057_TOTL.dfn']
#LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
#LOCUST_run='locust/run_2/' #added -DSOLCOL

#start by creating some axis objects
fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(ax10,ax11,ax12))=plt.subplots(4,3)


#sample point for second EP plot
r_sample_points=[1.797,2.0,2.1,1.3,1.5]
z_sample_points=[0.396,0.396,0.396,0.,0.]

#r_sample_points=[2.0511]
#z_sample_points=[0.1]

r_sample_points=[1.7835]
z_sample_points=[0.0143]

for r_sample_point,z_sample_point in zip(r_sample_points,z_sample_points):

    for ax in [ax10,ax11,ax12]:
        ax.cla()

    for radius,LOCUST_file,ASCOT_file,run_ID,colour in zip(radii,LOCUST_files,ASCOT_files,run_IDs,colours):

        for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]:
            ax.cla()

        #toggle colourbars
        colourbars=False
        colourbar_array=[]

        LOCUST_dfn=Distribution_Function('LOCUST fast ion density (0.1s) - limiter radius = {}'.format(str(radius)),'LOCUST',filename=LOCUST_run+LOCUST_file)
        TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density (0.1s) - limiter radius = {}'.format(str(radius)),filename=shot_number+run_ID+TRANSP_files_tail_FI)
        ASCOT_dfn=run_scripts.utils.ASCOT_output('ASCOT fast ion density (0.1s) - limiter radius = {}'.format(str(radius)),filename=ASCOT_run+ASCOT_file,datatype='distribution_function')

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

        vminmax=[0.,8.e11]

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


        TRANSP_mesh=TRANSP_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig,vminmax=vminmax)
        ASCOT_mesh=ASCOT_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax2,fig=fig,vminmax=vminmax)
        LOCUST_mesh=LOCUST_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax3,fig=fig,vminmax=vminmax)

        if colourbars is True:
            for ax,mesh in zip([ax1,ax2,ax3],[TRANSP_mesh,ASCOT_mesh,LOCUST_mesh]):
                cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
                colourbar_array.append(cbar)

        #sampled EP

        vminmax=[0.,4.5e7]

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

        ax1.scatter(r_transp,z_transp,color='m',s=5)
        ax1.scatter(r_ascot,z_ascot,color='g',s=5)
        ax1.scatter(r_locust,z_locust,color='y',s=5)

        ax2.scatter(r_transp,z_transp,color='m',s=5)
        ax2.scatter(r_ascot,z_ascot,color='g',s=5)
        ax2.scatter(r_locust,z_locust,color='y',s=5)

        ax3.scatter(r_transp,z_transp,color='m',s=5)
        ax3.scatter(r_ascot,z_ascot,color='g',s=5)
        ax3.scatter(r_locust,z_locust,color='y',s=5)

        axes=[index_rz_transp,slice(None),slice(None)]
        TRANSP_mesh=TRANSP_dfn.plot(axes=axes,ax=ax7,fig=fig,vminmax=vminmax) #plotting with slices and specifying the point plots from the non-integrated distribution function

        axes=[slice(None),slice(None),index_r_ascot,index_z_ascot]
        ASCOT_mesh=ASCOT_dfn.plot(axes=axes,ax=ax8,fig=fig,vminmax=vminmax)

        axes=[0,slice(None),slice(None),index_r_locust,index_z_locust]
        LOCUST_mesh=LOCUST_dfn.plot(axes=axes,ax=ax9,fig=fig,vminmax=vminmax)

        if colourbars is True:
            for ax,mesh in zip([ax7,ax8,ax9],[TRANSP_mesh,ASCOT_mesh,LOCUST_mesh]):
                cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
                colourbar_array.append(cbar)

        #EP

        axes=['E','V_pitch']

        vminmax=[0.,0.8e8]

        TRANSP_mesh=TRANSP_dfn.plot(axes=axes,ax=ax4,fig=fig,vminmax=vminmax)
        ASCOT_mesh=ASCOT_dfn.plot(axes=axes,ax=ax5,fig=fig,vminmax=vminmax)
        LOCUST_mesh=LOCUST_dfn.plot(axes=axes,ax=ax6,fig=fig,vminmax=vminmax)

        if colourbars is True:
            for ax,mesh in zip([ax4,ax5,ax6],[TRANSP_mesh,ASCOT_mesh,LOCUST_mesh]):
                cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
                colourbar_array.append(cbar)







        #E

        axes=['E']

        #either plotting at a single point
        
        #'''
        TRANSP_plot=TRANSP_dfn.dfn_integrate(space=False,energy=False) #crop and plot the TRANSP dfn
        dVOL=2.*constants.pi*LOCUST_dfn['R'][index_r_locust]*LOCUST_dfn['dR']*LOCUST_dfn['dZ'] #TRANSP volume elements dVOL are different sizes to ASCOT/LOCUST so need to do some scaling
        TRANSP_plot['dfn']=TRANSP_plot['dfn'][index_rz_transp,:]*dVOL #integrate over cell of interest and scale according to LOCUST volume elements
        ASCOT_dfn=process_output.dfn_crop(ASCOT_dfn,R=[r_sample_point],Z=[z_sample_point]) #crop and plot the ASCOT dfn
        LOCUST_dfn=process_output.dfn_crop(LOCUST_dfn,R=[r_sample_point],Z=[z_sample_point]) #crop and plot the LOCUST dfn
        #'''

        #or plotting over all current space
        '''
        TRANSP_plot=TRANSP_dfn.dfn_integrate(energy=False) #crop and plot the TRANSP dfn    
        '''

        vminmax=[0,np.amax(1.1*TRANSP_plot['dfn'])]

        ax10.plot(TRANSP_plot['E'],TRANSP_plot['dfn'],'r')
        ASCOT_dfn.plot(axes=axes,ax=ax10,fig=fig,colmap='b')
        LOCUST_dfn.plot(axes=axes,ax=ax10,fig=fig,colmap='g')

        ax10.set_title('')
        ax10.set_xlabel('E [eV]')
        ax10.set_ylabel('[#/eV]')
        ax10.set_ylim(vminmax)
        #ax10.set_ylim([0,np.max(TRANSP_plot['dfn'])])
        #ax10.set_ylim([0,np.max([np.max(dfn) for dfn in [TRANSP_plot['dfn'],ASCOT_dfn['dfn'],LOCUST_dfn['dfn']]])]) #find current max across arrays






        #pitch
        axes=['V_pitch']

        #either plotting at a single point
        
        #'''
        TRANSP_plot=TRANSP_dfn.dfn_integrate(space=False,pitch=False) #crop and plot the TRANSP dfn
        dVOL=2.*constants.pi*LOCUST_dfn['R'][index_r_locust]*LOCUST_dfn['dR']*LOCUST_dfn['dZ'] #TRANSP volume elements dVOL are different sizes to ASCOT/LOCUST so need to do some scaling
        TRANSP_plot['dfn']=TRANSP_plot['dfn'][index_rz_transp,:]*dVOL #integrate over cell of interest and scale according to LOCUST volume elements
        ASCOT_dfn=process_output.dfn_crop(ASCOT_dfn,R=[r_sample_point],Z=[z_sample_point]) #crop and plot the ASCOT dfn
        LOCUST_dfn=process_output.dfn_crop(LOCUST_dfn,R=[r_sample_point],Z=[z_sample_point]) #crop and plot the LOCUST dfn
        #'''

        #or plotting over all current space
        '''
        TRANSP_plot=TRANSP_dfn.dfn_integrate(pitch=False) #crop and plot the TRANSP dfn    
        '''

        vminmax=[0,np.amax(1.1*TRANSP_plot['dfn'])]

        ax11.plot(TRANSP_plot['V_pitch'],TRANSP_plot['dfn'],'r')
        ASCOT_dfn.plot(axes=axes,ax=ax11,fig=fig,colmap='b')
        LOCUST_dfn.plot(axes=axes,ax=ax11,fig=fig,colmap='g')

        ax11.set_title('')
        ax11.set_xlabel('V_pitch [V$_{\|\|}$/V]')
        ax11.set_ylabel('[#/dPitch]')
        ax11.set_ylim(vminmax)
        #ax11.set_ylim([0,np.max(TRANSP_plot['dfn'])])
        #ax11.set_ylim([0,np.max([np.max(dfn) for dfn in [TRANSP_plot['dfn'],ASCOT_dfn['dfn'],LOCUST_dfn['dfn']]])]) #find current max across arrays






             
        #set labels and plot

        for ax in [ax1,ax2,ax3]:
            ax.set_xlabel('R [m]')  
            ax.set_ylabel('Z [m]')  
            ax.set_xlim([np.min(equi['R_1D']),np.max(equi['R_1D'])])
            ax.set_ylim([np.min(equi['Z_1D']),np.max(equi['Z_1D'])])

        for ax in [ax4,ax5,ax6,ax7,ax8,ax9]:
            ax.set_title('')
            ax.set_xlabel('E [eV]')
            ax.set_ylabel('V$_{\|\|}$/V')
            ax.set_xlim([0,100000])
            ax.set_ylim([-0.8,0.8])

        plt.draw()
        plt.pause(0.0001)
        #fig.clear()


    ax10.legend(tuple(['TRANSP','ASCOT','LOCUST']))
    ax11.legend(tuple(['TRANSP','ASCOT','LOCUST']))
    
    plt.draw()
