#something to compare moments through time

import context
from classes.output_classes.moments import Moments as mom 
from classes.output_classes.particle_list import Final_Particle_List
import scipy
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
from mpl_toolkits import mplot3d #import 3D plotting axes
from mpl_toolkits.mplot3d import Axes3D

TRANSP_moments=['157418U69.CDF','157418U70.CDF','157418U71.CDF','157418U72.CDF','157418U73.CDF','157418U74.CDF']

LOCUST_moments=['LOCUST_04-12-2018_16-11-28.285.h5','LOCUST_03-12-2018_22-55-59.961.h5','LOCUST_03-12-2018_22-58-49.281.h5','LOCUST_03-12-2018_22-57-37.348.h5','LOCUST_04-12-2018_00-13-14.371.h5','LOCUST_04-12-2018_14-11-36.625.h5']
LOCUST_run='locust/run_1/'

LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
LOCUST_run='locust/run_2/' #added -DSOLCOL

from U6974_files_ASCOT import * #import all the ASCOT filenames
ASCOT_moments=ASCOT_files #rename this variable just because it's called something different in the ASCOT files file

radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
colours=['r-','g-','b-','m-','k-','c-']
TRANSP_moment_array=[]
LOCUST_moment_array=[]
ASCOT_moment_array=[]

fig,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8),(ax9,ax10,ax11,ax12))=plt.subplots(3,4) 

for LOCUST_moment,TRANSP_moment,ASCOT_moment,radius in zip(LOCUST_moments,TRANSP_moments,ASCOT_moments,radii):

    mom_transp=mom(ID='radius = '+radius,data_format='TRANSP',filename=TRANSP_moment)
    #mom_locust=mom(ID='LOCUST moments - radius = '+radius,data_format='LOCUST',filename=LOCUST_run+LOCUST_moment)
    mom_locust=0. #for now
    mom_ascot=mom(ID='radius = '+radius,data_format='ASCOT',filename=ASCOT_run+ASCOT_moment)

    locust_beam_power=1. #need to make sure that beam powers are match and deposition fraction matches too in both codes - LOCUST been assuming 1 so far
    for key in ['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source']:
        mom_transp[key]*=locust_beam_power/mom_transp['beam_source_captured']

    beam_depo=Final_Particle_List(ID=run_ID,data_format='ASCOT',filename=ASCOT_run+ASCOT_moment) #likewise we need to do the same for ASCOT too
    beam_power=1. #1 Watt beam power
    beam_energy=80000
    k=np.where(beam_depo['status_flag']>0)[0] #scale by all markers which actually contributed to the simulation
    BPCAP_ascot=np.sum(beam_energy*constants.species_charge*beam_depo['weight'][k])
    for key in []: #XXX this needs sorting
        mom_ascot[key]*=beam_power/(BPCAP_ascot)

    TRANSP_moment_array.append(mom_transp)
    LOCUST_moment_array.append(mom_locust)
    ASCOT_moment_array.append(mom_ascot)


for time in TRANSP_moment_array[0]['time']:
    axes=[[ax1,ax2,ax3,ax4],[ax5,ax6,ax7,ax8],[ax9,ax10,ax11,ax12]]
    axes=list(map(list,zip(*axes))) #stride the axes i.e [[a,b],[c,d]] becomes [[a,c],[b,d]]
    time_index_transp=(np.abs(TRANSP_moment_array[0]['time']-time)).argmin()
    #time_index_locust=(np.abs(LOCUST_moment_array[0]['time']-time)).argmin()

    for ax_column,quantity in zip(axes,['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source']):
            for radius,mom_transp,mom_locust,mom_ascot,colour in zip(radii,TRANSP_moment_array,LOCUST_moment_array,ASCOT_moment_array,colours):
                ax_column[0].plot(mom_transp['flux_pol_norm_sqrt'][time_index_transp],mom_transp[quantity][time_index_transp],colour)
                ax_column[0].set_title(quantity+'(t={})'.format(str(time)))
                #ax_column[1].plot(mom_locust['flux_pol_norm_sqrt'],mom_locust[quantity],'g-')
                #ax_column[2].plot(mom_ascot['flux_pol_norm_sqrt'],mom_ascot[quantity],'g-'))
    
    ax1.legend(radii)
    ax1.set_xlabel('TRANSP')
    ax5.set_xlabel('LOCUST')
    ax9.set_xlabel('ASCOT')
    #ax5.legend(radii)
    ax3.set_xlim([0.8,1.0])        
    plt.draw()
    plt.pause(0.1)
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]:
        ax.cla()
