#animates through limiter radii for a single output time

#go to comparisons and see how TRANSP stores it's time on its moments - and whether it's from 0 or 3 for example

import context
import constants
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
time=3.1 #time at which LOCUST and ASCOT moments are written out

TRANSP_moment_array=[]
LOCUST_moment_array=[]
ASCOT_moment_array=[]

#quantities=['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source','energy','energy_para','energy_perp','beam_source'] #shared by TRANSP and LOCUST
quantities=['NBI-heating-power(i1)','NBI-heating-power(e-)','energy','energy_para'] #shared by ASCOT, TRANSP and LOCUST

fig,(ax1,ax2,ax3,ax4)=plt.subplots(1,len(quantities)) 
axes=[ax1,ax2,ax3,ax4]

for LOCUST_moment,TRANSP_moment,ASCOT_moment,radius in zip(LOCUST_moments,TRANSP_moments,ASCOT_moments,radii):

    mom_locust=mom(ID='LOCUST radius = '+radius,data_format='LOCUST',filename=LOCUST_run+LOCUST_moment)
    locust_beam_power=1. #need to make sure that beam powers are matched and deposition fraction matches too in both codes - LOCUST been assuming 1 so far
    
    mom_transp=mom(ID='TRANSP radius = '+radius,data_format='TRANSP',filename=TRANSP_moment)
    time_index_transp=(np.abs(mom_transp['time']-time)).argmin() #find nearest timestep to one of interest
    
    mom_ascot=mom(ID='ASCOT radius = '+radius,data_format='ASCOT',filename=ASCOT_run+ASCOT_moment)
    beam_depo=Final_Particle_List(ID='',data_format='ASCOT',filename=ASCOT_run+ASCOT_moment) #likewise we need to do the same for ASCOT too
    beam_power=1. #1 Watt beam power
    beam_energy=80000
    k=np.where(beam_depo['status_flag']>0)[0] #scale by all markers which actually contributed to the simulation
    BPCAP_ascot=np.sum(beam_energy*constants.species_charge*beam_depo['weight'][k])

    for quantity,ax in zip(quantities,axes):

        mom_transp[quantity][time_index_transp]*=locust_beam_power/mom_transp['beam_source_captured'][time_index_transp] #normalise the TRANSP quantities according to 100% captured beam power
        mom_ascot[quantity]*=beam_power/(BPCAP_ascot) #XXX this needs sorting
        ax.plot(mom_transp['flux_pol_norm'][time_index_transp],mom_transp[quantity][time_index_transp],'r')
        ax.plot(mom_locust['flux_pol_norm'],mom_locust[quantity],'g')
        ax.plot(mom_ascot['flux_pol_norm'],mom_ascot[quantity],'b')

        ax.legend(('TRANSP','LOCUST','ASCOT'))
        ax.set_title(quantity)

    ax1.set_ylabel('radius={}'.format(radius))
    plt.draw()
    plt.pause(5)
    plt.title('radius = {radius}, time = {time}'.format(radius=radius,time=mom_transp['time'][time_index_transp]))

    for ax in axes:
        ax.cla()

