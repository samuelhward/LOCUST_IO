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

ascot_coulog=False

TRANSP_moments=['157418U69.CDF','157418U70.CDF','157418U71.CDF','157418U72.CDF','157418U73.CDF','157418U74.CDF']

LOCUST_moments=['LOCUST_04-12-2018_16-11-28.285.h5','LOCUST_03-12-2018_22-55-59.961.h5','LOCUST_03-12-2018_22-58-49.281.h5','LOCUST_03-12-2018_22-57-37.348.h5','LOCUST_04-12-2018_00-13-14.371.h5','LOCUST_04-12-2018_14-11-36.625.h5']
LOCUST_run='locust/run_1/'
if ascot_coulog:
    LOCUST_moments[0]='LOCUST_28-04-2020_23-59-02.590.h5'
    LOCUST_run='locust/U69_ascot_coulog/noSOLCOL/'

LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
LOCUST_run='locust/run_2/' #added -DSOLCOL
if ascot_coulog:
    LOCUST_moments[0]='LOCUST_28-04-2020_18-10-55.938.h5'
    LOCUST_run='locust/U69_ascot_coulog/SOLCOL/'

ASCOT_files=['ascot_freia_1470025.h5','ascot_freia_1470029.h5','ascot_freia_1470032.h5','ascot_freia_1470036.h5','ascot_freia_1470040.h5','ascot_freia_1470044.h5']
ASCOT_run='ascot/run_1/' #this is with old kinetic profiles which are not extrapolated, ORBITMETHOD=1

ASCOT_files=['ascot_freia_1470026.h5','ascot_freia_1470030.h5','ascot_freia_1470033.h5','ascot_freia_1470037.h5','ascot_freia_1470041.h5','ascot_freia_1470045.h5']
ASCOT_run='ascot/run_2/' #changed ORBITMETHOD to 4, added extrapolated kinetic profiles

ASCOT_files=['ascot_freia_1470027.h5','ascot_freia_1470031.h5','ascot_freia_1470034.h5','ascot_freia_1470038.h5','ascot_freia_1470042.h5','ascot_freia_1470046.h5']
ASCOT_run='ascot/run_3/' #changed ORBITMETHOD back to 1, keep extrapolated kinetic profiles

ASCOT_files=['ascot_freia_1470028.h5','ascot_freia_1480719.h5','ascot_freia_1470035.h5','ascot_freia_1470039.h5','ascot_freia_1470043.h5','ascot_freia_1470047.h5']
ASCOT_run='ascot/run_4/' #ORBITMETHOD 4 and non-extrapolated kinetic profiles

ASCOT_moments=ASCOT_files #rename this variable just because it's called something different in the ASCOT files file


radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
time=3.1 #time at which LOCUST and ASCOT moments are written out

rad=0 #XXX
radii=[radii[rad]]
ASCOT_moments=[ASCOT_moments[rad]]
LOCUST_moments=[LOCUST_moments[rad]]
TRANSP_moments=[TRANSP_moments[rad]]

TRANSP_moment_array=[]
LOCUST_moment_array=[]
ASCOT_moment_array=[]

quantities=['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source','torque-density(JxB-sweep)','J(NBCD)-raw'] #shared by TRANSP and LOCUST
#quantities=['NBI-heating-power(i1)','NBI-heating-power(e-)','energy','torque-density(JxB-sweep)'] #shared by ASCOT, TRANSP and LOCUST
#quantities=['NBI-heating-power(i1)','NBI-heating-power(e-)','energy','torque-density(JxB-sweep)'] #custom

fig,axes=plt.subplots(2,int(len(quantities)/2)) 
axes=list(axes)
axes=[_ for axis in axes for _ in axis]

for LOCUST_moment,TRANSP_moment,ASCOT_moment,radius in zip(LOCUST_moments,TRANSP_moments,ASCOT_moments,radii):

    mom_locust=mom(ID='LOCUST radius = '+radius,data_format='LOCUST',filename=LOCUST_run+LOCUST_moment)
    locust_beam_power=1. #need to make sure that beam powers are matched and deposition fraction matches too in both codes - LOCUST been assuming 1 so far
    
    mom_transp=mom(ID='TRANSP radius = '+radius,data_format='TRANSP',filename=TRANSP_moment)
    mom_transp.look()
    print(mom_transp['flux_pol'][-1])
    print(mom_transp['flux_pol'][0])
    time_index_transp=(np.abs(mom_transp['time']-time)).argmin() #find nearest timestep to one of interest
    
    mom_ascot=mom(ID='ASCOT radius = '+radius,data_format='ASCOT',filename=ASCOT_run+ASCOT_moment)

    beam_depo=Final_Particle_List(ID='',data_format='ASCOT',filename=ASCOT_run+ASCOT_moment) #likewise we need to do the same for ASCOT too
    beam_power=1. #1 Watt beam power
    beam_energy=80000
    k=np.where(beam_depo['status_flag']>0)[0] #scale by all markers which actually contributed to the simulation
    BPCAP_ascot=np.sum(beam_energy*constants.species_charge*beam_depo['weight'][k])

    for quantity,ax in zip(quantities,axes):

        legend=[]

        try:
            mom_transp[quantity][time_index_transp]*=locust_beam_power/mom_transp['beam_source_captured'][time_index_transp] #normalise the TRANSP quantities according to 100% captured beam power
            ax.plot(mom_transp['flux_pol_norm'][time_index_transp],mom_transp[quantity][time_index_transp],'r')
            legend.append('TRANSP')
        except:
            pass
        try:
            ax.plot(mom_locust['flux_pol_norm'],mom_locust[quantity],'g')
            legend.append('LOCUST')
        except:
            pass
        try:
            mom_ascot[quantity]*=beam_power/(BPCAP_ascot) #XXX this needs sorting
            #ax.plot(mom_ascot['flux_pol_norm'],mom_ascot[quantity],'b')
            #legend.append('ASCOT')
        except:
            pass

        ax.legend(legend,fontsize=10)
        ax.set_title(quantity)

    axes[0].set_ylabel('radius={}'.format(radius))
    plt.show() #XXX
    #plt.draw()
    #plt.pause(5)
    plt.title('radius = {radius}, time = {time}'.format(radius=radius,time=mom_transp['time'][time_index_transp]))

    for ax in axes:
        ax.cla()

