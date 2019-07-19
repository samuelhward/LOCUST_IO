#something to compare TRANSP and LOCUST moments

import context
from classes.output_classes.moments import Moments as mom
from classes.output_classes.particle_list import Final_Particle_List as fpl
import constants
import run_scripts.utils
import scipy
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
from mpl_toolkits import mplot3d #import 3D plotting axes
from mpl_toolkits.mplot3d import Axes3D

TRANSP_moments=['157418U69.CDF','157418U70.CDF','157418U71.CDF','157418U72.CDF','157418U73.CDF','157418U74.CDF']


LOCUST_beam_depo_tail='_ptcles.dat'
LOCUST_files=['F_04-12-2018_16-11-28.285_TOTL.dfn','F_03-12-2018_22-55-59.961_TOTL.dfn','F_03-12-2018_22-58-49.281_TOTL.dfn','F_03-12-2018_22-57-37.348_TOTL.dfn','F_04-12-2018_00-13-14.371_TOTL.dfn','F_04-12-2018_14-11-36.625_TOTL.dfn']
LOCUST_moments=['LOCUST_04-12-2018_16-11-28.285.h5','LOCUST_03-12-2018_22-55-59.961.h5','LOCUST_03-12-2018_22-58-49.281.h5','LOCUST_03-12-2018_22-57-37.348.h5','LOCUST_04-12-2018_00-13-14.371.h5','LOCUST_04-12-2018_14-11-36.625.h5']
LOCUST_run='locust/run_1/'

LOCUST_files=['F_04-12-2018_16-10-58.977_TOTL.dfn','F_04-12-2018_15-06-13.311_TOTL.dfn','F_04-12-2018_15-10-53.134_TOTL.dfn','F_04-12-2018_17-17-41.999_TOTL.dfn','F_04-12-2018_17-19-35.424_TOTL.dfn','F_05-12-2018_01-15-34.057_TOTL.dfn']
LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
LOCUST_run='locust/run_2/' #added -DSOLCOL

LOCUST_FINT='FINT.dat_'
LOCUST_run='locust/run_3/' #FINT runs with -DSOLCOL
#LOCUST_run='locust/run_4/' #FINT runs without -DSOLCOL

from U6974_files_ASCOT import * #import all the ASCOT filenames

radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles

colours=['r-','g-','b-','m-','k-','c-']

fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)

beam_source_loss_transp=[] #store the loss fraction at a given time for each limiter radius
beam_source_loss_locust=[]
beam_source_loss_ascot=[]

for LOCUST_moment,TRANSP_moment,ASCOT_file,radius,colour in zip(LOCUST_moments,TRANSP_moments,ASCOT_files,radii,colours):

    #read in the data

    mom_transp=mom(ID='radius = '+radius,data_format='TRANSP',filename=TRANSP_moment)
    mom_transp['time']-=3. #TRANSP runs on experiment time via OMFIT instead of simulation time like ASCOT and LOCUST
    mom_locust=run_scripts.utils.FINT_LOCUST(ID='LOCUST moments - radius = '+radius,filename=LOCUST_run+LOCUST_FINT+radius)
    mom_ascot=fpl(ID='ASCOT final particle list - radius = '+radius,data_format='ASCOT',filename=ASCOT_run+ASCOT_file) #get the prompt losses from the ASCOT particle list output

    
    #scale data to PINJ=PABS=1W
    
    beam_power=1. #this is desired power captured (=power injected in LOCUST)
    beam_energy=80000
    for quantity in ['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source','PFC_power']: #need to scale these up by 1/BCAP
        mom_transp[quantity]*=beam_power/mom_transp['beam_source_captured']

    j=np.where(mom_ascot['status_flag']>0)[0] #scale by all markers which actually contributed to the simulation
    captured_power_ascot=np.sum(beam_energy*constants.species_charge*mom_ascot['weight'][j]) #also need to scale ASCOT weights accordingly, since the markers should make up 100% of the injected weight (unlike in LOCUST and TRANSP)
    mom_ascot['weight']*=beam_power/captured_power_ascot

    i=np.where(mom_ascot['status_flag']==3)[0] #separate ASCOT markers which have hit PFC components    
    mom_ascot['PFC_power']=np.histogram(mom_ascot['time'][i],bins=len(mom_locust['time']),weights=mom_ascot['E'][i]*constants.species_charge*mom_ascot['weight'][i])[0] #place against same time base as LOCUST
    mom_ascot['PFC_power']=np.cumsum(mom_ascot['PFC_power'])

    #plot PFC loss power vs time

    ax1.plot(mom_locust['time'],mom_locust['PFC_power'],colour)
    ax2.plot(mom_transp['time'],beam_power*mom_transp['PFC_power'],colour)
    ax3.plot(mom_locust['time'],mom_ascot['PFC_power'],colour) #ASCOT and LOCUST on same time base

    #grab the PFC loss power power at a specific time for plotting vs radii later 

    time=3. #TRANSP time offset by 3s since TRANSP is real tokamak pulse time whereas LOCUST is simulation time
    time_index_transp=(np.abs(mom_transp['time']-time)).argmin() #find nearest timestep to one of interest
    time=0.1 #time at which LOCUST and ASCOT moments are written out
    time_index_locust=(np.abs(mom_locust['time']-time)).argmin()
    time_index_ascot=time_index_locust #ASCOT and LOCUST using the same timebase for PFC_power array

    beam_source_loss_transp.append(mom_transp['PFC_power'][time_index_transp])
    beam_source_loss_locust.append(mom_locust['PFC_power'][time_index_locust])
    beam_source_loss_ascot.append(mom_ascot['PFC_power'][time_index_locust])

#plot PFC loss power vs radii for each code

ax4.plot([1.05,1.1,1.2,1.3,1.4,1.5],beam_source_loss_locust,'g')
ax4.plot([1.05,1.1,1.2,1.3,1.4,1.5],beam_source_loss_transp,'r')
ax4.plot([1.05,1.1,1.2,1.3,1.4,1.5],beam_source_loss_ascot,'b')

ax1.legend(radii)
ax1.set_title('LOCUST')
ax1.set_ylim([0.,0.15])
ax1.set_ylabel('beam source loss [W]')

ax2.set_title('TRANSP')
ax2.set_ylim([0.,0.15])

ax3.set_title('ASCOT')
ax3.set_ylim([0.,0.15])

ax4.set_xlabel('radius')
ax4.set_ylabel('total PFC power [W] at time={}s'.format(time))
ax4.legend(['LOCUST','TRANSP','ASCOT'])

#ax4.set_xlabel('time [s]')
#ax4.set_ylabel('beam source capture fraction')

plt.show()
