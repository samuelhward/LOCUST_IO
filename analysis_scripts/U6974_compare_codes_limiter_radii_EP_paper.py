#look at RZ, EP and EP sampled plots for TRANSP vs LOCUST vs ASCOT then animate through the limiter radii

import sys
import numpy as np
import context
import copy
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.wall import Wall
from classes.output_classes.distribution_function import Distribution_Function
from classes.output_classes.particle_list import Final_Particle_List
import processing.utils 
import run_scripts.utils
import constants
import settings


#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])
cmap_default=matplotlib.cm.get_cmap('inferno_r')

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

wall_files='input.wall_2d_'
radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
run_IDs=['U69','U70','U71','U72','U73','U74']
shot_number='157418'
colours=['r-','g-','b-','m-','k-','c-']
ascot_coulog=True

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
if ascot_coulog:
    LOCUST_files[0]='F_28-04-2020_23-59-02.590_TOTL.dfn'
    LOCUST_moments[0]='LOCUST_28-04-2020_23-59-02.590.h5'
    LOCUST_run='locust/U69_ascot_coulog/noSOLCOL/'

LOCUST_files=['F_04-12-2018_16-10-58.977_TOTL.dfn','F_04-12-2018_15-06-13.311_TOTL.dfn','F_04-12-2018_15-10-53.134_TOTL.dfn','F_04-12-2018_17-17-41.999_TOTL.dfn','F_04-12-2018_17-19-35.424_TOTL.dfn','F_05-12-2018_01-15-34.057_TOTL.dfn']
LOCUST_moments=['LOCUST_04-12-2018_16-10-58.977.h5','LOCUST_04-12-2018_15-06-13.311.h5','LOCUST_04-12-2018_15-10-53.134.h5','LOCUST_04-12-2018_17-17-41.999.h5','LOCUST_04-12-2018_17-19-35.424.h5','LOCUST_05-12-2018_01-15-34.057.h5']
LOCUST_run='locust/run_2/' #added -DSOLCOL
if ascot_coulog:
    LOCUST_files[0]='F_28-04-2020_18-10-55.938_TOTL.dfn'
    LOCUST_moments[0]='LOCUST_28-04-2020_18-10-55.938.h5'
    LOCUST_run='locust/U69_ascot_coulog/SOLCOL/'


#plot just one radius (still need to comment out desired files where necessary)
rad=0
radii=[radii[rad]]
ASCOT_files=[ASCOT_files[rad]]
run_IDs=[run_IDs[rad]]
LOCUST_files=[LOCUST_files[rad]]

#start by creating some axis objects
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


    fig,ax=plt.subplots(ncols=2,sharex=True)
    axes=['E','V_pitch']

    number_bins=5
    vminmax=[1.e7,6.e7]
    lines=[]
    line_labels=[]
    for dfn,label,cmap in zip([TRANSP_dfn,ASCOT_dfn,LOCUST_dfn],['TRANSP','ASCOT','LOCUST'],[cmap_r,cmap_b,cmap_g]):
        dfn['E']/=1000. #get axes in keV - .plot integrates the DFNs using .transform with 'dE' so should not affect plotting
        mesh=dfn.plot(axes=axes,ax=ax[0],fig=fig,vminmax=vminmax,real_scale=False,fill=False,number_bins=number_bins,colmap=cmap,label=label)
        contours,_ = mesh.legend_elements()    
        lines.append(contours[0])
        line_labels.append(label)
    ax[0].legend(lines,line_labels)
    #ax.set_title('limiter radius = {}'.format(radii[0]))
    ax[0].set_title('Fast ion density $f$',fontsize=25)

    ASCOT_dfn_=ASCOT_dfn.transform(axes=axes)
    LOCUST_dfn_=LOCUST_dfn.transform(axes=axes)
    #interpolate ASCOT grid onto LOCUST grid
    interpolator=processing.utils.interpolate_2D(ASCOT_dfn_['E'],ASCOT_dfn_['V_pitch'],ASCOT_dfn_['dfn'],type='RBF',function='linear')
    V_pitch,E=np.meshgrid(LOCUST_dfn_['V_pitch'],LOCUST_dfn_['E'])
    ASCOT_dfn_['dfn']=interpolator(E,V_pitch)
    ASCOT_dfn_['E'],ASCOT_dfn_['V_pitch']=LOCUST_dfn_['E'],LOCUST_dfn_['V_pitch']
    DFN_diff=copy.deepcopy(LOCUST_dfn)
    DFN_diff.ID='LOCUST dfn - ASCOT dfn'
    DFN_diff['dfn']=np.nan_to_num(np.log10(np.abs((LOCUST_dfn_['dfn']-ASCOT_dfn_['dfn'])/LOCUST_dfn_['dfn'])),nan=-5.)
    DFN_diff['dfn'][DFN_diff['dfn']>1.e3]=-5.
    DFN_diff_mesh=DFN_diff.plot(fig=fig,ax=ax[1],axes=axes,transform=False,vminmax=[-5,2.5])
    ax[1].set_title('$log_{10}(f_{LOCUST}-f_{ASCOT})\slash f_{LOCUST}$',fontsize=25)
    ax[1].set_facecolor(settings.cmap_default(0.0))

    for ax_ in ax:
        ax_.set_xlabel('E [KeV]',fontsize=25)  
        ax_.set_ylabel('$\lambda$',fontsize=25)  
        ax_.set_xlim([np.min(DFN_diff['E']),np.max(DFN_diff['E'])])
        ax_.set_ylim([np.min(DFN_diff['V_pitch']),np.max(DFN_diff['V_pitch'])])
    if colourbars is True:
        for ax,mesh in zip([ax[1]],[DFN_diff_mesh]):
            cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
            colourbar_array.append(cbar)

    plt.show()