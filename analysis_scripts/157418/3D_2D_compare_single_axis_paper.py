import context
import glob
import copy
import support
import settings
import processing.utils
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib import cm
import pathlib
import constants
from scipy.io import readsav
from classes.output_classes.distribution_function import Distribution_Function as DFN
from classes.output_classes.particle_list import Final_Particle_List as FPL
from classes.output_classes.moments import Moments as Mom
from classes.input_classes.equilibrium import Equilibrium as EQ
from classes.input_classes.wall import Wall as Wall
import processing.utils


axis='E'
shot_number='157418'

#read Mike data
filepath_sav=support.dir_output_files / shot_number / 'SPIRAL_output' / 'all_variables_for_sam.sav'
sav=readsav(filepath_sav)

#integrate Mike data
if axis=='V_pitch':
    axis_to_integrate=1
    quantity_to_integrate='ehis'
    quantity_to_plot_against='phis'
    convert_to_eV=1000.
    xlabel='Pitch'
elif axis=='E':
    axis_to_integrate=0
    quantity_to_integrate='phis'
    quantity_to_plot_against='ehis'
    convert_to_eV=1.
    xlabel='Energy [keV]'

sav['surfhist_wc'][:,sav['ehis']>80.5]=0.
sav['surfhist_noc'][:,sav['ehis']>80.5]=0.
surfhist_wc_summed=np.sum(sav['surfhist_wc'],axis=axis_to_integrate)
surfhist_noc_summed=np.sum(sav['surfhist_noc'],axis=axis_to_integrate)
Z=np.nan_to_num((surfhist_wc_summed-surfhist_noc_summed)/surfhist_noc_summed) #normalise
#plot Mike data
fig,ax=plt.subplots(1)
cmap_orange=settings.colour_custom([230,74,25,1])
cmap_g=settings.colour_custom([76,175,80,1])
ax.plot(sav[quantity_to_plot_against],Z,color=cmap_orange(0.))

#remove silly numbers from array
def remove_outliers(array):
    array[array==np.inf]=0.#1.e2
    array[array==-np.inf]=0.#1.e2
    array[array>1.e10]=0.#1.e2
    array[array<-1.e10]=0.#1.e2
    return array

def plot_LOCUST_vs_SPIRAL(
    ax,
    fig,
    response_types=['response_hi_res_90_i3dr-1_n-3_hiPitchres_noDCTDC2_new_profiles_2D_wall','response_hi_res_90_i3dr-1_n-3_hiPitchres_noDCTDC2_new_profiles'],
    linestyles=['--','-'],
    axes='E',
):
    for linestyle,response_type in zip(linestyles,response_types):


        folder_2D='2D_full_slow'  
        folder_3D='3D_full_slow'

        DFN_2D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST'/ response_type / folder_2D).glob('*.dfn'))[0]
        DFN_3D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST'/ response_type / folder_3D).glob('*.dfn'))[0]
        eq_filename=support.dir_input_files / shot_number / 'LOCUST' / 'g157418.03000'

        DFN_2D=DFN(ID='2D',data_format='LOCUST',filename=DFN_2D_filename)
        DFN_3D=DFN(ID='3D',data_format='LOCUST',filename=DFN_3D_filename)
        DFN_2D['E']/=1000. #convert to keV
        DFN_3D['E']/=1000.
        eq=EQ(ID='157418 300ms equilibrium',data_format='GEQDSK',filename=eq_filename,GEQDSKFIX=1)

        rho_crop=0.77 #best fit
        eq.set(flux_pol_norm=(eq['flux_pol']-eq['flux_pol'][0])/(eq['flux_pol'][-1]-eq['flux_pol'][0]))
        eq.set(flux_tor_norm=(eq['flux_tor']-eq['flux_tor'][0])/(eq['flux_tor'][-1]-eq['flux_tor'][0]))
        eq.set(flux_tor_norm_sqrt=np.sqrt(np.abs(eq['flux_tor_norm'])))
        rho_interpolator=processing.utils.interpolate_1D(eq['flux_pol_norm'],eq['flux_tor_norm_sqrt'])
        eq.set(flux_tor_norm_sqrt_rz=processing.utils.flux_func_to_RZ(eq['flux_pol'],eq['flux_tor_norm_sqrt'],eq)) #calculate rho grid
        flux_tor_norm_sqrt_rz=processing.utils.value_at_RZ(*[dimension.T.flatten() for dimension in np.meshgrid(DFN_2D['R'],DFN_2D['Z'])],eq['flux_tor_norm_sqrt_rz'],eq).reshape(len(DFN_2D['R']),len(DFN_2D['Z']))
        DFN_2D['dfn'][...,flux_tor_norm_sqrt_rz<rho_crop]=0.
        DFN_3D['dfn'][...,flux_tor_norm_sqrt_rz<rho_crop]=0.
        DFN_2D['dfn'][:,DFN_2D['E']<10.,:,:,:]=0. #cut off DFN below 10keV
        DFN_3D['dfn'][:,DFN_3D['E']<10.,:,:,:]=0.
        DFN_2D['dfn'][:,DFN_2D['E']>80.5,:,:,:]=0. #cut off DFN above max SPIRAL energy 80.5keV (and seeming cut-off in data at ~82.5keV)
        DFN_3D['dfn'][:,DFN_3D['E']>80.5,:,:,:]=0.

        DFN_diff=copy.deepcopy(DFN_3D)
        DFN_diff['dfn']=DFN_3D.transform(axes=[axes])['dfn']-DFN_2D.transform(axes=[axes])['dfn']
        DFN_diff['dfn']=np.nan_to_num((DFN_diff['dfn'])/DFN_2D.transform(axes=[axes])['dfn']) #normalise
        DFN_diff['dfn']=remove_outliers(DFN_diff['dfn'])
        #DFN_diff['dfn']=(DFN_diff['dfn']-np.max(DFN_diff['dfn']))/(np.min(DFN_diff['dfn'])-np.max(DFN_diff['dfn'])) #normalise
        DFN_diff.plot(ax=ax,fig=fig,axes=[axes],colmap=cmap_g,transform=False,line_style=linestyle)

        if axis=='V_pitch':
            ax.set_xticks(np.linspace(-1,1,5))
            ax.set_xlim([-1,1])
        elif axis=='E':
            ax.set_xticks(np.linspace(10,80,8))
            ax.set_xlim([10,80])

plot_LOCUST_vs_SPIRAL(fig=fig,ax=ax,axes=axis)

legend=[]
legend.append('SPIRAL')
legend.append('LOCUST 2D wall')
legend.append('LOCUST 3D wall')
ax.legend(legend)
#ax.set_yticks([-0.6,-0.4,-0.2,0.])
ax.set_xlabel(xlabel)
ax.set_ylabel('$\delta f$')
ax.set_title('')
plt.show()