#3D_2D_compare_individual.py
 
"""
Samuel Ward
22/04/2019
----
compare results from 2D and 3D runs
---
usage:
 
todo:
    make into command-line arguement script
    add comments
    add ability to plot all 3 types of field
notes:


i3dr=+1

cp = np.cos(mode_number*phi_point*i3dr-90) = cos(n phi - 90 ) = sin(n phi) = - sin(n phi) for n=-3
sp = np.sin(mode_number*phi_point*i3dr-90) = sin(n phi - 90 ) = -cos(n phi) = - cos(n phi) for n=-3


i3dr=-1

cp = np.cos(-mode_number*phi_point-90) = cos (n phi + 90) = -sin(n phi) = sin (n phi) for n=-3
sp = np.sin(-mode_number*phi_point-90) =-sin( n phi + 90) = -cos(n phi) = - cos (n phi) for n=-3



Re*cp - Im*sp --> - ReCp - ImSp

so at phi = pi/2 we should see the same thing!

if real component is small we should also get the same answer - and they are! the only time the real component is large is for the toroidal component....
which is the same in this case! hence i3dr is a very very small correction

likewise if imaginary componetn ~ 0 then this will give the same thing

-13% for response_hi_res_90_i3dr-1

-11% for response_hi_res_90_i3dr-1_n-3_noDCTDC2

-4% for response_hi_res_90_i3dr+1_n-3_noDCTDC2


our evaluate function is only a bit different to LOCUST's


MVZ has just picked a part of space with symmetry unfortunately, so we can't tell the difference between two switches
---
"""

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

cmap_reds=matplotlib.cm.get_cmap('Reds_r')
cmap_blues=matplotlib.cm.get_cmap('Blues_r')
cmap_IDL=matplotlib.cm.get_cmap('CMRmap')

from classes.output_classes.distribution_function import Distribution_Function as DFN
from classes.output_classes.particle_list import Final_Particle_List as FPL
from classes.output_classes.moments import Moments as Mom
from classes.input_classes.equilibrium import Equilibrium as EQ
from classes.input_classes.wall import Wall as Wall
import processing.utils


shot_number='157418'

#read Mike data
filepath_sav=support.dir_output_files / shot_number / 'SPIRAL_output' / 'all_variables_for_sam.sav'
sav=readsav(filepath_sav)

#response_type='response'    
#response_type='response_lo_res_90_i3dr-1_n-3_noDCTDC2_hiPitchres_notune'
#response_type='vacuum'
#response_type='response_hi_res'
#response_type='response_hi_res_-90_i3dr1'
#response_type='response_hi_res_90_i3dr-1_n-3' #old best fit
#response_type='response_hi_res_90_i3dr-1_n-3_hiPitchres' #increased resolution in pitch and remove BRELAX
#response_type='response_hi_res_90_i3dr+1_n-3' #need to run

#latest runs
#response_type='response_hi_res_90_i3dr-1_n-3_noDCTDC2_2D_wall'
#response_type='response_hi_res_90_i3dr-1_n-3_noDCTDC2'  #~40% max loss
#response_type='response_hi_res_90_i3dr+1_n-3_noDCTDC2_2D_wall' 
#response_type='response_hi_res_90_i3dr+1_n-3_noDCTDC2'  #~20% max loss - recheck
#response_type='response_hi_res_90_i3dr+1_n-3_2D_wall'
#response_type='response_hi_res_90_i3dr-1_n-3_2D_wall'

response_type='response_hi_res_90_i3dr-1_n-3_hiPitchres_noDCTDC2_new_profiles_2D_wall'
response_type='response_hi_res_90_i3dr-1_n-3_hiPitchres_noDCTDC2_new_profiles' #new best fit

# notes
# BIG difference in 2D dfn between 3D and 2D wall...much higher 2D density with 3D wall

folder_2D='2D_23ms'  
folder_3D='3D_23ms'
folder_2D='2D_full_slow'  
folder_3D='3D_full_slow'

DFN_2D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST'/ response_type / folder_2D).glob('*.dfn'))[0]
DFN_3D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST'/ response_type / folder_3D).glob('*.dfn'))[0]
DFN_2D_split_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST' / response_type / folder_2D).glob('ptcl_cache.dat'))[0]
DFN_3D_split_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST' / response_type / folder_3D).glob('ptcl_cache.dat'))[0]
MOM_2D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST' / response_type / folder_2D).glob('*.h5'))[0]
MOM_3D_filename=list(pathlib.Path(support.dir_output_files / shot_number / 'LOCUST' / response_type / folder_3D).glob('*.h5'))[0]
eq_filename=support.dir_input_files / shot_number / 'LOCUST' / 'g157418.03000'
wall_filename=eq_filename

DFN_2D=DFN(ID='2D',data_format='LOCUST',filename=DFN_2D_filename)
DFN_3D=DFN(ID='3D',data_format='LOCUST',filename=DFN_3D_filename)
DFN_2D['E']/=1000. #convert to keV
DFN_3D['E']/=1000.
#XXX commented out whilst FPL is too big DFN_2D_split=FPL(ID='2D -DSPLIT',data_format='LOCUST',filename=DFN_2D_split_filename,compression=True,coordinates=['R','Z','phi','V_R','V_tor','V_Z','status_flag','time','PFC_intercept','additional_flag1'])
#XXX commented out whilst FPL is too big DFN_3D_split=FPL(ID='3D -DSPLIT',data_format='LOCUST',filename=DFN_3D_split_filename,compression=True,coordinates=['R','Z','phi','V_R','V_tor','V_Z','status_flag','time','PFC_intercept','additional_flag1'])
MOM_2D=Mom(ID='2D moments',data_format='LOCUST',filename=MOM_2D_filename)
MOM_3D=Mom(ID='3D moments',data_format='LOCUST',filename=MOM_3D_filename)
eq=EQ(ID='157418 300ms equilibrium',data_format='GEQDSK',filename=eq_filename,GEQDSKFIX=1)
wall=Wall(ID='157418 2D wall',data_format='GEQDSK',filename=wall_filename)

#crop DFN to within rho=sqrt(toroidal_flux)=0.7
#'''
rho_crop=0.8  #
rho_crop=0.73 #point at which new DFN grid discretisation cannot be blamed
rho_crop=0.79
rho_crop=0.78 
rho_crop=0.76
rho_crop=0.77 #best fit 
eq.set(flux_tor_norm=(eq['flux_tor']-eq['flux_tor'][0])/(eq['flux_tor'][-1]-eq['flux_tor'][0]))
eq.set(flux_tor_norm_rz=processing.utils.flux_func_to_RZ(eq['flux_pol'],np.sqrt(np.abs(eq['flux_tor_norm'])),eq)) #calculate rho grid
#eq.set(flux_tor_norm_rz=processing.utils.LCFS_crop(eq['flux_tor_norm_rz'],eq,eq,crop_value=0.0,outside=True)) #remove DFN outisde LCFS
flux_tor_norm_rz=processing.utils.value_at_RZ(*[dimension.T.flatten() for dimension in np.meshgrid(DFN_2D['R'],DFN_2D['Z'])],eq['flux_tor_norm_rz'],eq).reshape(len(DFN_2D['R']),len(DFN_2D['Z']))
DFN_2D.set(flux_tor_norm_rz=flux_tor_norm_rz)
#DFN_2D.plot(axes=['R','Z'],key='flux_tor_norm_rz',vminmax=[0.7,1.],real_scale=True,LCFS=eq) #XXX
DFN_3D.set(flux_tor_norm_rz=flux_tor_norm_rz)
DFN_2D['dfn'][...,flux_tor_norm_rz<rho_crop]=0.
DFN_3D['dfn'][...,flux_tor_norm_rz<rho_crop]=0.
DFN_2D['dfn'][:,DFN_2D['E']<10.,:,:,:]=0. #cut off DFN below 10keV
DFN_3D['dfn'][:,DFN_3D['E']<10.,:,:,:]=0.
DFN_2D['dfn'][:,DFN_2D['E']>80.5,:,:,:]=0. #cut off DFN above max SPIRAL energy 80.5keV (and seeming cut-off in data at ~82.5keV)
DFN_3D['dfn'][:,DFN_3D['E']>80.5,:,:,:]=0.
sav['surfhist_wc'][:,sav['ehis']>80.5]=0.
sav['surfhist_noc'][:,sav['ehis']>80.5]=0.
#'''

fig,((ax1))=plt.subplots(1) #plot the difference

vminmax=[0,1]
number_bins_diff=7#8

#compare with Mike's data
X=copy.deepcopy(sav['ehis'])
Y=copy.deepcopy(sav['phis'])
Z=copy.deepcopy(sav['surfhist_wc']-sav['surfhist_noc'])
#Z/=np.mean(sav['surfhist_noc']) #normalise
Z=(Z-np.max(Z))/(np.min(Z)-np.max(Z)) #normalise
#Z/=np.sum(Z)*(sav['ehis'][1]-sav['ehis'][0])*(sav['phis'][1]-sav['phis'][0]) #normalise
X,Y=np.meshgrid(X,Y) #dfn is r,z so need to swap order here
#mesh=ax1.contourf(X,Y,Z,levels=np.linspace(vminmax[0],vminmax[1],num=number_bins_diff),colors=settings.cmap_plasma(np.linspace(vminmax[0],vminmax[1],num=number_bins_diff)),edgecolor='none',linewidth=0,antialiased=True,vmin=vminmax[0],vmax=vminmax[1])
mesh=ax1.contourf(X,Y,Z,cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True,levels=np.linspace(vminmax[0],vminmax[1],num=number_bins_diff))
SPIRAL_contours,_ = mesh.legend_elements()
axes=['R','Z']
axes=['E','V_pitch']
DFN_diff=copy.deepcopy(DFN_3D)
DFN_diff['dfn']=(DFN_3D['dfn']-DFN_2D['dfn'])
print(f"total difference in FI fraction LOCUST ={100*(DFN_3D.transform(axes=['N'])['dfn']-DFN_2D.transform(axes=['N'])['dfn'])/DFN_2D.transform(axes=['N'])['dfn']}")
print(f"total difference in FI fraction SPIRAL ={100*(np.sum(sav['surfhist_wc'])-np.sum(sav['surfhist_noc']))/(np.sum(sav['surfhist_noc']))}")

#DFN_diff['dfn']/=DFN_diff.transform(axes=['N'])['dfn'] #normalise
DFN_diff=DFN_diff.transform(axes=axes)

#DFN_diff['dfn']/=np.mean(DFN_2D.transform(axes=axes)['dfn']) #normalise
DFN_diff['dfn']=(DFN_diff['dfn']-np.max(DFN_diff['dfn']))/(np.min(DFN_diff['dfn'])-np.max(DFN_diff['dfn'])) #normalise
#eq.plot(ax=ax1,fig=fig,key='flux_tor_norm_rz',vminmax=[0.7,1.],fill=False,number_bins=2) #XXX
DFN_diff_mesh=DFN_diff.plot(ax=ax1,fig=fig,axes=axes,colmap=settings.cmap_w,fill=False,number_bins=number_bins_diff,vminmax=vminmax,transform=False)
LOCUST_contours,_ = DFN_diff_mesh.legend_elements()    
cbar=fig.colorbar(mesh,orientation='vertical')
#for t in cbar.ax.get_yticklabels():
#     t.set_fontsize(20)
cbar.add_lines(DFN_diff_mesh)
ax1.set_xlim(10,80)
ax1.set_ylim(-1,1)
ax1.set_xlabel('Energy [keV]')
ax1.set_ylabel(r'Pitch')
ax1.set_xticks(np.linspace(10,80,8))
ax1.set_yticks(np.linspace(-1.,1.,5))
ax1.legend([LOCUST_contours[0],SPIRAL_contours[-1]],['LOCUST','SPIRAL'],loc='center right',frameon=False,labelcolor='w')
ax1.set_title(r'Normalised fast-ion density reduction due to $n$=3 RMP [a.u]',fontsize=20,pad=15)
plt.show()

#################################
 
##################################################################
 
###################################################################################################
