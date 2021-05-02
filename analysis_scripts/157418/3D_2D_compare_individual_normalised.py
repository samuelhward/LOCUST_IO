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

response_type='response'    
#response_type='vacuum'
#response_type='response_hi_res'
response_type='response_hi_res_90_i3dr-1'
#response_type='response_hi_res_-90_i3dr1'
#response_type='response_hi_res_90_i3dr-1_n-3_noDCTDC2' #XXX latest re-runs
##response_type='response_hi_res_90_i3dr+1_n-3_noDCTDC2' #XXX latest re-runs
#response_type='response_hi_res_90_i3dr+1_n-3' #XXX currently running on 8

response_type='response_hi_res_90_i3dr-1_n-3_noDCTDC2' #latest runs
#response_type='response_hi_res_90_i3dr+1_n-3_noDCTDC2'


#folder_2D='2D_23ms'  
#folder_3D='3D_23ms'
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
#XXX commented out whilst FPL is too big DFN_2D_split=FPL(ID='2D -DSPLIT',data_format='LOCUST',filename=DFN_2D_split_filename,compression=True,coordinates=['R','Z','phi','V_R','V_tor','V_Z','status_flag','time','PFC_intercept','additional_flag1'])
#XXX commented out whilst FPL is too big DFN_3D_split=FPL(ID='3D -DSPLIT',data_format='LOCUST',filename=DFN_3D_split_filename,compression=True,coordinates=['R','Z','phi','V_R','V_tor','V_Z','status_flag','time','PFC_intercept','additional_flag1'])
MOM_2D=Mom(ID='2D moments',data_format='LOCUST',filename=MOM_2D_filename)
MOM_3D=Mom(ID='3D moments',data_format='LOCUST',filename=MOM_3D_filename)
eq=EQ(ID='157418 300ms equilibrium',data_format='GEQDSK',filename=eq_filename,GEQDSKFIX=1)
wall=Wall(ID='157418 2D wall',data_format='GEQDSK',filename=wall_filename)


#calculate some missing quantities and edit the data
#XXX commented out whilst FPL is too big DFN_2D_split['E']=.5*constants.species_mass*(DFN_2D_split['V_R']**2+DFN_2D_split['V_tor']**2+DFN_2D_split['V_Z']**2)/constants.charge_e
#XXX commented out whilst FPL is too big DFN_3D_split['E']=.5*constants.species_mass*(DFN_3D_split['V_R']**2+DFN_3D_split['V_tor']**2+DFN_3D_split['V_Z']**2)/constants.charge_e
#DFN_2D_split['V_pitch']=processing.utils.pitch_calc_2D(DFN_2D_split,eq)
#DFN_3D_split['V_pitch']=processing.utils.pitch_calc_2D(DFN_3D_split,eq)


#these status flags count as a loss in LOCUST

'''
         if( nint(ph(i,8))== -1 .or.                                          &
             nint(ph(i,8))== -3 .or.                                          &
             nint(ph(i,8))== -4 .or.                                          &
             nint(ph(i,8))== -5 .or.                                          &
             nint(ph(i,8))== -6 .or.                                          &
             nint(ph(i,8))== -8 .or.                                          &
             nint(ph(i,8))== -9 .or.                                          &
             nint(ph(i,8))==-10 .or.                                          &
             nint(ph(i,8))==-11 .or.                                          &
             nint(ph(i,8))==-12 .or.                                          &
             nint(ph(i,8))==-14 .or.                                          &
             nint(ph(i,8))==-15 )then
'''

#need to scale up weights of particle list to match 1W injected power
#XXX commented out whilst FPL is too big i=np.where((DFN_2D_split['status_flag']<=0)&(DFN_2D_split['status_flag']>=-15.))[0] #calculate the PFC power flux
#i=np.where((DFN_2D_split['status_flag']==-5) | (DFN_2D_split['status_flag']==-8.) | (DFN_2D_split['status_flag']==-9.) | (DFN_2D_split['status_flag']==-14.) | (DFN_2D_split['status_flag']==0.))[0] #calculate the PFC power flux
#XXX commented out whilst FPL is too big number_of_particles_2D=len(DFN_2D_split['status_flag'][i]) 
#XXX commented out whilst FPL is too big i=np.where((DFN_3D_split['status_flag']<=0)&(DFN_3D_split['status_flag']>=-15.))[0]
#i=np.where((DFN_3D_split['status_flag']==-5) | (DFN_3D_split['status_flag']==-8.) | (DFN_3D_split['status_flag']==-9.) | (DFN_3D_split['status_flag']==-14.) | (DFN_3D_split['status_flag']==0.))[0] #calculate the PFC power flux
#XXX commented out whilst FPL is too big number_of_particles_3D=len(DFN_3D_split['status_flag'][i])
#XXX commented out whilst FPL is too big Pdep_2D_desired=1.
#XXX commented out whilst FPL is too big Pdep_3D_desired=1.
#XXX commented out whilst FPL is too big Einj_2D=80000
#XXX commented out whilst FPL is too big Einj_3D=80000
#XXX commented out whilst FPL is too big Pdep_2D_simulation=Einj_2D*number_of_particles_2D*constants.charge_e
#XXX commented out whilst FPL is too big Pdep_3D_simulation=Einj_3D*number_of_particles_3D*constants.charge_e
#XXX commented out whilst FPL is too big Pdep_scale_factor_2D=Pdep_2D_desired/Pdep_2D_simulation
#XXX commented out whilst FPL is too big Pdep_scale_factor_3D=Pdep_3D_desired/Pdep_3D_simulation

#XXX commented out whilst FPL is too big i=np.where((DFN_2D_split['PFC_intercept']<0)) #look at all markers which were lost from plasma
#XXX commented out whilst FPL is too big DFN_2D_split['PFC_power']=np.histogram(DFN_2D_split['time'][i],bins=100,weights=DFN_2D_split['E'][i]*constants.charge_e*Pdep_scale_factor_2D)[0] #place against same time base as LOCUST
#XXX commented out whilst FPL is too big DFN_2D_split['PFC_power']=np.cumsum(DFN_2D_split['PFC_power'])
#XXX commented out whilst FPL is too big i=np.where((DFN_3D_split['PFC_intercept']<0))
#XXX commented out whilst FPL is too big DFN_3D_split['PFC_power']=np.histogram(DFN_3D_split['time'][i],bins=100,weights=DFN_3D_split['E'][i]*constants.charge_e*Pdep_scale_factor_3D)[0] #place against same time base as LOCUST
#XXX commented out whilst FPL is too big DFN_3D_split['PFC_power']=np.cumsum(DFN_3D_split['PFC_power'])

#XXX commented out whilst FPL is too big DFN_2D_split['weight']=DFN_2D_split['E']*constants.charge_e*Pdep_scale_factor_2D
#XXX commented out whilst FPL is too big DFN_3D_split['weight']=DFN_3D_split['E']*constants.charge_e*Pdep_scale_factor_3D

#crop DFN to within rho=sqrt(toroidal_flux)=0.7
#'''
rho_crop=0.7
eq.set(flux_tor_norm=(eq['flux_tor']-eq['flux_tor'][0])/(eq['flux_tor'][-1]-eq['flux_tor'][0]))
eq.set(flux_tor_norm_rz=processing.utils.flux_func_to_RZ(eq['flux_pol'],np.sqrt(np.abs(eq['flux_tor_norm'])),eq)) #calculate rho grid
eq.set(flux_tor_norm_rz=processing.utils.LCFS_crop(eq['flux_tor_norm_rz'],eq,eq,crop_value=0.0,outside=True))
#eq.plot(key='flux_tor_norm_rz',vminmax=[0.7,1.],fill=False,number_bins=2) #XXX
flux_tor_norm_rz=processing.utils.value_at_RZ(*[dimension.T.flatten() for dimension in np.meshgrid(DFN_2D['R'],DFN_2D['Z'])],eq['flux_tor_norm_rz'],eq).reshape(len(DFN_2D['R']),len(DFN_2D['Z']))
DFN_2D.set(flux_tor_norm_rz=flux_tor_norm_rz)
#DFN_2D.plot(axes=['R','Z'],key='flux_tor_norm_rz',vminmax=[0.7,1.],real_scale=True,LCFS=eq) #XXX
DFN_3D.set(flux_tor_norm_rz=flux_tor_norm_rz)
DFN_2D['dfn'][...,flux_tor_norm_rz<rho_crop]=0.
DFN_3D['dfn'][...,flux_tor_norm_rz<rho_crop]=0.
DFN_2D['dfn'][:,DFN_2D['E']<1.e4,:,:,:]=0. #cut off DFN below 10keV
DFN_3D['dfn'][:,DFN_3D['E']<1.e4,:,:,:]=0.
#'''



fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)

'''
mesh_2D=DFN_2D.plot(fig=fig,ax=ax1,axes=['E','V_pitch'],vminmax=[0,30000.])
mesh_3D=DFN_3D.plot(fig=fig,ax=ax2,axes=['E','V_pitch'],vminmax=[0,30000.])
ax1.set_xlim(10000,85000)
ax1.set_ylim(-1,1)
ax2.set_xlim(10000,85000)
ax2.set_ylim(-1,1)
'''

DFN_2D.plot(fig=fig,ax=ax1,axes=['R','Z'],real_scale=True,limiters=wall)
DFN_3D.plot(fig=fig,ax=ax2,axes=['R','Z'],real_scale=True,limiters=wall)
#for ax,mesh in zip([ax1,ax2],[mesh_2D,mesh_3D]):
#    cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')


#DFN_3D_split.plot(fig=fig,ax=ax3,axes=['time'],status_flags=['PFC_intercept_3D'],colmap=settings.cmap_k,weight=True,number_bins=300)#,style='scatter') #as a function of energy
#DFN_2D_split.plot(fig=fig,ax=ax3,axes=['time'],status_flags=['PFC_intercept_3D'],colmap=settings.cmap_b,weight=True,number_bins=300)#,style='scatter') 
#DFN_2D_split.plot(fig=fig,ax=ax3,axes=['R','Z'],status_flags=['PFC_intercept_3D'],colfield='E',colmap=cmap_blues,weight=True,style='scatter',real_scale=True,limiters=wall) #scatter in RZ space 
#DFN_3D_split.plot(fig=fig,ax=ax3,axes=['R','Z'],status_flags=['PFC_intercept_3D'],colfield='E',colmap=cmap_reds,weight=True,style='scatter',real_scale=True,limiters=wall) 
#XXX commented out whilst FPL is too big DFN_2D_split.plot(fig=fig,ax=ax3,axes=['R','Z'],status_flags=['PFC_intercept_3D'],weight=True,style='histogram',real_scale=True,number_bins=500,limiters=wall) #heat in RZ space 
#XXX commented out whilst FPL is too big DFN_3D_split.plot(fig=fig,ax=ax3,axes=['R','Z'],status_flags=['PFC_intercept_3D'],weight=True,style='histogram',real_scale=True,number_bins=500,limiters=wall) 
#DFN_3D_split.plot(fig=fig,ax=ax3,axes=['phi','Z'],status_flags=['PFC_intercept_3D'],colfield='E',colmap=cmap_reds,weight=True,style='scatter') 
#DFN_2D_split.plot(fig=fig,ax=ax3,axes=['phi','Z'],status_flags=['PFC_intercept_3D'],colfield='E',colmap=cmap_blues,weight=True,style='scatter') #in Z phi space 
#DFN_2D_split.plot(fig=fig,ax=ax3,axes=['R'],status_flags=['PFC_intercept_3D'],colmap=settings.cmap_g,weight=True,style='histogram',number_bins=100) #heat in R space 
#DFN_3D_split.plot(fig=fig,ax=ax3,axes=['R'],status_flags=['PFC_intercept_3D'],colmap=settings.cmap_r,weight=True,style='histogram',number_bins=100) 


#moment_to_plot='density'
#MOM_2D.plot(key=moment_to_plot,fig=fig,ax=ax4,colmap=settings.cmap_b)
#MOM_3D.plot(key=moment_to_plot,fig=fig,ax=ax4,colmap=settings.cmap_k)
#ax4.set_xlim([0.6,1.1])
#ax4.set_ylim([0.,8e11])
#plt.show()

#XXX commented out whilst FPL is too big plt.plot(100*DFN_2D_split['PFC_power']/Pdep_2D_desired,color='b') #100 to turn into percentage
#XXX commented out whilst FPL is too big plt.plot(100*DFN_3D_split['PFC_power']/Pdep_3D_desired,color='k')
#plt.title('% PFC power/Pdep')
#plt.legend(('2D','3D'))
#plt.show()


# read SPIRAL data
filepath_sav=support.dir_output_files / shot_number / 'SPIRAL_output' / 'all_variables_for_sam.sav'
sav=readsav(filepath_sav)

# check SPIRAL data
fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2) #plot the difference
mesh1=ax1.contourf(sav['ehis'],sav['phis'],sav['surfhist_wc']-sav['surfhist_noc'],cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True)
cbar1=fig.colorbar(mesh1,ax=ax1,orientation='vertical')
mesh2=ax2.contourf(sav['ehis'],sav['phis'],sav['surfhist_noc'],cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True)
cbar2=fig.colorbar(mesh2,ax=ax2,orientation='vertical')

# check LOCUST data

DFN_diff=copy.deepcopy(DFN_3D)
DFN_diff['dfn']=DFN_3D['dfn']-DFN_2D['dfn']
DFN_diff_mesh=DFN_diff.plot(ax=ax3,fig=fig,axes=['E','V_pitch'],colmap=settings.cmap_inferno,fill=True,number_bins=10)
DFN_2D_mesh=DFN_2D.plot(ax=ax4,fig=fig,axes=['E','V_pitch'],colmap=settings.cmap_inferno,fill=True,number_bins=10)
cbar3=fig.colorbar(DFN_diff_mesh,ax=ax3,orientation='vertical')
cbar4=fig.colorbar(DFN_2D_mesh,ax=ax4,orientation='vertical')
plt.show()

# check LOCUST+SPIRAL data
fig,ax=plt.subplots(1)
mesh1=ax1.contour(sav['ehis'],sav['phis'],sav['surfhist_noc'],cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True,levels=5)
DFN_2D_mesh=DFN_2D.plot(ax=ax,fig=fig,axes=['E','V_pitch'],colmap=settings.cmap_inferno,fill=False,number_bins=5)
ax.set_xlim(10,80)
plt.show()


#remove silly numbers from array
def remove_outliers(array):
    array[array==np.inf]=1.e2
    array[array==-np.inf]=1.e2
    array[array>1.e10]=1.e2
    array[array<-1.e10]=1.e2
    return array

fig,((ax1))=plt.subplots(1) #plot the difference

vminmax=[-1.5,0]
number_bins_diff=10

#compare with Mike's data
X=copy.deepcopy(sav['ehis'])
Y=copy.deepcopy(sav['phis'])
Z=np.nan_to_num(np.abs(sav['surfhist_wc']-sav['surfhist_noc'])/sav['surfhist_noc']) #normalise
print(f"spiral total fractional change = {np.sum(np.abs(sav['surfhist_wc']-sav['surfhist_noc']))/np.sum(sav['surfhist_noc'])}")
Z=remove_outliers(Z)
Z=np.log10(Z)
Z=remove_outliers(Z)
print(f'spiral max = {np.max(Z)}')
print(f'spiral min = {np.min(Z)}')
print(f'spiral mean = {np.mean(Z)}')#print(np.mean(np.nan_to_num(np.abs(sav['surfhist_wc']-sav['surfhist_noc'])/sav['surfhist_noc'],nan=1e-5)))
#Z=(Z-np.max(Z))/(np.min(Z)-np.max(Z)) #normalise
#Z/=np.sum(Z)*(sav['ehis'][1]-sav['ehis'][0])*(sav['phis'][1]-sav['phis'][0]) #normalise
X,Y=np.meshgrid(X,Y) #dfn is r,z so need to swap order here
#mesh=ax1.contourf(X,Y,Z,levels=np.linspace(vminmax[0],vminmax[1],num=number_bins_diff),colors=settings.cmap_plasma(np.linspace(vminmax[0],vminmax[1],num=number_bins_diff)),edgecolor='none',linewidth=0,antialiased=True,vmin=vminmax[0],vmax=vminmax[1])
mesh=ax1.contourf(X,Y,Z,cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True,levels=np.linspace(vminmax[0],vminmax[1],num=number_bins_diff))
#mesh=ax1.contourf(X,Y,Z,cmap=settings.cmap_inferno,edgecolor='none',linewidth=0,antialiased=True)
SPIRAL_contours,_ = mesh.legend_elements()
axes=['R','Z']
axes=['E','V_pitch']
DFN_diff=copy.deepcopy(DFN_3D)
DFN_diff['dfn']=DFN_3D['dfn']-DFN_2D['dfn']
print(f"locust total fractional change = {np.sum(np.abs(DFN_3D['dfn']-DFN_2D['dfn']))/np.sum(DFN_2D['dfn'])}")
#DFN_diff['dfn']/=DFN_diff.transform(axes=['N'])['dfn'] #normalise
DFN_diff=DFN_diff.transform(axes=axes)
DFN_diff['dfn']=np.nan_to_num(np.abs(DFN_diff['dfn'])/DFN_2D.transform(axes=axes)['dfn']) #normalise
DFN_diff['dfn']=remove_outliers(DFN_diff['dfn'])
DFN_diff['dfn']=np.log10(DFN_diff['dfn'])
DFN_diff['dfn']=remove_outliers(DFN_diff['dfn'])
print(f"locust max = {np.max(DFN_diff['dfn'])}")
print(f"locust min = {np.min(DFN_diff['dfn'])}")
print(f"locust mean = {np.mean(DFN_diff['dfn'])}")
#import sys
#sys.exit()
#DFN_diff['dfn']=(DFN_diff['dfn']-np.max(DFN_diff['dfn']))/(np.min(DFN_diff['dfn'])-np.max(DFN_diff['dfn'])) #normalise
DFN_diff['E']/=1000 #convert to keV for plotting
#eq.plot(ax=ax1,fig=fig,key='flux_tor_norm_rz',vminmax=[0.7,1.],fill=False,number_bins=2) #XXX
DFN_diff_mesh=DFN_diff.plot(ax=ax1,fig=fig,axes=axes,colmap=settings.cmap_inferno,fill=False,number_bins=number_bins_diff,vminmax=vminmax,transform=False)
#DFN_diff_mesh=DFN_diff.plot(ax=ax1,fig=fig,axes=axes,colmap=settings.cmap_inferno,fill=False,transform=False)
LOCUST_contours,_ = DFN_diff_mesh.legend_elements()    
cbar=fig.colorbar(mesh,ax=ax1,orientation='vertical')
cbar.add_lines(DFN_diff_mesh)
ax1.set_xlim(10,80)
ax1.set_ylim(-1,1)
ax1.set_xlabel('Energy [keV]')
ax1.set_ylabel(r'Pitch')
ax1.set_xticks(np.linspace(10,80,8))
ax1.set_yticks(np.linspace(-1.,1.,5))
ax1.legend([LOCUST_contours[0],SPIRAL_contours[-1]],['LOCUST','SPIRAL'],fontsize=20,loc='center right',framealpha=0.4)
ax1.set_title(r'Fast-ion density reduction due to $n$=3 magnetic perturbation')
fig.set_tight_layout(True)  
plt.show()

axis='V_pitch'
axis='E'
if axis=='V_pitch':
    axis_to_integrate=1
    quantity_to_integrate='ehis'
    quantity_to_plot_against='phis'
    convert_to_eV=1000.
elif axis=='E':
    axis_to_integrate=0
    quantity_to_integrate='phis'
    quantity_to_plot_against='ehis'
    convert_to_eV=1.

fig,((ax1))=plt.subplots(1) #plot the difference in E space
#Z=copy.deepcopy(sav['surfhist_wc']-sav['surfhist_noc'])

surfhist_wc_summed=np.sum(sav['surfhist_wc'],axis=axis_to_integrate)
surfhist_noc_summed=np.sum(sav['surfhist_noc'],axis=axis_to_integrate)
Z=np.nan_to_num(np.abs(surfhist_wc_summed-surfhist_noc_summed)/surfhist_noc_summed) #normalise
Z=remove_outliers(Z)
Z=np.log10(Z)
Z=remove_outliers(Z)
#Z=(Z-np.max(Z))/(np.min(Z)-np.max(Z)) #normalise
#Z/=np.sum(Z)*(sav['ehis'][1]-sav['ehis'][0])*(sav['phis'][1]-sav['phis'][0])*convert_to_eV #normalise
ax1.plot(sav[quantity_to_plot_against],Z,color='red')

DFN_diff=copy.deepcopy(DFN_3D)
DFN_diff['dfn']=DFN_3D.transform(axes=[axis])['dfn']-DFN_2D.transform(axes=[axis])['dfn']
DFN_diff['dfn']=np.nan_to_num(np.abs(DFN_diff['dfn'])/DFN_2D.transform(axes=[axis])['dfn']) #normalise
DFN_diff['dfn']=remove_outliers(DFN_diff['dfn'])
DFN_diff['dfn']=np.log10(DFN_diff['dfn'])
DFN_diff['dfn']=remove_outliers(DFN_diff['dfn'])
DFN_diff['E']/=1000 #convert to keV for plotting
#DFN_diff['dfn']=(DFN_diff['dfn']-np.max(DFN_diff['dfn']))/(np.min(DFN_diff['dfn'])-np.max(DFN_diff['dfn'])) #normalise
DFN_diff_mesh=DFN_diff.plot(ax=ax1,fig=fig,axes=[axis],colmap=settings.cmap_k,transform=False)

plt.title('')
plt.show()

#################################
 
##################################################################
 
###################################################################################################
