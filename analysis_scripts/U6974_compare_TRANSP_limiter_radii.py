#compare multiple TRANSP dfns with the different limiter profiles - 3 x N grid, one row with RZ, one with EP and one EP but sampled

import sys
import numpy as np
import context
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.wall import Wall
from classes.output_classes.distribution_function import Distribution_Function
import run_scripts.utils

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
TRANSP_runs=['U69','U70','U71','U72','U73','U74']
TRANSP_ID='157418'
wall_file='LCFS_DIII-D.dat_'
tail='_fi_1_gc.cdf'

#start by creating some axis objects
fig,((ax1,ax2,ax3,ax4,ax5,ax6),(ax7,ax8,ax9,ax10,ax11,ax12),(ax13,ax14,ax15,ax16,ax17,ax18))=plt.subplots(3,6)

axes=[[ax1,ax2,ax3,ax4,ax5,ax6],[ax7,ax8,ax9,ax10,ax11,ax12],[ax13,ax14,ax15,ax16,ax17,ax18]]
axes=list(map(list, zip(*axes)))

r=2.3
z=0.

for radius,TRANSP_run,axis_row in zip(radii,TRANSP_runs,axes):

    wall=Wall(radius,'LOCUST_2D',wall_file+radius)

    TRANSP_file=TRANSP_ID+TRANSP_run+tail    
    TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density (0.1s) - limiter radius = {}'.format(str(radius)),filename=TRANSP_file)

    #RZ

    axes=['R','Z']
    TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,some_equilibrium=wall,limiters=True,ax=axis_row[0],fig=fig)

    #EP

    axes=['E','V_pitch']
    TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,some_equilibrium=wall,limiters=True,ax=axis_row[1],fig=fig)

    #sampled EP

    #for point we have chosen, find closest dfn bin in real space
    diff_r=np.abs(TRANSP_dfn['R2D']-r)**2 #irregular grid so need to minimise distance to every point
    diff_z=np.abs(TRANSP_dfn['Z2D']-z)**2
    index_rz_transp=(diff_r+diff_z).argmin()
    r_sample=TRANSP_dfn['R2D'][index_rz_transp]
    z_sample=TRANSP_dfn['Z2D'][index_rz_transp]

    axes=[index_rz_transp,slice(None),slice(None)]
    TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,some_equilibrium=wall,limiters=True,ax=axis_row[2],fig=fig)

    #point we're sampling from

    axis_row[0].scatter(r_sample,z_sample,color='m',s=0.5)


    #set titles and axes labels
    #axis_row[0].set_title('') #keep as file ID - can act as row title
    axis_row[1].set_title('')
    axis_row[2].set_title('sample point - r={r}z={z}'.format(r=r_sample,z=z_sample))

    axis_row[0].set_xlabel('R [m]')
    axis_row[0].set_ylabel('Z [m]')
    axis_row[1].set_xlabel('Energy [eV]')
    axis_row[1].set_ylabel('Pitch')
    axis_row[2].set_xlabel('Energy [eV]')
    axis_row[2].set_ylabel('Pitch')


plt.show()
