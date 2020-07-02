#plot_FIDA.py
 
"""
Samuel Ward
11/03/2020
----
reads FIDA data from fida_gui.pro and plots
---
usage:

---
"""

def get_FIDA(filepath):
    """
    read FIDA data from fida_gui.pro

    notes:
    args:
        filepath - 
    """
    
    with open(filepath,'r') as file:

        data={}

        lines=file.readlines()
        del(lines[0]) #delete header

        while True:
            if ';' not in lines[1]:
                data[str(lines[0].split()[1:])]=lines[0].split()[0]
                del(lines[0])
            else: 
                del(lines[0])
                break

        column_names=[]
        for column_name in [column.strip() for column in lines[0].split(';')]:
            data[column_name]=[]
            column_names.append(column_name)
        del(lines[0])

        for line in lines:
            split_line=line.split()
            for counter,column in enumerate(split_line):
                data[column_names[counter]].append(float(column))
        for key in data:
            data[key]=np.asarray(data[key])
    
    return data


import context
import support
import pathlib
import settings
import numpy as np
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium as equi

#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])
cmap_grey=settings.colour_custom([97,97,97,1])

channel=5#5 #4567
filepath_measured=support.dir_output_files / 'FIDASIM' / f'29034_chn_{channel}.dat'
filepath_FIDASIM=support.dir_output_files / 'FIDASIM' / f'F_05-11-2019_23-14-13_chn_{channel}.dat'

data_FIDASIM=get_FIDA(filepath_FIDASIM)
data_measured=get_FIDA(filepath_measured)

fig,axes=plt.subplots(1,2)
sum_signals=np.zeros(len(data_FIDASIM['Wavelength [nm]']))
for key,value in data_FIDASIM.items():
    if value.ndim>=1 and 'Wavelength' not in key:
        axes[0].plot(data_FIDASIM['Wavelength [nm]'],value,label=f'{key}',linewidth=settings.plot_linewidth,zorder=0)
        sum_signals+=value
        print(key)
axes[0].plot(data_FIDASIM['Wavelength [nm]'],sum_signals,linestyle='--',label='sum',linewidth=settings.plot_linewidth,zorder=0)
axes[0].errorbar(data_measured['Wavelength [nm]'],data_measured['Radiance [photons/(s nm m^2 sr)]'],data_measured['Uncertainty [photons/(s nm m^2 sr)]'],label='measurements',fmt='.',color=cmap_grey(0.0),linewidth=settings.plot_linewidth,zorder=10)
axes[0].set_xlim([657,663])
axes[0].set_ylim([1.e15,1.e19])
axes[0].set_xlabel('Wavelength [nm]',fontsize=25)
axes[0].set_ylabel('Radiance [photons/(s nm m^2 sr)]',fontsize=25)
axes[0].set_title('')
axes[0].axvline(660.7,color=cmap_grey(0.0))
axes[0].axvline(661.5,color=cmap_grey(0.0))
axes[0].set_yscale('log')
axes[0].legend(fontsize=20)

file_eq=pathlib.Path('TRANSP') / '29034W04' / 'g29034' #grab the wall and equilibrium used in these simulations
eq=equi('',data_format='GEQDSK',filename=file_eq)

filepath_radial_profile_TRANSP=support.dir_output_files / 'FIDASIM' / '29034_W04_fi_10_radial_profile.dat'
filepath_radial_profile_LOCUST=support.dir_output_files / 'FIDASIM' / '29034_F_05-11-2019_23-14-13_radial_profile.dat'

data_radial_TRANSP=get_FIDA(filepath_radial_profile_TRANSP)
data_radial_LOCUST=get_FIDA(filepath_radial_profile_LOCUST)

axes[1].plot(data_radial_TRANSP['Radius [m]'],data_radial_TRANSP['FIDASIM integrated intensity [photons/(s m^2 sr)]'],label='TRANSP',color=cmap_r(0),linewidth=settings.plot_linewidth,zorder=0)
axes[1].plot(data_radial_LOCUST['Radius [m]'],data_radial_LOCUST['FIDASIM integrated intensity [photons/(s m^2 sr)]'],label='LOCUST',color=cmap_g(0),linewidth=settings.plot_linewidth,zorder=0)
axes[1].errorbar(data_radial_TRANSP['Radius [m]'],data_radial_TRANSP['Integrated intensity [photons/(s m^2 sr)]'],data_radial_TRANSP['Uncertainty [photons/(s m^2 sr)]'],label='measurements',fmt='.',color=cmap_grey(0),linewidth=settings.plot_linewidth,zorder=10)
axes[1].set_xlim([0.8,1.5])
axes[1].set_ylim([0,2.e16])
axes[1].set_xlabel('Radius [m]',fontsize=25)
axes[1].set_ylabel('Integrated intensity [photons/(s m^2 sr)]',fontsize=25)
axes[1].set_title('')
axes[1].legend(fontsize=20)
axes[1].axvline(np.max(eq['lcfs_r']),linewidth=settings.plot_linewidth,color=settings.plot_colour_LCFS)

plt.show()

#################################
 
##################################################################
 
###################################################################################################
