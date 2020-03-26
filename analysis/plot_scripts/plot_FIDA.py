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

#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])

channel=5#5 #4567
filepath_measured=support.dir_output_files / 'FIDASIM' / f'29034_chn_{channel}.dat'
filepath_FIDASIM=support.dir_output_files / 'FIDASIM' / f'F_05-11-2019_23-14-13_chn_{channel}.dat'

data_FIDASIM=get_FIDA(filepath_FIDASIM)
data_measured=get_FIDA(filepath_measured)

fig,ax=plt.subplots(1)
sum_signals=np.zeros(len(data_FIDASIM['Wavelength [nm]']))
for key,value in data_FIDASIM.items():
    if value.ndim>=1 and 'Wavelength' not in key:
        ax.plot(data_FIDASIM['Wavelength [nm]'],value,label=f'{key}')
        sum_signals+=value

        print(key)
ax.errorbar(data_measured['Wavelength [nm]'],data_measured['Radiance [photons/(s nm m^2 sr)]'],data_measured['Uncertainty [photons/(s nm m^2 sr)]'],label='measurements')
ax.plot(data_FIDASIM['Wavelength [nm]'],sum_signals,'.',label='sum')
ax.set_xlim([657,663])
ax.set_ylim([1.e15,1.e19])
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Radiance [photons/(s nm m^2 sr)]')
ax.set_title('')
plt.yscale('log')
ax.legend()
plt.show()




filepath_radial_profile_TRANSP=support.dir_output_files / 'FIDASIM' / '29034_W04_fi_10_radial_profile.dat'
filepath_radial_profile_LOCUST=support.dir_output_files / 'FIDASIM' / '29034_F_05-11-2019_23-14-13_radial_profile.dat'

data_radial_TRANSP=get_FIDA(filepath_radial_profile_TRANSP)
data_radial_LOCUST=get_FIDA(filepath_radial_profile_LOCUST)

fig,ax=plt.subplots(1)
ax.errorbar(data_radial_TRANSP['Radius [m]'],data_radial_TRANSP['Integrated intensity [photons/(s m^2 sr)]'],data_radial_TRANSP['Uncertainty [photons/(s m^2 sr)]'],label='measurements')
ax.plot(data_radial_TRANSP['Radius [m]'],data_radial_TRANSP['FIDASIM integrated intensity [photons/(s m^2 sr)]'],label='TRANSP')
ax.plot(data_radial_LOCUST['Radius [m]'],data_radial_LOCUST['FIDASIM integrated intensity [photons/(s m^2 sr)]'],label='LOCUST')
ax.set_xlim([0.8,1.5])
ax.set_ylim([0,2.e16])
ax.set_xlabel('Radius [m]')
ax.set_ylabel('Integrated intensity [photons/(s m^2 sr)]')
ax.set_title('')
ax.legend()
plt.show()

#################################
 
##################################################################
 
###################################################################################################
