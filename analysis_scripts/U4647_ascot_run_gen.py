import context
from classes.input_classes.temperature import Temperature
from classes.input_classes.number_density import Number_Density
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition
import run_scripts.utils
import numpy as np

file_ti='profile_Ti.dat'
file_te='profile_Te.dat'

file_ni='profile_ne.dat' #assume same density since Zeff=1
file_ne='profile_ne.dat'

file_eq='g157418.03000'

file_bd='U4647_ptcles.dat' #these are generated in the LOCUSTU4647 directory via the LOCUST_IO run scripts in there (from concatenating the TRANSP birth CDF files)



Ti=Temperature(ID='ion temperature',data_format='LOCUST',filename=file_ti)
Te=Temperature(ID='electron temperature',data_format='LOCUST',filename=file_te)

ni=Number_Density(ID='ion density',data_format='LOCUST',filename=file_ni)
ne=Number_Density(ID='electron density',data_format='LOCUST',filename=file_ne)

eq=Equilibrium(ID='equilibrium',data_format='GEQDSK',filename=file_eq)

bd=Beam_Deposition(ID='beam deposition',data_format='LOCUST_weighted',filename=file_bd)

bd.set(X=bd['R']*np.cos(bd['phi']),Y=bd['R']*np.sin(bd['phi']))

#add extra bits onto the end of the kinetic profiles to stop any divide by zero errors - as recommended by Emmi
Ti['T']=np.append(Ti['T'],np.zeros(20)+1.e1)
Te['T']=np.append(Te['T'],np.zeros(20)+1.e1)
ni['n']=np.append(ni['n'],np.zeros(20)+1.e5)
ne['n']=np.append(ne['n'],np.zeros(20)+1.e5)
Ti['flux_pol_norm']=np.append(Ti['flux_pol_norm'],np.linspace(1.001,1.5,20)) #extend out to psi=1.5
Te['flux_pol_norm']=np.append(Te['flux_pol_norm'],np.linspace(1.001,1.5,20))
ni['flux_pol_norm']=np.append(ni['flux_pol_norm'],np.linspace(1.001,1.5,20))
ne['flux_pol_norm']=np.append(ne['flux_pol_norm'],np.linspace(1.001,1.5,20))

for quantity in bd.data.keys(): #reduce by factor of 10 to 500k
    if type(bd[quantity])!=float:
        bd[quantity]=bd[quantity][1::10]

rot=np.zeros(len(Ti['T']))

run_scripts.utils.ASCOT_run_gen(temperature_i=Ti,temperature_e=Te,density_i=ni,density_e=ne,rotation_toroidal=rot,equilibrium=eq,beam_deposition=bd,guiding_centre=True)