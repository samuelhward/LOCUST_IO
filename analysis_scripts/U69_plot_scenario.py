#plot_scenario

import context,pathlib
import settings
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.number_density import Number_Density
from classes.input_classes.temperature import Temperature
import matplotlib.pyplot as plt

filename_eq = pathlib.Path('',"g157418.03000")
eq = Equilibrium("", data_format="GEQDSK", filename=filename_eq, GEQDSKFIX=1)
filename_dens_e=pathlib.Path("profile_ne.dat")
dens_e=Number_Density('','LOCUST',filename_dens_e)
filename_temp_e=pathlib.Path("profile_Te.dat")
filename_temp_i=pathlib.Path("profile_Ti.dat")
temp_e=Temperature('','LOCUST',filename_temp_e)
temp_i=Temperature('','LOCUST',filename_temp_i)

fig,(ax0,ax)=plt.subplots(1,2)
eq.plot(ax=ax0,fig=fig)
axes=[ax]
position_y_ax=0
for colour,kinetic_profile,plotting_variable in zip([settings.cmap_g,settings.cmap_r,settings.cmap_b],[dens_e,temp_e,temp_i],['electron density','electron temperature','ion temperature']):
    axes.append(axes[0].twinx())
    axes[-1].tick_params(axis='y', labelcolor=colour(0.5))
    axes[-1].set_ylabel(plotting_variable,color=colour(0.5)) 
    axes[-1].spines['right'].set_position(('outward', position_y_ax))
    position_y_ax+=60      
    kinetic_profile.plot(label=kinetic_profile.ID,colmap=colour,ax=axes[-1],fig=fig)
    axes[-1].set_title('')
    axes[-1].legend()
del(axes[0])

plt.show()

