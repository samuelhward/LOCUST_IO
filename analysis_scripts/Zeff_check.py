import context
from classes.input_classes.number_density import Number_Density as nd
import matplotlib.pyplot as plt
import processing.utils
import pathlib

nd1=nd('','LOCUST',pathlib.Path('LOCUST')/'profile_ni1.dat',species='ions')
nd2=nd('','LOCUST',pathlib.Path('LOCUST')/'profile_ni2.dat',species='ions')


#ratio of densities should be constant for constant Zeff
print(nd1['n'][0:50]/nd2['n'][0:50])
plt.plot(nd1['n'][0:50]/nd2['n'][0:50])
plt.show()

print(processing.utils.Zeff_calc(density=[nd1['n'][0],nd2['n'][0]],charge=[2,7]))