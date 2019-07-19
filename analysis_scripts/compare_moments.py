#something to compare TRANSP and LOCUST moments

import context
from classes.output_classes.moments import Moments as mom 
import scipy
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
from mpl_toolkits import mplot3d #import 3D plotting axes
from mpl_toolkits.mplot3d import Axes3D

#U46
file_transp='157418U46.CDF'
file_locust='LOCUST_24-09-2018_11-55-41.801.h5'

#U48
file_transp='157418U48.CDF'
file_locust='LOCUST_24-09-2018_12-00-14.111.h5'

#U50
file_transp='157418U50.CDF'
file_locust='LOCUST_24-09-2018_11-55-49.392.h5'

mom_transp=mom(ID='TRANSP U46 moments - t = 0.1s',data_format='TRANSP',filename=file_transp)
mom_locust=mom(ID='LOCUST U46 moments - t = 0.1s',data_format='LOCUST',filename=file_locust)

locust_beam_power=1. #need to make sure that beam powers are match and deposition fraction matches too in both codes - LOCUST been assuming 1 so far
mom_transp['density']*=locust_beam_power/mom_transp['beam_source_captured']
mom_transp['NBI-heating-power(i1)']*=locust_beam_power/mom_transp['beam_source_captured']
mom_transp['NBI-heating-power(e-)']*=locust_beam_power/mom_transp['beam_source_captured']
mom_transp['beam_source']*=locust_beam_power/mom_transp['beam_source_captured']

time=3.1
time_index=(np.abs(mom_transp['time']-time)).argmin()

quantities=['density','NBI-heating-power(i1)','NBI-heating-power(e-)','beam_source']

fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)
axes=[ax1,ax2,ax3,ax4]

for quantity,ax in zip(quantities,axes):

    legend=['transp','locust']
    ax.plot(mom_transp['flux_pol_norm_sqrt'][time_index],mom_transp[quantity][time_index],'r-')
    ax.plot(mom_locust['flux_pol_norm_sqrt'],mom_locust[quantity],'g-')
    ax.set_title(quantity)

plt.legend(legend)
plt.show()