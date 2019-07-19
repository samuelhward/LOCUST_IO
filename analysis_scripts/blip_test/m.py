#ascot/FO/flat_prof_shifted/ascot_freia_1452884.h5
#ascot/GC/flat_prof_shifted_high_stats/ascot_freia_1700856.h5

import context
from classes.output_classes.moments import Moments as mom
import matplotlib.pyplot as plt
import numpy as np

m_a=mom('','ASCOT','ascot/FO/flat_prof_shifted/ascot_freia_1452884.h5')
m_l=mom('','LOCUST','locust/GC/flat_prof_shifted_high_stats/LOCUST_19-04-2019_12-17-32.730.h5')

fig,ax=plt.subplots(1)
#m_a.plot(key='NBI-heating-power(e-)',ax=ax,fig=fig,colmap='b')
m_l.plot(key='NBI-heating-power(e-)',ax=ax,fig=fig,colmap='g')
plt.show()

fig,ax=plt.subplots(1)
#m_a.plot(key='NBI-heating-power(i1)',ax=ax,fig=fig,colmap='b')
m_l.plot(key='NBI-heating-power(i1)',ax=ax,fig=fig,colmap='g')
plt.show()

fig,ax=plt.subplots(1)
m_a['density']/=-np.max(m_a['density'])/np.max(-m_a['NBI-heating-power(e-)'])
m_a.plot(key='density',ax=ax,fig=fig,colmap='b')
m_a.plot(key='NBI-heating-power(e-)',ax=ax,fig=fig,colmap='b')
m_l.plot(key='density',ax=ax,fig=fig,colmap='g')
plt.show()
