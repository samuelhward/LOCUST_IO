#XXX unfinished



import numpy as np
import context 
from classes.output_classes.moments import Moments as mom
from classes.input_classes.equilibrium import Equilibrium as equi
from classes.input_classes.wall import Wall as wall
from run_scripts.utils import  TRANSP_output as output_T
from run_scripts.utils import TRANSP_output_FI as dfn_t

from W03_files import * #grab all the filenames in this folder - assumes structure LOCUST_IO/input_files/ASCOT/29034W0?/*


mom_










dfn_LOCUST=dfn('',data_format='LOCUST',filename=file_LOCUST)
dfn_TRANSP=dfn_t('',filename=file_TRANSP)
#dfn_ASCOT=dfn('',data_format='ASCOT',filename=file_ASCOT)

eq=equi('',data_format='GEQDSK',filename=file_eq,GEQDSKFIX=1)
wall_MAST=wall('',data_format='UFILE',filename=file_wall)
output_TRANSP=output_T('',filename=file_TRANSP_output)

N_LOCUST=dfn_LOCUST.transform(axes=['N'])['dfn']
N_TRANSP=dfn_TRANSP.dfn_integrate()['dfn']
#N_ASCOT=dfn_ASCOT..transform(axes=['N'])['dfn']

BPCAP_index=np.abs(output_TRANSP['TIME']-time_slice).argmin()
BPCAP=output_TRANSP['BPCAP'][BPCAP_index]
dfn_TRANSP['dfn']/=N_TRANSP #normalise distribution functions against eachother
dfn_LOCUST['dfn']/=N_LOCUST

import matplotlib.pyplot as plt

fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)

axes=['E','V_pitch']
ax1.set_xlim(0,80000)
ax2.set_xlim(0,80000)
ax1.set_ylim(-1.1,1.1)
ax2.set_ylim(-1.1,1.1)
real_scale=False
limiters=False#wall_MAST
LCFS=False#eq
dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale)
dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax2,LCFS=LCFS,limiters=limiters,real_scale=real_scale)
#dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax3,LCFS=eq,limiters=wall_MAST,real_scale=real_scale)
axes=['E']
dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax4,colmap='g')
dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax4,colmap='r')
#dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax4,colmap='b')

plt.show()