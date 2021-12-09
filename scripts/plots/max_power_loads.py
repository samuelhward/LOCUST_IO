#script for plotting max power loads for d


"""
https://docs.paraview.org/en/latest/UsersGuide/selectingData.html
https://www.cfd-online.com/Forums/paraview/173602-get-maximum-minimum-values-field-variable-make-probes-them.html
select a view
hit "v"  or Edit > Find Data
select from Cells, and VTK file, and PFC flux and is max
hit run query
then you should be able to see the value...
"""

import numpy as np 
import context
import settings

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])


#n=3+6 90kAt
n_3_6=np.array([
[26, 0.184], 
[46, 0.183],
[66, 0.5327],
[86, 0.337],
[106, 0.322],
[126, 0.351],
]) 

#n=4+5 90kAt
n_4_5=np.array([
[26, 0.665],
[41, 0.428], #on lower panels between first wall and outer baffle
[56, 0.443], #on panels higher up
[71, 0.691],
[86, 0.466], #on panels higher up
[101, 0.498], #on panels higher up
])

import matplotlib.pyplot as plt 



fig,ax=plt.subplots(1)

ax.plot(n_3_6[:,0],n_3_6[:,1],linestyle='solid',label=r'$n=3+6$',color='r')
ax.plot(n_4_5[:,0],n_4_5[:,1],linestyle='solid',label=r'$n=4+5$',color=cmap_b(0.))
ax.axhline(np.mean(n_3_6,axis=0)[-1],linestyle='dashed',label=r'$n=3+6$ (mean)',color='r')
ax.axhline(np.mean(n_4_5,axis=0)[-1],linestyle='dashed',label=r'$n=4+5$ (mean)',color=cmap_b(0.))
ax.legend(fontsize=20)

ax.set_xlabel(r'Absolute phase shift of RMP ($\Phi_{\mathrm{m}}$) [deg]',fontsize=20)
ax.set_ylabel(r'Maximum first wall power load [$\mathrm{MWm}^{-2}$]',fontsize=20)

print(np.mean(n_3_6,axis=0)[-1]-np.mean(n_4_5,axis=0)[-1])

plt.show()

