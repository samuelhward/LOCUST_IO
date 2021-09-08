#plot_marker_number_scan.py
 
"""
Samuel Ward
17/08/20
----
look at effects of changing number of markers
---
 
notes:         
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
    import matplotlib.pyplot as plt
    import os
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    cwd=pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    sys.path.append(str(cwd.parents[1]))
    import templates.plot_mod
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)
    
################################################################## 
#Main 

import scan_marker_number_launch as batch_data

Pinj=33.e6
outputs=list(templates.plot_mod.get_output_files(batch_data,'fpl'))

fig = plt.figure()
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(221)
ax3 = fig.add_subplot(212)

fig_paper = plt.figure(constrained_layout=True)
ax_paper_total = fig_paper.add_subplot(111)
ax_paper_total.spines['top'].set_color('none')
ax_paper_total.spines['bottom'].set_color('none')
ax_paper_total.spines['left'].set_color('none')
ax_paper_total.spines['right'].set_color('none')
ax_paper_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_paper_1 = fig_paper.add_subplot(121)
ax_paper_2 = fig_paper.add_subplot(122)

for counter,(output,col_val) in enumerate(zip(outputs,np.linspace(0,0.9,len(batch_data.args_batch['LOCUST_run__dir_output'])))):
    if output: 
        output.set(theta=processing.utils.angle_pol(R_major=0.639277243e1,R=output['R'],Z=output['Z'],Z_major=0.597384943))
        output['E']/=1000. #convert to keV
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        number_markers=len(output['weight'])
        PFC_power=1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron
        print(f'number_markers={number_markers},FG={output["FG"][0]},f={output["f"]}')
        output['weight']*=1.e6*output['f']*output['E']*constants.charge_e/1.e3 #plot loss power of markers in kW
        output.plot(fig=fig,ax=ax1,axes=['R'],fill=False,label=number_markers,colmap=settings.cmap_inferno,colmap_val=col_val,number_bins=200,weight=True)
        output.plot(fig=fig,ax=ax2,axes=['E'],fill=False,label=output.ID,colmap=settings.cmap_inferno,colmap_val=col_val,number_bins=200,weight=True)
        ax3.scatter(np.log2(number_markers),100*PFC_power/Pinj,color=settings.cmap_inferno(col_val),label=number_markers,s=20)
        ax_paper_1.scatter(np.log2(number_markers),100*PFC_power/Pinj,color=settings.cmap_inferno(col_val),label=number_markers,s=20)
        output.plot(fig=fig_paper,ax=ax_paper_2,axes=['time'],fill=False,label=output.ID,colmap=settings.cmap_inferno,colmap_val=col_val,number_bins=500,weight=True)

#ax3.legend()
ax1.set_xlabel('R [m]')
ax2.set_xlabel('Marker energy [keV]')
ax3.set_xlabel('$\mathrm{log}_{2}(n_{\mathrm{markers}})$ ')
ax1.set_ylabel('Loss power [kW]')
ax2.set_ylabel('Loss power [kW]')
ax3.set_ylabel('Loss power [%]')
for ax in [ax1,ax2,ax3,ax_paper_2]: ax.set_title('')
ax_paper_1.set_title('')
ax_paper_1.set_xlabel(r'$\mathrm{log}_{2}(n_{\mathrm{markers}})$ ')
ax_paper_2.set_xlabel('Time [s]')
ax_paper_1.set_ylabel('Loss power')
ax_paper_1.text(0.05, 0.05, '[%]', ha='left', va='center', rotation='horizontal',transform=ax_paper_1.transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.)) 
ax_paper_2.text(0.05, 0.05, '[kW]', ha='left', va='center', rotation='horizontal',transform=ax_paper_2.transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.))
#ax_paper_total.set_ylabel('Loss power',labelpad=100)

plt.show()

#################################
 
##################################################################
 
###################################################################################################