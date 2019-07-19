from scipy.io import netcdf as ncdf
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
import copy





#PLOT PROMPT LOSS MAPS FOR EACH RUN

import context
from classes.output_classes.particle_list import Final_Particle_List as fpl
from processing import plot_output as plou 

extensions=['U60','U55','U59','U57','U58']
fluxes=[1.01,1.2,1.3,1.4,1.5]
colors=['r','g','b','y','m']

file_numbers=[2,3] #only these times have fast ion loss data unfortunately
times=['3.010','3.015']

#seperate plots
for time,file_number in zip(times,file_numbers):	
	fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)
	axes=[ax1,ax2,ax3,ax4,ax5,ax6]
	for extension,ax,flux,color in zip(extensions,axes,fluxes,colors):
		loss_list=fpl(ID=extension+'   XBMBND='+str(flux)+'   time='+time+'s',data_format='TRANSP',filename='157418'+extension+'_fi_'+str(file_number)+'_gc.cdf')
		plou.plot_final_particle_list(loss_list,type='scatter',axes=['R','Z'],fig=fig,ax=ax,colmap=color,status_flags=['all_losses'],real_scale=True)
	plt.show()

for time,file_number in zip(times,file_numbers):	
	fig,ax1=plt.subplots(1,1)
	for extension,flux,color in zip(extensions,fluxes,colors):
		loss_list=fpl(ID=extension+'   XBMBND='+str(flux)+'   time='+time+'s',data_format='TRANSP',filename='157418'+extension+'_fi_'+str(file_number)+'_gc.cdf')
		plou.plot_final_particle_list(loss_list,type='scatter',axes=['R','Z'],fig=fig,ax=ax1,colmap=color,status_flags=['all_losses'],real_scale=True)
	
	legend=[]
	for extension,flux in zip(extensions,fluxes):
		legend+=[extension+'   XBMBND='+str(flux)+'   time='+time+'s']
	plt.legend(legend)
	ax1.set_title('loss maps for fixed limiter TRANSP runs')
	plt.show()

'''
#PLOT TOTAL LOSSES OVER TIME FOR VARIOUS RUNS

extensions=['U60','U55','U59','U57']
quantities=['TIME','BSORB','BSORBPR','BPLIM']
fluxes=[1.01,1.2,1.3,1.4]#flux_max = 1.01,1.2 1.3 1.4 respectively
data={}

for extension in extensions:
    filename='157418'+extension+'.CDF'
    file=ncdf.netcdf_file(filename,'r')
    for quantity in quantities:
        data[extension+'_'+quantity]=copy.deepcopy(file.variables[quantity].data)
    file.close()
    file.close()


plot_quantity='BPLIM'

for counter,extension in enumerate(extensions):
    plt.plot(data[extension+'_TIME'],data[extension+'_'+plot_quantity],color=cmap_default(counter/len(extensions)),label=plot_quantity+' - '+str(fluxes[counter]))

plt.legend()
plt.show()

'''