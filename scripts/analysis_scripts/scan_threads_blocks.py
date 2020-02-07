#scan_threads_blocks.py

'''
Samuel Ward
03/02/2020
----
plot script for analysing the block/thread/TIMAX/GPU generation timing runs
---
notes: 

---
'''


##################################################################
#Preamble

try:
    import sys
    import subprocess
    import pathlib
    import shlex
    import copy
    import numpy as np
    import ast
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
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    from classes.output_classes.rundata import Rundata
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/rundata.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

#grab information from the batch script including parameters, dirs etc.

sys.path.append("..")
try:
    import batch_scripts.profiling_batch
except:
    pass

#find all the output filepaths, times and errors across all runs for all GPUs
output_filepaths={}
output_times_total={}
output_times_total_error={}
output_times_total_mask={}
output_times_total_error_mask={}
for machine in ['TITAN']:#batch_scripts.profiling_batch.GPU_card_dispatch.keys(): 
    
    type_gpu=batch_scripts.profiling_batch.GPU_card_dispatch[machine]
    output_filepaths[type_gpu]=[] #initialise lists for holding filepaths, times and errors
    output_times_total[type_gpu]=[]
    output_times_total_error[type_gpu]=[]

    for dir_output in batch_scripts.profiling_batch.args_batch['LOCUST_run__dir_output']: #within each GPU folder the path to each output is the same
        dir_output=pathlib.Path(dir_output.strip("\'"))
        dir_output=dir_output.parents[1] / type_gpu / dir_output.parts[-1] #edit the full path to change GPU type 
        output_rundata_filepaths=list(dir_output.glob('rundata*')) #get all filenames for runs corresponding to this choice of parameters
        output_filepaths[type_gpu].append(output_rundata_filepaths)  

        if output_rundata_filepaths:
            _times=[]
            for output_filepath in output_rundata_filepaths:
                run_data=Rundata(ID='',data_format='LOCUST',filename=output_filepath)
                if run_data['time_total']: _times.append(run_data['time_total'])
            output_times_total[type_gpu].append(np.mean(_times))
            output_times_total_error[type_gpu].append(np.std(_times))
        else: 
            output_times_total[type_gpu].append(None)
            output_times_total_error[type_gpu].append(None)

    #populate masks for these arrays
    output_times_total_mask[type_gpu]=[True if _ is None or _ is np.nan else False for _ in output_times_total[type_gpu]] 
    output_times_total_error_mask[type_gpu]=[True if _ is None or _ is np.nan else False for _ in output_times_total_error[type_gpu]]
    #at this point we just have one flattened array - reshape into 3D array across threads, blocks and TIMAX dimensions
    output_times_total[type_gpu]=np.ma.masked_array(output_times_total[type_gpu],mask=output_times_total_mask[type_gpu]).reshape(len(batch_scripts.profiling_batch.parameter__threads),len(batch_scripts.profiling_batch.parameter__blocks),len(batch_scripts.profiling_batch.parameter__timax))
    output_times_total_error[type_gpu]=np.ma.masked_array(output_times_total_error[type_gpu],mask=output_times_total_error_mask[type_gpu]).reshape(len(batch_scripts.profiling_batch.parameter__threads),len(batch_scripts.profiling_batch.parameter__blocks),len(batch_scripts.profiling_batch.parameter__timax))

    #output_times_total[type_gpu]=np.ma.masked_array(output_times_total[type_gpu],mask=output_times_total_mask[type_gpu]).reshape(len(batch_scripts.profiling_batch.parameter__timax),len(batch_scripts.profiling_batch.parameter__blocks),len(batch_scripts.profiling_batch.parameter__threads))
    #output_times_total_error[type_gpu]=np.ma.masked_array(output_times_total_error[type_gpu],mask=output_times_total_error_mask[type_gpu]).reshape(len(batch_scripts.profiling_batch.parameter__timax),len(batch_scripts.profiling_batch.parameter__blocks),len(batch_scripts.profiling_batch.parameter__threads))


for counter,type_gpu in enumerate(output_filepaths.keys()):

    #plotting stuff
    import matplotlib.pyplot as plt
    fig,(ax1)=plt.subplots(1)
    legend=[]
    plot_styles=[settings.cmap_r,settings.cmap_b,settings.cmap_g]
    index_timax=3 #choose which slice we want to take
    index_threads=slice(None)
    index_blocks=slice(None)

    #go through and scale by number of markers in total
    for thread,num_threads in enumerate(batch_scripts.profiling_batch.parameter__threads):
        for block,num_blocks in enumerate(batch_scripts.profiling_batch.parameter__blocks):
            for timax,time_timax in enumerate(batch_scripts.profiling_batch.parameter__timax):
                if not output_times_total[type_gpu].mask[thread,block,timax]:
                    output_times_total[type_gpu][thread,block,timax]/=(num_threads*num_blocks)
                    print('{} {} {} {} {}'.format(
                    output_times_total[type_gpu][thread,block,timax],
                    num_threads,num_blocks,
                    batch_scripts.profiling_batch.parameter__timax[timax],
                    output_times_total[type_gpu].mask[thread,block,timax]))

    legend+=[type_gpu]
    Y,X=np.meshgrid(batch_scripts.profiling_batch.parameter__blocks,batch_scripts.profiling_batch.parameter__threads)
    Z=output_times_total[type_gpu][index_threads,index_blocks,index_timax]

    #ax1.set_facecolor(colmap(np.amin(dfn_copy[key])))
    #mesh=ax1.contour(X,Y,Z,levels=np.logspace(0,6,num=6),colors=plot_styles[counter](np.linspace(0,1,num=6)),edgecolor='none',linewidth=0,antialiased=True)
    #mesh=ax1.contour(X,Y,Z,levels=np.logspace(6,16,num=10,base=2),edgecolor='none',linewidth=0,antialiased=True)
    mesh=ax1.contourf(np.log2(X),np.log2(Y),np.log2(Z.astype(np.float64)),edgecolor='none',linewidth=0,antialiased=True)#,levels=np.linspace(-12,-3,10))
    cbar=fig.colorbar(mesh,ax=ax1,orientation='horizontal')
    cbar.set_label('log_2(time/marker) (s)')

    #if settings.plot_contour_labels:
    #    ax.clabel(mesh,inline=1,fontsize=10)

    #plt.errorbar(np.ma.masked_array(output_times_total, mask=output_times_total_mask),plot_styles[counter])
    ax1.set_xlabel('log_2(threads)')
    ax1.set_ylabel('log_2(blocks)')
    plt.show()

#################################

##################################################################

###################################################################################################