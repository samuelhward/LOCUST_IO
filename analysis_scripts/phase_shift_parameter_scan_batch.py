#phase_shift_parameter_scan_batch.py
"""
Samuel Ward
07/11/2019
----
LOCUST_IO batch for scanning RMP phase shift
---
usage:
    see README.md for usage
 
notes:         
---
"""


import context
from classes.input_classes.perturbation import Perturbation as pert  
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import pathlib
import run_scripts.batch
import support

phase_min=0 #in degrees
phase_max=180
number_phases=10

phase_shifts=np.linspace(phase_min,phase_max,number_phases)
n2=[False]*number_phases #True #read, plot and dump n=2 harmonic?
n6=[True]*number_phases #read, plot and dump n=6 harmonic?
response=[False]*number_phases #include plasma response?
response_tag = ['response']*number_phases if response else ['vacuum']*number_phases  
ideal=[True]*number_phases #ideal or resistive?
ideal_tag = ['ideal']*number_phases if ideal else ['resistive']*number_phases  
data_format_input=['MARSF_bplas']*number_phases #data format of data source
data_format_output=['LOCUST']*number_phases #data format to dump data to

LOCUST_run__dir_LOCUST=[str(support.dir_locust)]*number_phases
LOCUST_run__system_name=[]*number_phases
LOCUST_run__repo_URL=[]*number_phases
LOCUST_run__commit_hash=[]*number_phases
LOCUST_run__settings_prec_mod={}
#set here
LOCUST_run__settings_prec_mod=[LOCUST_run__settings_prec_mod for number_simulations in range(number_phases)]
LOCUST_run__flags={}
#set here
LOCUST_run__settings_prec_mod=[LOCUST_run__settings_prec_mod for number_simulations in range(number_phases)]

RMP_batch_run=run_scripts.batch.Batch(
        phase_shifts=phase_shifts,
        n2=n2,
        n6=n6,
        response=response,
        response_tag=response_tag,
        ideal=ideal,
        ideal_tag=ideal_tag,
        data_format_input=data_format_input,
        data_format_output=data_format_output,
        LOCUST_run__dir_LOCUST=LOCUST_run__dir_LOCUST,
        LOCUST_run__system_name=LOCUST_run__system_name,
        LOCUST_run__repo_URL=LOCUST_run__repo_URL,
        LOCUST_run__commit_hash=LOCUST_run__commit_hash,
        LOCUST_run__settings_prec_mod=LOCUST_run__settings_prec_mod,
        LOCUST_run__flags=LOCUST_run__flags
        )

RMP_batch_run.launch(workflow_filepath='phase_shift_parameter_scan_run.py',system_name='TITAN')   

#################################

##################################################################

###################################################################################################