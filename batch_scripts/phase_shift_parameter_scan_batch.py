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
import settings

phase_min=0 #in degrees
phase_max=180
number_phases=5
phase_shift_degrees=np.linspace(phase_min,phase_max,number_phases)
phase_shift=phase_shift_degrees*2.*np.pi/360 #in radians

n2=[True]*number_phases #True #read, plot and dump n=2 harmonic?
n6=[False]*number_phases #read, plot and dump n=6 harmonic?
response=[True]*number_phases #include plasma response?
response_tag = ['response']*number_phases if response else ['vacuum']*number_phases  
ideal=[True]*number_phases #ideal or resistive?
ideal_tag = ['ideal']*number_phases if ideal else ['resistive']*number_phases  
data_format_input=['MARSF_bplas']*number_phases #data format of data source
data_format_output=['LOCUST']*number_phases #data format to dump data to

LOCUST_run__dirs_LOCUST=['locust_{}'.format(phase) for phase in phase_shift_degrees]
LOCUST_run__dir_LOCUST=[support.dir_locust / LOCUST_run__dir_LOCUST for LOCUST_run__dir_LOCUST in LOCUST_run__dirs_LOCUST]
LOCUST_run__system_name=[None]*number_phases
LOCUST_run__repo_URL=[None]*number_phases
LOCUST_run__commit_hash=[None]*number_phases

LOCUST_run__settings_prec_mod={}
LOCUST_run__settings_prec_mod['threadsPerBlock']=64
LOCUST_run__settings_prec_mod['blocksPerGrid']=256
LOCUST_run__settings_prec_mod['file_eqm']="'LOCUST_GEQDSK'"
LOCUST_run__settings_prec_mod['nmde']=1#2
LOCUST_run__settings_prec_mod['nnum']="[-2]"#"[-2,-6]"
LOCUST_run__settings_prec_mod['icoll']=0
LOCUST_run__settings_prec_mod['omega']='[0.0e0_gpu]'#'[0.0e0_gpu,0.0e0_gpu]'
LOCUST_run__settings_prec_mod['phase']='[0.0e0_gpu]'#'[0.0e0_gpu,0.0e0_gpu]'
LOCUST_run__settings_prec_mod['i3dr']=-1
LOCUST_run__settings_prec_mod=[LOCUST_run__settings_prec_mod]*number_phases
for phase,LOCUST_run__setting_prec_mod in zip(phase_shift_degrees,LOCUST_run__settings_prec_mod):
    LOCUST_run__setting_prec_mod['root']="'/tmp/locust_RMP_scan/{}'".format(phase)

LOCUST_run__flags={}
LOCUST_run__flags['LEIID']=8
LOCUST_run__flags['TOKAMAK']=2
LOCUST_run__flags['BRELAX']=True
LOCUST_run__flags['UNBOR']=100
LOCUST_run__flags['OPENMESH']=True
LOCUST_run__flags['OPENTRACK']=True
LOCUST_run__flags['PFCMOD']=True
LOCUST_run__flags['TOKHEAD']=True
LOCUST_run__flags['JXB2']=True
LOCUST_run__flags['PROV']=True
LOCUST_run__flags['PITCHCUR']=True
LOCUST_run__flags['EGSET']='100000.0D0'
LOCUST_run__flags['EBASE']=True
LOCUST_run__flags['UHST']=True
LOCUST_run__flags['LNLBT']=True
LOCUST_run__flags['GEQDSKFIX1']=True
LOCUST_run__flags['BP']=True
LOCUST_run__flags['TIMAX']='0.01D0'
LOCUST_run__flags['B3D']=True
LOCUST_run__flags['B3D_EX']=True
LOCUST_run__flags['SPLIT']=True
LOCUST_run__flags['NOINF']=True

LOCUST_run__flags=[LOCUST_run__flags]*number_phases

RMP_batch_run=run_scripts.batch.Batch(
        phase_shift=phase_shift,
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