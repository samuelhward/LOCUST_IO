#resolution_parameter_scan_batch.py
"""
Samuel Ward
18/11/2019
----
LOCUST_IO batch for scanning perturbation grid resolution
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
import copy

number_simulations=10
study_name='convergence_test'

phase_min=0 #in degrees
phase_max=0
phase_shift_degrees=np.zeros(number_simulations)
phase_shift=phase_shift_degrees*2.*np.pi/360 #in radians

n2=[True]*number_simulations #True #read, plot and dump n=2 harmonic?
n6=[False]*number_simulations #read, plot and dump n=6 harmonic?
response=True #include plasma response?
response_tag = ['response']*number_simulations if response else ['vacuum']*number_simulations  
response=[response]*number_simulations #include plasma response?
ideal=False #ideal or resistive?
ideal_tag = ['ideal']*number_simulations if ideal else ['resistive']*number_simulations  
ideal=[ideal]*number_simulations #ideal or resistive?
data_format_input=['MARSF_bplas']*number_simulations #data format of data source
data_format_output=['LOCUST']*number_simulations #data format to dump data to

perturbation_nR=np.linspace(40,800,number_simulations)
perturbation_nZ=np.linspace(60,1200,number_simulations)

#define file structure here, uses <parameter_name>_<value>_.....format
dir_input_files=[str(support.dir_input_files / study_name / ('nR_1D_{nR_1D}_nZ_1D_{nZ_1D}_response_{response}'.format(nR_1D=nR_1D,nZ_1D=nZ_1D,response=response_tag_))) for nR_1D,nZ_1D,response_tag_ in zip(perturbation_nR,perturbation_nZ,response)]
dir_output_files=[str(support.dir_output_files / study_name / ('nR_1D_{nR_1D}_nZ_1D_{nZ_1D}_response_{response}'.format(nR_1D=nR_1D,nZ_1D=nZ_1D,response=response_tag_))) for nR_1D,nZ_1D,response_tag_ in zip(perturbation_nR,perturbation_nZ,response)]

LOCUST_run__dirs_LOCUST=['locust_{}'.format(phase) for phase in phase_shift_degrees]
LOCUST_run__dir_LOCUST=[support.dir_locust / LOCUST_run__dir_LOCUST for LOCUST_run__dir_LOCUST in LOCUST_run__dirs_LOCUST]
LOCUST_run__system_name=[None]*number_simulations
LOCUST_run__repo_URL=[None]*number_simulations
LOCUST_run__commit_hash=[None]*number_simulations

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
LOCUST_run__settings_prec_mod['niter']=1
LOCUST_run__settings_prec_mod=[copy.deepcopy(LOCUST_run__settings_prec_mod) for _ in range(number_simulations)]
for counter,phase in enumerate(phase_shift_degrees):
    LOCUST_run__settings_prec_mod[counter]['root']="'/tmp/{study_name}/{phase}'".format(phase=phase,study_name=study_name)

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
LOCUST_run__flags['TIMAX']='0.003D0' #about the time taken to saturate prompt losses
LOCUST_run__flags['B3D']=True
LOCUST_run__flags['B3D_EX']=True
LOCUST_run__flags['SPLIT']=True
LOCUST_run__flags['NOINF']=True

LOCUST_run__flags=[LOCUST_run__flags]*number_simulations

resolution_batch_run=run_scripts.batch.Batch(
        phase_shift=phase_shift,
        perturbation_nR=perturbation_nR,
        perturbation_nZ=perturbation_nZ,
        n2=n2,
        n6=n6,
        response=response,
        response_tag=response_tag,
        ideal=ideal,
        ideal_tag=ideal_tag,
        data_format_input=data_format_input,
        data_format_output=data_format_output,
        dir_input_files=dir_input_files,
        dir_output_files=dir_output_files,
        LOCUST_run__dir_LOCUST=LOCUST_run__dir_LOCUST,
        LOCUST_run__system_name=LOCUST_run__system_name,
        LOCUST_run__repo_URL=LOCUST_run__repo_URL,
        LOCUST_run__commit_hash=LOCUST_run__commit_hash,
        LOCUST_run__settings_prec_mod=LOCUST_run__settings_prec_mod,
        LOCUST_run__flags=LOCUST_run__flags
        )

resolution_batch_run.launch(workflow_filepath='phase_shift_parameter_scan_run.py',system_name='TITAN')   

#################################

##################################################################

###################################################################################################