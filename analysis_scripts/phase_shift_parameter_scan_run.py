#phase_shift_parameter_scan_run.py
"""
Samuel Ward
07/11/2019
----
LOCUST_IO workflow for scanning RMP phase shift
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
import run_scripts.workflow
import run_scripts.LOCUST_run
import run_scripts.utils
import support
import settings
import sys
import subprocess
import shlex


class RMP_scan_workflow(run_scripts.workflow.Workflow):
    def __init__(self,
        phase_shift,
        n2,
        n6,
        response,
        response_tag,
        ideal,
        ideal_tag,
        data_format_input,
        data_format_output,
        LOCUST_run__dir_LOCUST,
        LOCUST_run__system_name,
        LOCUST_run__repo_URL,
        LOCUST_run__commit_hash,
        LOCUST_run__settings_prec_mod,
        LOCUST_run__flags,
        *args,
        **kwargs
        ):
        super().__init__()

        #attach settings to class
        self.n2=n2
        self.n6=n6
        self.phase_shift=phase_shift
        self.response=response
        self.response_tag=response_tag
        self.ideal=ideal
        self.ideal_tag=ideal_tag
        self.data_format_input=data_format_input
        self.data_format_output=data_format_output
        self.dir_input_files=support.dir_input_files / 'RMP_phase_scan' / ('{phase_shift}_{ideal}_{response}'.format(phase_shift=self.phase_shift,ideal=self.ideal_tag,response=self.response_tag))
        self.dir_output_files=support.dir_output_files / 'RMP_phase_scan' / ('{phase_shift}_{ideal}_{response}'.format(phase_shift=self.phase_shift,ideal=self.ideal_tag,response=self.response_tag))
        self.dir_cache_files=support.dir_cache_files / 'RMP_phase_scan'
        self.LOCUST_run__dir_LOCUST=pathlib.Path(LOCUST_run__dir_LOCUST) / ('locust_'.format(self.phase_shift))
        self.LOCUST_run__dir_input=self.dir_input_files
        self.LOCUST_run__dir_output=self.dir_output_files
        self.LOCUST_run__dir_cache=self.dir_cache_files
        self.LOCUST_run__system_name=LOCUST_run__system_name
        self.LOCUST_run__repo_URL=LOCUST_run__repo_URL
        self.LOCUST_run__commit_hash=LOCUST_run__commit_hash
        self.LOCUST_run__settings_prec_mod=LOCUST_run__settings_prec_mod
        self.LOCUST_run__flags=LOCUST_run__flags

        self.add_command(command_name='mkdir',command_function=self.create_directories,position=1) 
        self.add_command(command_name='get_input',command_function=self.get_input_files,position=2) 
        self.add_command(command_name='get_3D',command_function=self.generate_3D_fields,position=3) 
        self.add_command(command_name='run_code',command_function=self.run_LOCUST,position=4) 
        self.add_command(command_name='cleanup',command_function=self.cleanup,position=5) 

    def create_directories(self,*args,**kwargs):
        """
        notes:
        """

        for direc in [
                        self.dir_input_files,
                        self.dir_output_files,
                        self.dir_cache_files
                        ]:
            if not direc.is_dir(): direc.mkdir(parents=True)

    def get_input_files(self,*args,**kwargs):
        """
        notes:
            assumes already ran ASCOT_2_LOCUST at this point
        """

        for file in (support.dir_input_files/'LOCUST').glob('*'):
            try:
                subprocess.run(shlex.split('cp {file} {inputdir}'.format(file=str(file),inputdir=str(self.dir_input_files))),shell=False)
            except:
                pass
        subprocess.run(shlex.split('mv {file_old} {file_new}'.format(file_old=str(self.dir_input_files/'g033143.02730'),file_new=str(self.dir_input_files/'LOCUST_GEQDSK'))),shell=False)
        subprocess.run(shlex.split('mv {file_old} {file_new}'.format(file_old=str(self.dir_input_files/'profile_Ti1.dat'),file_new=str(self.dir_input_files/'profile_Ti.dat'))),shell=False)

    def generate_3D_fields(self,*args,**kwargs):

        filepaths=[]
        harmonics=[]
        if self.n2:
            filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_mars_data')
            harmonics.append('n2')
        if self.n6:
            filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_n_6')
            harmonics.append('n6')

        for filepath,harmonic in zip(filepaths,harmonics): #cycle through harmonics
            perturbation=pert(ID='harmonic - {harmonic}, phase_shift - {phase_shift}'.format(harmonic=harmonic,phase_shift=self.phase_shift),data_format=self.data_format_input,filename=filepath,response=self.response,ideal=self.ideal,phase=self.phase_shift,bcentr=1.75660107,rmaxis=1.70210874)
            perturbation.dump_data(data_format=self.data_format_output,filename=self.dir_input_files/'BPLASMA_{harmonic}.dat'.format(harmonic=harmonic))

    def run_LOCUST(self,*args,**kwargs):
        run_LOCUST=run_scripts.LOCUST_run.LOCUST_run(dir_LOCUST=self.LOCUST_run__dir_LOCUST,
            dir_input=self.LOCUST_run__dir_input,
            dir_output=self.LOCUST_run__dir_output,
            dir_cache=self.LOCUST_run__dir_cache,
            system_name=self.LOCUST_run__system_name,
            repo_URL=self.LOCUST_run__repo_URL,
            commit_hash=self.LOCUST_run__commit_hash,
            settings_prec_mod=self.LOCUST_run__settings_prec_mod,
            flags=self.LOCUST_run__flags)
        run_LOCUST.run()

    def cleanup(self,*args,**kwargs):
        subprocess.run(shlex.split('rm -r {directory}'.format(directory=str(self.dir_input_files))),shell=False)

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run RMP scan')
    
    parser.add_argument('--phase_shift',type=float,action='store',dest='phase_shift',help="")
    parser.add_argument('--n2',type=run_scripts.utils.string_2_bool,action='store',dest='n2',help="")
    parser.add_argument('--n6',type=run_scripts.utils.string_2_bool,action='store',dest='n6',help="")
    parser.add_argument('--response',type=run_scripts.utils.string_2_bool,action='store',dest='response',help="")
    parser.add_argument('--response_tag',type=str,action='store',dest='response_tag',help="")
    parser.add_argument('--ideal',type=run_scripts.utils.string_2_bool,action='store',dest='ideal',help="")
    parser.add_argument('--ideal_tag',type=str,action='store',dest='ideal_tag',help="")
    parser.add_argument('--data_format_input',type=str,action='store',dest='data_format_input',help="")
    parser.add_argument('--data_format_output',type=str,action='store',dest='data_format_output',help="")
    parser.add_argument('--LOCUST_run__dir_LOCUST',type=str,action='store',dest='LOCUST_run__dir_LOCUST',help="",default=support.dir_locust)
    parser.add_argument('--LOCUST_run__system_name',type=str,action='store',dest='LOCUST_run__system_name',help="",default='TITAN')
    parser.add_argument('--LOCUST_run__repo_URL',type=str,action='store',dest='LOCUST_run__repo_URL',help="",default=settings.repo_URL_LOCUST)
    parser.add_argument('--LOCUST_run__commit_hash',type=str,action='store',dest='LOCUST_run__commit_hash',help="",default=None)
    parser.add_argument('--LOCUST_run__settings_prec_mod',nargs='+',type=str,action='store',dest='LOCUST_run__settings_prec_mod',help="",default={})
    parser.add_argument('--LOCUST_run__flags',nargs='+',type=str,action='store',dest='LOCUST_run__flags',help="",default={})

    args=parser.parse_args()

    #provide some extra parsing steps to dict-like input arguments
    args.LOCUST_run__settings_prec_mod=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__settings_prec_mod)
    args.LOCUST_run__flags=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__flags)

    this_run=RMP_scan_workflow(
        phase_shift=args.phase_shift,
        n2=args.n2,
        n6=args.n6,
        response=args.response,
        response_tag=args.response_tag,
        ideal=args.ideal,
        ideal_tag=args.ideal_tag,
        data_format_input=args.data_format_input,
        data_format_output=args.data_format_output,
        LOCUST_run__dir_LOCUST=args.LOCUST_run__dir_LOCUST,
        LOCUST_run__system_name=args.LOCUST_run__system_name,
        LOCUST_run__repo_URL=args.LOCUST_run__repo_URL,
        LOCUST_run__commit_hash=args.LOCUST_run__commit_hash,
        LOCUST_run__settings_prec_mod=args.LOCUST_run__settings_prec_mod,
        LOCUST_run__flags=args.LOCUST_run__flags
        )
    this_run.generate_3D_fields()

#################################

##################################################################

###################################################################################################