#resolution_parameter_scan_run.py
"""
Samuel Ward
07/11/2019
----
LOCUST_IO workflow for scanning scanning perturbation grid resolution
---
usage:
    see README.md for usage
 
notes:         
---
"""


import context
from classes.input_classes.perturbation import Perturbation as pert  
from classes.input_classes.equilibrium import Equilibrium as equi
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


class resolution_scan_workflow(run_scripts.workflow.Workflow):
    def __init__(self,
        phase_shift,
        perturbation_nR,
        perturbation_nZ,
        n2,
        n6,
        response,
        response_tag,
        ideal,
        ideal_tag,
        data_format_input,
        data_format_output,
        dir_input_files,
        dir_output_files,
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
        self.phase_shift=phase_shift
        self.perturbation_nR=perturbation_nR
        self.perturbation_nZ=perturbation_nZ        
        self.n2=n2
        self.n6=n6
        self.response=response
        self.response_tag=response_tag
        self.ideal=ideal
        self.ideal_tag=ideal_tag
        self.data_format_input=data_format_input
        self.data_format_output=data_format_output
        self.dir_input_files=pathlib.Path(dir_input_files)
        self.dir_output_files=pathlib.Path(dir_output_files)
        self.dir_cache_files=support.dir_cache_files / 'resolution_scan'
        self.LOCUST_run__dir_LOCUST=LOCUST_run__dir_LOCUST
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

        equilibrium=equi(ID='',data_format='GEQDSK',filename=self.dir_input_files/'LOCUST_GEQDSK') #open up the equilibrium to extract some quantities for reading perturbation
        filepaths=[]
        harmonics=[]
        mode_numbers=[]
        if self.n2:
            filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_mars_data')
            harmonics.append('n2')
            mode_numbers.append(2)
        if self.n6:
            filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_n_6')
            harmonics.append('n6')
            mode_numbers.append(6)

        for filepath,harmonic,mode_number in zip(filepaths,harmonics,mode_numbers): #cycle through harmonics
            perturbation=pert(ID='harmonic - {harmonic}, phase_shift - {phase_shift}'.format(harmonic=harmonic,phase_shift=self.phase_shift),
                data_format=self.data_format_input,filename=filepath,mode_number=mode_number,response=self.response,ideal=self.ideal,phase=self.phase_shift,
                bcentr=equilibrium['bcentr'],rmaxis=equilibrium['rmaxis'],nR_1D=self.perturbation_nR,nZ_1D=self.perturbation_nZ)
            perturbation.dump_data(data_format=self.data_format_output,filename=self.dir_input_files/'BPLASMA_{harmonic}'.format(harmonic=harmonic))

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
        subprocess.run(shlex.split('rm -r {directory}'.format(directory=self.LOCUST_run__settings_prec_mod['root'])),shell=False)

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run resolution scan')
    
    parser.add_argument('--phase_shift',type=float,action='store',dest='phase_shift',help="")
    parser.add_argument('--perturbation_nR',type=int,action='store',dest='perturbation_nR',help="")
    parser.add_argument('--perturbation_nZ',type=int,action='store',dest='perturbation_nZ',help="")
    parser.add_argument('--n2',type=run_scripts.utils.string_2_bool,action='store',dest='n2',help="")
    parser.add_argument('--n6',type=run_scripts.utils.string_2_bool,action='store',dest='n6',help="")
    parser.add_argument('--response',type=run_scripts.utils.string_2_bool,action='store',dest='response',help="")
    parser.add_argument('--response_tag',type=str,action='store',dest='response_tag',help="")
    parser.add_argument('--ideal',type=run_scripts.utils.string_2_bool,action='store',dest='ideal',help="")
    parser.add_argument('--ideal_tag',type=str,action='store',dest='ideal_tag',help="")
    parser.add_argument('--data_format_input',type=str,action='store',dest='data_format_input',help="")
    parser.add_argument('--data_format_output',type=str,action='store',dest='data_format_output',help="")
    parser.add_argument('--dir_input_files',type=str,action='store',dest='dir_input_files',help="")
    parser.add_argument('--dir_output_files',type=str,action='store',dest='dir_output_files',help="")
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

    this_run=resolution_scan_workflow(
        phase_shift=args.phase_shift,
        perturbation_nR=args.perturbation_nR,
        perturbation_nZ=args.perturbation_nZ,
        n2=args.n2,
        n6=args.n6,
        response=args.response,
        response_tag=args.response_tag,
        ideal=args.ideal,
        ideal_tag=args.ideal_tag,
        data_format_input=args.data_format_input,
        data_format_output=args.data_format_output,
        dir_input_files=args.dir_input_files,
        dir_output_files=args.dir_output_files,
        LOCUST_run__dir_LOCUST=args.LOCUST_run__dir_LOCUST,
        LOCUST_run__system_name=args.LOCUST_run__system_name,
        LOCUST_run__repo_URL=args.LOCUST_run__repo_URL,
        LOCUST_run__commit_hash=args.LOCUST_run__commit_hash,
        LOCUST_run__settings_prec_mod=args.LOCUST_run__settings_prec_mod,
        LOCUST_run__flags=args.LOCUST_run__flags
        )
    this_run.run()

#################################

##################################################################

###################################################################################################