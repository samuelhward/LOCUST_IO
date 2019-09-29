#LOCUST_run.py

'''
Samuel Ward
28/09/2019
----
tools to run the LOCUST code
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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import context
except:
    raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_environment.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_build.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_edit_var
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_edit_var.py could not be imported!\nreturning\n")
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

class LOCUST_run:
    """
    class to run LOCUST

    notes:
    usage:
    """ 

    def __init__(self,system_name='TITAN',repo_URL=settings.repo_URL_LOCUST,commit_hash=None,dir_LOCUST=support.dir_locust,dir_output=support.dir_output_files):
        """
        notes:
            most information stored in LOCUST_run.environment and LOCUST_run.build, most init args are to init these instances
        args:
            system_name - string identifier to choose from selection of environments stored as class attributes 
            build - associated LOCUST_build
            repo_URL - SSH address of remote target git repo 
            commit_hash - commit hash identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
            dir_LOCUST - directory to temporarily store source code
            dir_output - directory to write results to 
        """
        
        files_in_dir_LOCUST=dir_LOCUST.glob('*') #check for files in dir_LOCUST
        if [file_in_dir_LOCUST for file_in_dir_LOCUST in files_in_dir_LOCUST]:
            print("ERROR: new LOCUST_run found files in target dir_LOCUST - please specify empty dir!\nreturning\n")
            return
        files_in_dir_output=dir_output.glob('*') #check for files in dir_output
        if [file_in_dir_output for file_in_dir_output in files_in_dir_output]:
            print("ERROR: new LOCUST_run found files in target dir_output - please specify empty dir!\nreturning\n")
            return

        self.dir_LOCUST=dir_LOCUST
        self.dir_output=dir_output

        try:
            self.dir_LOCUST.mkdir() 
        except:
            print("WARNING: LOCUST_run unable to mkdir dir_LOCUST - already exists?")
        try:
            self.dir_output.mkdir()
        except:
            print("WARNING: LOCUST_run unable to mkdir dir_output - already exists?")
        self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name=system_name) #create an environment for running
        self.build=run_scripts.LOCUST_build.LOCUST_build(system_name=system_name,repo_URL=repo_URL,commit_hash=commit_hash)

    def run(self,settings_prec_mod,flags):
        """
        notes:
        args:
            settings_prec_mod - dict denoting variable names and values to set them to within prec_mod
            flags - dict denoting compile flags
        """

        #retrieve code
        self.build.clone(directory=self.dir_LOCUST)
        self.dir_LOCUST=self.dir_LOCUST / 'locust' #from now on all code is within second folder due to clone

        #source code edits
        run_scripts.LOCUST_edit_var.LOCUST_edit_var(filename_in=self.dir_LOCUST / 'prec_mod.f90',filename_out=self.dir_LOCUST / 'prec_mod.f90',root="\'{}\'".format(str(self.dir_output)))
        run_scripts.LOCUST_edit_var.LOCUST_edit_var(filename_in=self.dir_LOCUST / 'prec_mod.f90',filename_out=self.dir_LOCUST / 'prec_mod.f90',**settings_prec_mod)

        #compile
        self.build.add_flags(**flags) 
        self.build.make(directory=self.dir_LOCUST,clean=False)

        #run LOCUST
        command=' '.join([self.environment.create_command(),'; ./locust'])
        try:
            with subprocess.Popen(command,shell=True,cwd=str(self.dir_LOCUST)) as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT
                pass
        except subprocess.CalledProcessError as err:
            raise(err)

        #retain copy of prec_mod.f90
        subprocess.run(shlex.split('mv {prec_mod} {dir_output}'.format(prec_mod=str(self.dir_LOCUST/'prec_mod.f90'),dir_output=str(self.dir_output))),shell=False)
        #remove source code folder
        #subprocess.run(shlex.split('rm -rf {}'.format(str(self.dir_LOCUST.parents[0]))),shell=False) #delete folder up since added new child folder when cloning

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run LOCUST from command line in Python!')

    parser.add_argument('--sys_name',type=str,action='store',default='TITAN',dest='system_name',help="computer system currently running on e.g. \'TITAN\'",required=False)
    parser.add_argument('--repo_URL',type=str,action='store',default=settings.repo_URL_LOCUST,dest='repo_URL',help="",required=False)
    parser.add_argument('--commit',type=str,action='store',default=None,dest='commit_hash',help="",required=False)
    parser.add_argument('--dir_LOCUST',type=str,action='store',default=support.dir_locust,dest='dir_LOCUST',help="",required=False)
    parser.add_argument('--dir_out',type=str,action='store',default=support.dir_output_files,dest='dir_output',help="",required=False)
    parser.add_argument('--settings_prec_mod',type=str,action='store',default={},dest='settings_prec_mod',help="",required=False)
    parser.add_argument('--flags',type=str,action='store',default={},dest='flags',help="",required=False)
    
    args=parser.parse_args()
    this_run=LOCUST_run(system_name=args.system_name,repo_URL=args.repo_URL,commit_hash=args.commit_hash,dir_LOCUST=args.dir_LOCUST,dir_output=args.dir_output)

    #NEED WAY OF PARSING PREC SETTINGS AND FLAGS HERE

    flags={}
    flags['TOKAMAK']=8
    flags['BRELAX']=None
    flags['UNBOR']=100
    flags['LEIID']=8
    flags['OPENMESH']=None
    flags['OPENTRACK']=None
    flags['PFCMOD']=None
    flags['TOKHEAD']=None
    flags['JXB2']=None
    flags['PROV']=None
    flags['GEQDSKFIX1']=None
    flags['PITCHCUR']=None
    flags['EGSET']='100000.0D0'
    flags['EBASE']=None
    flags['UHST']=None
    flags['LNLBT']=None
    flags['B3D']=None
    flags['B3D_EX']=None
    flags['SPLIT']=None
    flags['NOINF']=None
    flags['SMALLEQ']=None
    flags['BP']=None
    flags['TIMAX']='0.023D0'
    
    settings_prec_mod={}
    settings_prec_mod['threadsPerBlock']=987234
    #this_run.run(args.settings_prec_mod,args.flags)
    this_run.run(settings_prec_mod,flags)

#################################

##################################################################

###################################################################################################