#RMP_study_run.py

'''
Samuel Ward
14/10/2019
----
workflow designed for RMP studies (with static RMP field)
---
notes: 
    since RMP field is static, total field is combined before LOCUST simulation
    at first care was taken to decide as many parameters as possible outside of this script (best practice to be able to reuse this script) however I ran out of time and started deriving information herein
    all filepaths are absolute at the beginning
todo:
    XXX might need to call get_3D_fields twice, once for fundamental and another for harmonic - just make two identical functions but each referring to a different harmonic variable that is passed via parameters__toroidal_mode_numbers_fundamental, parameters__toroidal_mode_numbers_harmonics
    XXX or equally, call once but instead of parameters__toroidal_mode_numbers=[1] you could maybe have parameters__toroidal_mode_numbers=[[1,4],[2,6]] for two simulations with two harmonics each - this would also need to be appleid to RMP_study__filepaths_3D_field_U__batch

    XXX copy all BPLASMA original files to look the same, then re-run my workflow and see if output is only one coil row
    
    XXX for resolution scans, dXR and dXZ need to be varied

    XXX way to check all interpolated quantities is to overplot on the same axes for each flux_grid and they should all still overlap fine

    XXX add some progress messages to workflow stages
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
    import run_scripts.workflow
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/workflow.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/environment.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.LOCUST_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_run.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from run_scripts.MARS_builder_run import MARS_builder_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/MARS_builder_run.py could not be imported!\nreturning\n")
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
    from classes.input_classes.temperature import Temperature
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/temperature.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.number_density import Number_Density
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/number_density.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.rotation import Rotation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/rotation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.perturbation import Perturbation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
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

#need to make some additions to the LOCUST_run workflow component
class LOCUST_run_RMP(run_scripts.LOCUST_run.LOCUST_run):
    def get_inputs(self,*args,**kwargs):
        """
        LOCUST_run stage for moving inputs to correct location

        notes:
            edited for RMP workflow to move instead of copy
        """

        #copy input and cache files to correct location
        for file in self.dir_input.glob('*'): #move all input files to correct location
            subprocess.run(shlex.split('mv {file} {inputdir}'.format(file=str(file),inputdir=str(self.root / settings.username / self.tokhead / settings.LOCUST_dir_inputfiles_default))),shell=False)
        for file in self.dir_cache.glob('*'): #move all cache files to correct location
            subprocess.run(shlex.split('mv {file} {cachedir}'.format(file=str(file),cachedir=str(self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default))),shell=False)


class RMP_study_run(run_scripts.workflow.Workflow):
    """
    defines a workflow which runs a 3D field RMP run

    notes:
    usage:
    """ 

    def __init__(self,
                parameters__database,
                parameters__sheet_name_kinetic_prof,
                parameters__sheet_name_rotation,
                parameters__var_name_rotation,
                parameters__kinetic_prof_n, #variables defining the parameter set
                parameters__kinetic_prof_tF_tE,
                parameters__kinetic_prof_Pr,
                parameters__toroidal_mode_number,
                parameters__phase_upper,
                parameters__phase_middle,
                parameters__phase_lower,
                parameters__rotation_upper,
                parameters__rotation_middle,
                parameters__rotation_lower, 
                parameters__parameter_string,
                LOCUST_run__dir_LOCUST, #variables used by LOCUST_run (called within this workflow)
                LOCUST_run__dir_input,
                LOCUST_run__dir_output,
                LOCUST_run__dir_cache,
                LOCUST_run__system_name,
                LOCUST_run__repo_URL,
                LOCUST_run__commit_hash,
                LOCUST_run__settings_prec_mod,
                LOCUST_run__flags,
                NEMO_run__dir_NEMO,
                MARS_read__tail_U,
                MARS_read__tail_M,
                MARS_read__tail_L,
                MARS_read__settings,
                MARS_read__flags,
                RMP_study__name, #variables defining controlling this workflow
                RMP_study__dir_input_database,
                RMP_study__filepath_kinetic_profiles, 
                RMP_study__filepath_equilibrium,
                RMP_study__filepath_additional_data,
                RMP_study__filepaths_3D_field_U,
                RMP_study__filepaths_3D_field_M,
                RMP_study__filepaths_3D_field_L,
                IDS__shot,
                IDS__run,
                IDS__username,
                IDS__imasdb,
                *args,
                **kwargs):
        """
        notes:
        """


        super().__init__()


        #################################
        #read args needed for workflow steps

        self.parameters__database=parameters__database #variables defining the parameter set
        self.parameters__sheet_name_kinetic_prof=parameters__sheet_name_kinetic_prof
        self.parameters__sheet_name_rotation=parameters__sheet_name_rotation
        self.parameters__var_name_rotation=parameters__var_name_rotation
        self.parameters__kinetic_prof_n=parameters__kinetic_prof_n
        self.parameters__kinetic_prof_tF_tE=parameters__kinetic_prof_tF_tE
        self.parameters__kinetic_prof_Pr=parameters__kinetic_prof_Pr
        self.parameters__toroidal_mode_number=parameters__toroidal_mode_number
        self.parameters__phase_upper=parameters__phase_upper
        self.parameters__phase_middle=parameters__phase_middle
        self.parameters__phase_lower=parameters__phase_lower
        self.parameters__rotation_upper=parameters__rotation_upper
        self.parameters__rotation_middle=parameters__rotation_middle
        self.parameters__rotation_lower=parameters__rotation_lower
        self.parameters__parameter_string=parameters__parameter_string
        self.LOCUST_run__dir_LOCUST=pathlib.Path(LOCUST_run__dir_LOCUST) #data needed for LOCUST_run
        self.LOCUST_run__dir_input=pathlib.Path(LOCUST_run__dir_input)
        self.LOCUST_run__dir_output=pathlib.Path(LOCUST_run__dir_output)
        self.LOCUST_run__dir_cache=pathlib.Path(LOCUST_run__dir_cache)
        self.LOCUST_run__system_name=LOCUST_run__system_name
        self.LOCUST_run__repo_URL=LOCUST_run__repo_URL
        self.LOCUST_run__commit_hash=LOCUST_run__commit_hash
        self.LOCUST_run__settings_prec_mod=LOCUST_run__settings_prec_mod
        self.LOCUST_run__flags=LOCUST_run__flags
        self.NEMO_run__dir_NEMO=NEMO_run__dir_NEMO
        self.MARS_read__tail_U=MARS_read__tail_U
        self.MARS_read__tail_M=MARS_read__tail_M
        self.MARS_read__tail_L=MARS_read__tail_L
        self.MARS_read__settings=MARS_read__settings
        self.MARS_read__flags=MARS_read__flags
        self.RMP_study__name=RMP_study__name #data needed for additional steps
        self.RMP_study__dir_input_database=pathlib.Path(RMP_study__dir_input_database)
        self.RMP_study__filepath_kinetic_profiles=pathlib.Path(RMP_study__filepath_kinetic_profiles)
        self.RMP_study__filepath_equilibrium=pathlib.Path(RMP_study__filepath_equilibrium)
        self.RMP_study__filepath_additional_data=pathlib.Path(RMP_study__filepath_additional_data)
        self.RMP_study__filepaths_3D_field_U=pathlib.Path(RMP_study__filepaths_3D_field_U)
        self.RMP_study__filepaths_3D_field_M=pathlib.Path(RMP_study__filepaths_3D_field_M)
        self.RMP_study__filepaths_3D_field_L=pathlib.Path(RMP_study__filepaths_3D_field_L)
        self.IDS__shot=IDS__shot
        self.IDS__run=IDS__run
        self.IDS__username=IDS__username
        self.IDS__imasdb=IDS__imasdb

        #################################
        #derive some information

        #################################
        #add workflow stages

        self.add_command(command_name='mkdir',command_function=self.setup_RMP_study_dirs,position=1) #add new workflow stages
        self.add_command(command_name='get_kinetic',command_function=self.get_kinetic_profiles,position=2) #add new workflow stages
        self.add_command(command_name='get_3D',command_function=self.get_3D_fields,position=3) #add new workflow stages
        self.add_command(command_name='get_others',command_function=self.get_other_input_files,position=4) #add new workflow stages
        self.add_command(command_name='create_IDS',command_function=self.create_new_IDS,position=5) #add new workflow stages
        self.add_command(command_name='populate_IDS',command_function=self.inputs_to_IMAS,position=6) #add new workflow stages
        #self.add_command(command_name='run_NEMO',command_function=self.NEMO_run,position=7) #add new workflow stages
        #self.add_command(command_name='get_beam_deposition',command_function=self.get_beam_deposition,position=8) #add new workflow stages

    def setup_RMP_study_dirs(self,*args,**kwargs):
        """
        notes:
        """

        for direc in [
                        self.LOCUST_run__dir_input,
                        self.LOCUST_run__dir_output,
                        self.LOCUST_run__dir_cache,
                        self.LOCUST_run__dir_LOCUST
                        ]:
            if not direc.is_dir(): direc.mkdir(parents=True)

    def get_kinetic_profiles(self,*args,**kwargs):
        """
        prepare kinetic profiles for LOCUST
        notes:
            extracts from excel spreadsheet provided by Yueqiang
        """

        temperatures=[]
        densities=[]

        for species in ['electrons','deuterium']: #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804,'tritium','helium','hydrogen','tungsten','helium3']:
            try:
                densities.append(Number_Density(ID='',data_format='EXCEL1',species=species,filename=self.RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=self.parameters__sheet_name_kinetic_prof))
            except:
                pass
        for species in ['electrons','ions']:
            temperatures.append(Temperature(ID='',data_format='EXCEL1',species=species,filename=self.RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=self.parameters__sheet_name_kinetic_prof))
        for temperature in temperatures: temperature.dump_data(data_format='LOCUST',filename=self.LOCUST_run__dir_input / 'temperature_{}'.format(temperature.properties['species']))
        for density in densities: density.dump_data(data_format='LOCUST',filename=self.LOCUST_run__dir_input / 'density_{}'.format(density.properties['species']))
        rotation=Rotation(ID='',data_format='EXCEL1',filename=self.RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=self.parameters__sheet_name_kinetic_prof,rotation_name=self.parameters__var_name_rotation,sheet_name_rotation=self.parameters__sheet_name_rotation)
        rotation.dump_data(data_format='LOCUST',filename=self.LOCUST_run__dir_input / 'rotation') 

    def get_3D_fields(self,*args,**kwargs):
        """
        prepare 3D fields for LOCUST

        notes:
            copies 3D field files from the target folder into destination and renames file correctly before running mars_builder to combine them into single field
            suited for parameter scans involving static RMP fields
        """

        for file,TAIL,section_to_remove in zip(
                        [self.RMP_study__filepaths_3D_field_U,self.RMP_study__filepaths_3D_field_M,self.RMP_study__filepaths_3D_field_L],
                        [self.MARS_read__tail_U,self.MARS_read__tail_M,self.MARS_read__tail_L],
                        ['_cU_','_cM_','_cL_']
                        ): #move all input files to correct location and add appropriate TAIL

            destination=(str(self.LOCUST_run__dir_input / file.parts[-1])+TAIL).replace(section_to_remove,'_')
            HEAD=pathlib.Path((str(self.LOCUST_run__dir_input / file.parts[-1])).replace(section_to_remove,'_')) #while we are here grab the original filepath without added TAIL
            subprocess.run(shlex.split('cp {file} {inputdir}'.format(file=str(file),inputdir=destination)),shell=False)

        field_builder=MARS_builder_run(filepath_input=HEAD,dir_output=self.LOCUST_run__dir_input,dir_MARS_builder=self.LOCUST_run__dir_cache / 'MARS_builder',system_name='TITAN',settings_mars_read=self.MARS_read__settings,flags=self.MARS_read__flags)
        field_builder.run()

        for file in self.LOCUST_run__dir_input.glob(str(HEAD.parts[-1])+'*'):
            subprocess.run(shlex.split('rm {file}'.format(file=str(file))),shell=False) #delete old perturbation input files

    def get_other_input_files(self,*args,**kwargs):
        """
        fetch 2D equilibrium, tokamak wall description and cross section data
        notes:
        todo:
            make sure that beam depositions that are eventually generated are also copied over too
        """

        #equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=self.RMP_study__filepath_equilibrium)
        subprocess.run(shlex.split('cp {file} {inputdir}'.format(file=str(self.RMP_study__filepath_equilibrium),inputdir=str(self.LOCUST_run__dir_input / 'LOCUST_GEQDSK'))),shell=False)
        for file in self.RMP_study__filepath_additional_data.glob('*'):
            subprocess.run(shlex.split('cp {file} {inputdir}'.format(file=str(file),inputdir=str(self.LOCUST_run__dir_input / file.parts[-1]))),shell=False)

    def create_new_IDS(self,*args,**kwargs):
        """
        notes:
        """

        try:
            import imas 
        except:
            raise ImportError("ERROR: LOCUST_run_RMP.create_new_IDS could not import IMAS module!\nreturning\n")
            return

        NBI_IDS=imas.ids(130011,1) #take NBI geometry from sample public IDS
        NBI_IDS.open_env('public','ITER','3')
        NBI_IDS.nbi.get()

        new_IDS=imas.ids(self.IDS__shot,self.IDS__run) #initialise new blank IDS
        new_IDS.create_env(settings.username,settings.imasdb,'3')
        new_IDS_ctx=new_IDS.nbi.getPulseCtx() #retain file handle for later
        new_IDS.nbi=copy.deepcopy(NBI_IDS.nbi) #copy and overwrite NBI portion of IDS
        new_IDS.nbi.setPulseCtx(new_IDS_ctx) #reset file handle
        new_IDS.nbi.put()
        new_IDS.close()
        NBI_IDS.close()

    def inputs_to_IMAS(self,*args,**kwargs):
        """
        write current input data to IMAS

        notes:
            dumps kinetic profiles, perturbation and equilibrium 
            allows NEMO to be ran to generate beam deposition
        """

        try:
            import imas 
        except:
            raise ImportError("ERROR: LOCUST_run_RMP.inputs_to_IMAS could not import IMAS module!\nreturning\n")
            return

        data_elements_table={} #A,Z
        data_elements_table['deuterium']=[2.,1.]
        data_elements_table['electrons']=[0.,0.]
        #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804 data_elements_table['tritium']=[3.,1.]
        #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804 data_elements_table['helium']=[4.,2.]
        #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804 data_elements_table['hydrogen']=[1.,1.]
        #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804 data_elements_table['tungsten']=[184.,74.]
        #XXX HACK WITH SINGLE SPECIES FOR NOW - WAITING FOR https://jira.iter.org/browse/IMAS-2804 data_elements_table['helium3']=[3.,2.]

    #equilibrium
    #try:
        for filename_equilibrium in self.LOCUST_run__dir_input.glob('*GEQDSK*'):
            equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=filename_equilibrium)
            equilibrium.B_calc()
            equilibrium.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run)
    #except:
        print("ERROR: inputs_to_IMAS could not dump equilibrium to IDS!")

    #densities
    #try:
        for filename_density in self.LOCUST_run__dir_input.glob('density*'):
            species=str(filename_density.parts[-1]).split('_')[-1]
            density=Number_Density(ID='',data_format='LOCUST',filename=filename_density,species=species,A=data_elements_table[species][0],Z=data_elements_table[species][1])
            density['n']*=1000 #XXX debug
            density['flux_pol']=density['flux_pol_norm']*equilibrium['sibry']
            density.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run,species=species,A=data_elements_table[species][0],Z=data_elements_table[species][1])
    #except:
        print("ERROR: inputs_to_IMAS could not dump densities to IDS!")

    #temperatures
    #try:
        for filename_temperature in self.LOCUST_run__dir_input.glob('temperature*'):
            species=str(filename_temperature.parts[-1]).split('_')[-1]
            if species=='ions':#dump same ion temperature profile for each ion species
                for species_ in data_elements_table.keys():
                    if species_!='electrons': #loop through all non-electronic ion species and set equal temperature
                        temperature=Temperature(ID='',data_format='LOCUST',filename=filename_temperature,species=species_,A=data_elements_table[species_][0],Z=data_elements_table[species_][1])
                        temperature['T']*=0.01 #XXX debug
                        temperature['flux_pol']=temperature['flux_pol_norm'] #XXX hack for now - eventually need way of calculating
                        temperature.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run,species=species_,A=data_elements_table[species_][0],Z=data_elements_table[species_][1])
            else:
                temperature=Temperature(ID='',data_format='LOCUST',filename=filename_temperature,species=species)
                temperature['T']*=0.01 #XXX debug
                temperature['flux_pol']=temperature['flux_pol_norm']*equilibrium['sibry']
                temperature.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run,species=species,A=data_elements_table[species][0],Z=data_elements_table[species][1])
    #except:
        print("ERROR: inputs_to_IMAS could not dump temperatures to IDS!")

    #rotation
    #try:
        for filename_rotation in self.LOCUST_run__dir_input.glob('rotation'):
            for species_ in data_elements_table.keys(): #dump same rotation profile for each ion species
                if species_!='electrons': #loop through all non-electronic ion species and set equal temperature
                    rotation=Rotation(ID='',data_format='LOCUST',filename=filename_rotation,species=species_,A=data_elements_table[species_][0],Z=data_elements_table[species_][1])
                    rotation['flux_pol']=temperature['flux_pol_norm']*equilibrium['sibry'] #XXX hack for now - eventually need way of calculating
                    rotation.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run,A=data_elements_table[species_][0],Z=data_elements_table[species_][1])
    #except:
        print("ERROR: inputs_to_IMAS could not dump rotation to IDS!")

    #perturbation
    #try:
        for filename_perturbation in self.LOCUST_run__dir_input.glob('*BPLASMA*'):
            perturbation=Perturbation(ID='',data_format='LOCUST',filename=filename_perturbation,mode_number=self.parameters__toroidal_mode_number)
            perturbation.dump_data(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run)        
    #except:
        print("ERROR: inputs_to_IMAS could not dump perturbations to IDS!")
    
        #dump additional data required by NEMO to IMAS here
        
        new_IDS=imas.ids(self.IDS__shot,self.IDS__run) #initialise new blank IDS
        new_IDS.open_env(settings.username,settings.imasdb,'3')
        new_IDS.equilibrium.get()
        new_IDS.core_profiles.get()

        new_IDS.equilibrium.time_slice[0].boundary.geometric_axis.r=6.2 #hard-coded major radius from excel file
        new_IDS.equilibrium.time_slice[0].boundary.minor_radius=2. #hard-coded minor radius from excel file
        _,b_field_tor,_=equilibrium.B_calc_point(R=[equilibrium['rmaxis'].item(0)],Z=[equilibrium['zmaxis'].item(0)]) #dump direction of magnetic field
        new_IDS.equilibrium.time_slice[0].global_quantities.magnetic_axis.b_field_tor=float(b_field_tor[0])

        rho_tor=np.sqrt(np.abs(equilibrium['flux_tor'])*2.*np.pi/(np.pi*np.abs(equilibrium['bcentr']))) #calculate rho_tor on equilibrium flux grid - rho_tor = sqrt(b_flux_tor/(pi*b0)) where I think b_flux_tor is in [Wb]
        new_IDS.equilibrium.time_slice[0].profiles_1d.rho_tor=rho_tor #dump rho_tor on this grid to equilibrium IDS
        interpolator_rho_tor=processing.utils.interpolate_1D(equilibrium['flux_pol'],rho_tor) #interpolate this onto flux grid of kinetic profiles
        new_IDS.core_profiles.profiles_1d[0].grid.rho_tor=interpolator_rho_tor(temperature['flux_pol']) #XXX REQUIRE THIS

        _,Vr=run_scripts.utils.read_kinetic_profile_data_excel_1(filepath=self.RMP_study__filepath_kinetic_profiles,x='Fp',y='Vr',sheet_name=self.parameters__sheet_name_kinetic_prof) #volume profile on kinetic profile flux grid
        hro=run_scripts.utils.read_kinetic_profile_data_excel_1(filepath=self.RMP_study__filepath_kinetic_profiles,y='hro',sheet_name=self.parameters__sheet_name_kinetic_prof)
        dV=Vr*hro 
        interpolator_volume=processing.utils.interpolate_1D(temperature['flux_pol'],dV) #now need to interpolate volume profile to equilibrium flux grid
        new_IDS.equilibrium.time_slice[0].profiles_1d.volume=interpolator_volume(equilibrium['flux_pol'])

        #dump (R, Z) grids of toroidal flux coordinate ρ but normalised, since it is renormalised in NEMO anyway (and excel data does not contain un-normalised ρ)
        #this is now all done when dumping GEQDSK    
        flux_tor_2D=processing.utils.flux_func_to_RZ(equilibrium['flux_pol'],equilibrium['flux_tor'],equilibrium) 
        flux_tor_2D=(equilibrium['flux_tor']-np.min(equilibrium['flux_tor']))/(np.max(equilibrium['flux_tor'])-np.min(equilibrium['flux_tor'])) #normalise according to the 1D flux profile since this extends only to LCFS
        #XXX at this point check whether flux_tor stored in GEQDSK is similar to that stored in excel spreadsheet 
        #XXX flux_pol_norm,flux_tor_norm_sqrt=run_scripts.utils.read_kinetic_profile_data_excel_1(filepath=RMP_study__filepath_kinetic_profiles,x='Fp',y='x',sheet_name=parameters__sheet_name_kinetic_prof)
        #XXX flux_tor_norm=flux_tor_norm_sqrt**2

        #(R, Z) coordinates of the tokamak wall and/or SOL radius [m]
        #XXX read wall from GEQDSK here and dump to NEMO-friendly wall here

        new_IDS.core_profiles.put()
        new_IDS.equilibrium.put()
        new_IDS.close()

#XXX works up to here

    def NEMO_run(self,*args,**kwargs):
        """
        notes:
            assumes you have already cloned and compiled the nemo source and actor into dir_NEMO 
            takes most code from nemo source run_nemo.py by Mireille Schneider
        """

        # --------------------------------------------
        # PYTHON WRAPPER TO CALL PHYSICS CODE
        # --------------------------------------------
        # STRATEGY: USE PYUAL ONLY TO CALL ACTORS,
        # FOR THE REST USE PYTHON IMAS API'S DIRECTLY
        # --------------------------------------------

        # NEEDED MODULES
        import os,imas,sys

        # IMPORT MODULE(S) FOR SPECIFIC PHYSICS CODE(S)
        actor_path = os.path.join(os.getenv('KEPLER'), 'imas/src/org/iter/imas/python')
        list_of_actors = ['nemo']
        for name in list_of_actors:
          sys.path[:0] = [os.path.join(actor_path,name)]
          globals()[name] = getattr(__import__(name), name)

        # LOCAL DATABASE ENVIRONMENT
        user = os.getenv('USER')
        tokamakname = 'iter'
        version = os.getenv('IMAS_VERSION')[0]

        # OPEN INPUT DATAFILE TO GET DATA FROM IMAS SCENARIO DATABASE
        input = imas.ids(self.IDS__shot,self.IDS__run,0,0)
        input.open_env(self.IDS__username,self.IDS__imasdb,'3')

        # READ INPUT IDS'S AND CLOSE INPUT FILE
        input.equilibrium.get()
        input.core_profiles.get()
        input.nbi.get()
        input.distribution_sources.get()
        input.distributions.get()
        input.close()

        # DURATION TO SIMULATE (S)
        dt_required = 2.0

        # IF DISTRIBUTIONS HAS NEVER BEEN FILLED
        # (IF THIS IS THE FIRST TIME STEP OF THE SIMULATION)
        # THEN HOMOGENEOUS_TIME HAS NEVER BEEN FILLED
        # => SYSTEMATICALLY FILLED HERE THEN
        input.distributions.ids_properties.homogeneous_time = 1

        # OPEN OUTPUT OBJECT, IN VIEW OF SAVING RESULTS TO LOCAL DB
        output = imas.ids(self.IDS__shot,self.IDS__run,0,0)

        # EXECUTE PHYSICS CODE
        fokker_flag = 0
        nmarker = 50000
        output.distributions = nemo(input.equilibrium,input.core_profiles,input.nbi,
                                    input.distribution_sources,fokker_flag,nmarker,'input/input_nemo_sa_imas.xml')

        # CREATE OUTPUT DATAFILE AND EXPORT RESULTS TO LOCAL DATABASE
        output.create_env(self.IDS__username,self.IDS__imasdb,version)
        output.distribution_sources.put()
        output.close()

    def get_beam_deposition(self,*args,**kwargs):
        """
        notes:
        """

        beam_deposition=Beam_Deposition(data_format='IDS',shot=self.IDS__shot,run=self.IDS__run)
        beam_deposition.dump_data(data_format='LOCUST',filename=self.LOCUST_run__dir_input / 'ptcles.dat') 


'''
    def setup_directories(self,*args,**kwargs):
        """
        notes:
        """

    def run_LOCUST(self,*args,**kwargs)
        import locust_run here

        this_run=LOCUST_run_RMP(system_name=args.system_name,repo_URL=args.repo_URL,commit_hash=args.commit_hash,dir_locust=args.dir_locust,dir_input=args.dir_input,dir_output=args.dir_output,dir_cache=args.dir_cache,settings_prec_mod=settings_prec_mod,flags=flags)
        
        #default run order is: mkdir-get_code-make-get_input-run_code-get_output-cleanup
        #add new step to retrieve cache here
        def retrieve_cache():
            something

        this_run.add_command(retrieve_cache)

        this_run.run()

'''

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run RMP_study from command line in Python!')
    
    parser.add_argument('--parameters__database',type=str,action='store',dest='parameters__database',help="",default=None)
    parser.add_argument('--parameters__sheet_name_kinetic_prof',type=str,action='store',dest='parameters__sheet_name_kinetic_prof',help="",default=None)
    parser.add_argument('--parameters__sheet_name_rotation',type=str,action='store',dest='parameters__sheet_name_rotation',help="",default=None)
    parser.add_argument('--parameters__var_name_rotation',type=str,action='store',dest='parameters__var_name_rotation',help="",default=None)
    parser.add_argument('--parameters__kinetic_prof_n',type=str,action='store',dest='parameters__kinetic_prof_n',help="",default=None)
    parser.add_argument('--parameters__kinetic_prof_tF_tE',type=str,action='store',dest='parameters__kinetic_prof_tF_tE',help="",default=None)
    parser.add_argument('--parameters__kinetic_prof_Pr',type=str,action='store',dest='parameters__kinetic_prof_Pr',help="",default=None)
    parser.add_argument('--parameters__toroidal_mode_number',type=int,action='store',dest='parameters__toroidal_mode_number',help="",default=None)
    parser.add_argument('--parameters__phase_upper',type=str,action='store',dest='parameters__phase_upper',help="",default=None)
    parser.add_argument('--parameters__phase_middle',type=str,action='store',dest='parameters__phase_middle',help="",default=None)
    parser.add_argument('--parameters__phase_lower',type=str,action='store',dest='parameters__phase_lower',help="",default=None)
    parser.add_argument('--parameters__rotation_upper',type=str,action='store',dest='parameters__rotation_upper',help="",default=None)
    parser.add_argument('--parameters__rotation_middle',type=str,action='store',dest='parameters__rotation_middle',help="",default=None)
    parser.add_argument('--parameters__rotation_lower',type=str,action='store',dest='parameters__rotation_lower',help="",default=None)
    parser.add_argument('--parameters__parameter_string',type=str,action='store',dest='parameters__parameter_string',help="",default=None)
    parser.add_argument('--LOCUST_run__dir_LOCUST',type=str,action='store',dest='LOCUST_run__dir_LOCUST',help="",default=support.dir_locust)
    parser.add_argument('--LOCUST_run__dir_input',type=str,action='store',dest='LOCUST_run__dir_input',help="",default=support.dir_input_files)
    parser.add_argument('--LOCUST_run__dir_output',type=str,action='store',dest='LOCUST_run__dir_output',help="",default=support.dir_output_files)
    parser.add_argument('--LOCUST_run__dir_cache',type=str,action='store',dest='LOCUST_run__dir_cache',help="",default=support.dir_cache_files)
    parser.add_argument('--LOCUST_run__system_name',type=str,action='store',dest='LOCUST_run__system_name',help="",default='TITAN')
    parser.add_argument('--LOCUST_run__repo_URL',type=str,action='store',dest='LOCUST_run__repo_URL',help="",default=settings.repo_URL_LOCUST)
    parser.add_argument('--LOCUST_run__commit_hash',type=str,action='store',dest='LOCUST_run__commit_hash',help="",default=None)
    parser.add_argument('--LOCUST_run__settings_prec_mod',nargs='+',type=str,action='store',dest='LOCUST_run__settings_prec_mod',help="",default={})
    parser.add_argument('--LOCUST_run__flags',nargs='+',type=str,action='store',dest='LOCUST_run__flags',help="",default={})
    parser.add_argument('--NEMO_run__dir_NEMO',type=str,action='store',dest='NEMO_run__dir_NEMO',help="",default=support.dir_nemo)
    parser.add_argument('--MARS_read__tail_U',type=str,action='store',dest='MARS_read__tail_U',help="",default=None)
    parser.add_argument('--MARS_read__tail_M',type=str,action='store',dest='MARS_read__tail_M',help="",default=None)
    parser.add_argument('--MARS_read__tail_L',type=str,action='store',dest='MARS_read__tail_L',help="",default=None)
    parser.add_argument('--MARS_read__settings',nargs='+',type=str,action='store',dest='MARS_read__settings',help="",default={})
    parser.add_argument('--MARS_read__flags',nargs='+',type=str,action='store',dest='MARS_read__flags',help="",default={})
    parser.add_argument('--RMP_study__name',type=str,action='store',dest='RMP_study__name',help="",default=None)
    parser.add_argument('--RMP_study__dir_input_database',type=str,action='store',dest='RMP_study__dir_input_database',help="",default=None)
    parser.add_argument('--RMP_study__filepath_kinetic_profiles',type=str,action='store',dest='RMP_study__filepath_kinetic_profiles',help="",default=None)
    parser.add_argument('--RMP_study__filepath_equilibrium',type=str,action='store',dest='RMP_study__filepath_equilibrium',help="",default=None)
    parser.add_argument('--RMP_study__filepath_additional_data',type=str,action='store',dest='RMP_study__filepath_additional_data',help="",default=None)
    parser.add_argument('--RMP_study__filepaths_3D_field_U',type=str,action='store',dest='RMP_study__filepaths_3D_field_U',help="",default=None)
    parser.add_argument('--RMP_study__filepaths_3D_field_M',type=str,action='store',dest='RMP_study__filepaths_3D_field_M',help="",default=None)
    parser.add_argument('--RMP_study__filepaths_3D_field_L',type=str,action='store',dest='RMP_study__filepaths_3D_field_L',help="",default=None)
    parser.add_argument('--IDS__shot',type=int,action='store',dest='IDS__shot',help="",default=1)
    parser.add_argument('--IDS__run',type=int,action='store',dest='IDS__run',help="",default=1)
    parser.add_argument('--IDS__username',type=str,action='store',dest='IDS__username',help="",default=None)
    parser.add_argument('--IDS__imasdb',type=str,action='store',dest='IDS__imasdb',help="",default=None)

    args=parser.parse_args()

    #provide some extra parsing steps to dict-like input arguments
    args.LOCUST_run__settings_prec_mod=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__settings_prec_mod)
    args.LOCUST_run__flags=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__flags)
    args.MARS_read__settings=run_scripts.utils.command_line_arg_parse_dict(args.MARS_read__settings)
    args.MARS_read__flags=run_scripts.utils.command_line_arg_parse_dict(args.MARS_read__flags)

    this_run=RMP_study_run(
    parameters__database=args.parameters__database,
    parameters__sheet_name_kinetic_prof=args.parameters__sheet_name_kinetic_prof,
    parameters__sheet_name_rotation=args.parameters__sheet_name_rotation,
    parameters__var_name_rotation=args.parameters__var_name_rotation,
    parameters__kinetic_prof_n=args.parameters__kinetic_prof_n,
    parameters__kinetic_prof_tF_tE=args.parameters__kinetic_prof_tF_tE,
    parameters__kinetic_prof_Pr=args.parameters__kinetic_prof_Pr,
    parameters__toroidal_mode_number=args.parameters__toroidal_mode_number,
    parameters__phase_upper=args.parameters__phase_upper,
    parameters__phase_middle=args.parameters__phase_middle,
    parameters__phase_lower=args.parameters__phase_lower,
    parameters__rotation_upper=args.parameters__rotation_upper,
    parameters__rotation_middle=args.parameters__rotation_middle,
    parameters__rotation_lower=args.parameters__rotation_lower,
    parameters__parameter_string=args.parameters__parameter_string,
    LOCUST_run__dir_LOCUST=args.LOCUST_run__dir_LOCUST,
    LOCUST_run__dir_input=args.LOCUST_run__dir_input,
    LOCUST_run__dir_output=args.LOCUST_run__dir_output,
    LOCUST_run__dir_cache=args.LOCUST_run__dir_cache,
    LOCUST_run__system_name=args.LOCUST_run__system_name,
    LOCUST_run__repo_URL=args.LOCUST_run__repo_URL,
    LOCUST_run__commit_hash=args.LOCUST_run__commit_hash,
    LOCUST_run__settings_prec_mod=args.LOCUST_run__settings_prec_mod,
    LOCUST_run__flags=args.LOCUST_run__flags,
    NEMO_run__dir_NEMO=args.NEMO_run__dir_NEMO,
    MARS_read__tail_U=args.MARS_read__tail_U,
    MARS_read__tail_M=args.MARS_read__tail_M,
    MARS_read__tail_L=args.MARS_read__tail_L,
    MARS_read__settings=args.MARS_read__settings,
    MARS_read__flags=args.MARS_read__flags,
    RMP_study__name=args.RMP_study__name,
    RMP_study__dir_input_database=args.RMP_study__dir_input_database,
    RMP_study__filepath_kinetic_profiles=args.RMP_study__filepath_kinetic_profiles,
    RMP_study__filepath_equilibrium=args.RMP_study__filepath_equilibrium,
    RMP_study__filepath_additional_data=args.RMP_study__filepath_additional_data,
    RMP_study__filepaths_3D_field_U=args.RMP_study__filepaths_3D_field_U,
    RMP_study__filepaths_3D_field_M=args.RMP_study__filepaths_3D_field_M,
    RMP_study__filepaths_3D_field_L=args.RMP_study__filepaths_3D_field_L,
    IDS__shot=args.IDS__shot,
    IDS__run=args.IDS__run,
    IDS__username=args.IDS__username,
    IDS__imasdb=args.IDS__imasdb
        )

    this_run.run()

#################################

##################################################################

###################################################################################################