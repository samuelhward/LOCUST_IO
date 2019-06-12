#GENERATE_DUMMY_IDS.py

'''
Samuel Ward
07/06/2019
----
Run this file to generate an IDS filled with example LOCUST inputs for testing
---
notes:
    some parts are not fully generated and require external file
---
'''



##################################################################
#Preamble

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

import numpy as np
import imas

from classes.input_classes.beam_deposition import Beam_Deposition as BD
from classes.input_classes.equilibrium import Equilibrium as EQ
from classes.input_classes.perturbation import Perturbation as P
from classes.input_classes.wall import Wall as W

##################################################################
#Main

def generate_test_IDS(shot,run,use_core_profiles=True,
    use_distribution_sources=True,use_equilibrium=True,
    use_mhd_linear=False,use_wall=False):
    """
    generate test IDS filled with dummy LOCUST input data

    notes:
    """

    test_IDS=imas.ids(shot,run)
    test_IDS.open_env('wards2','locust_imas','3')

    if use_core_profiles:

        temperature_e=np.linspace(1000,1,100)
        temperature_i=np.linspace(5000,1,100)
        density=np.linspace(1e19,1e6,100)
        flux_pol=np.linspace(0,1,100)

        test_IDS.core_profiles.get()
        test_IDS.core_profiles.code.version='1'
        test_IDS.core_profiles.code.name="LOCUST_IO"
        test_IDS.core_profiles.ids_properties.homogeneous_time=1   
        test_IDS.core_profiles.ids_properties.comment='a comment'
        test_IDS.core_profiles.time=np.array([0.0])

        test_IDS.core_profiles.profiles_1d.resize(1)
        test_IDS.core_profiles.profiles_1d[0].time=0.0
        test_IDS.core_profiles.profiles_1d[0].ion.resize(1)
        test_IDS.core_profiles.profiles_1d[0].ion[0].temperature=temperature_i
        test_IDS.core_profiles.profiles_1d[0].ion[0].density=density
        test_IDS.core_profiles.profiles_1d[0].grid.psi=flux_pol
        test_IDS.core_profiles.profiles_1d[0].electrons.temperature=temperature_e
        test_IDS.core_profiles.profiles_1d[0].electrons.density=density

        test_IDS.core_profiles.put()

    if use_distribution_sources:

        beam=BD(ID='LOCUST_IO test beam deposition')
        beam.set(R=np.random.uniform(size=10))
        beam.set(phi=np.random.uniform(size=10))
        beam.set(Z=np.random.uniform(size=10))
        beam.set(V_R=np.random.uniform(size=10))
        beam.set(V_tor=np.random.uniform(size=10))
        beam.set(V_Z=np.random.uniform(size=10))
        beam.set(weight=np.random.uniform(size=10))
        beam.dump_data(data_format='IDS',shot=shot,run=run)

    if use_equilibrium:

        filename_equilibrium='sample_GEQDSK.file'
        equil=EQ(ID='LOCUST_IO test equilibrium',filename=filename_equilibrium,data_format='GEQDSK')
        equil.dump_data(data_format='IDS',shot=shot,run=run)

    if use_mhd_linear:
        filename_mhd_linear='BPLASMA_n3'
        pert=P(ID='LOCUST_IO test mhd_linear')

    if use_wall:
        filename_wall=''
        wall=W(ID='LOCUST_IO test wall')

    test_IDS.close()
    
if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='generate dummy LOCUST input data in IMAS')

    parser.add_argument('--shot',type=int,action='store',dest='shot',help="IMAS shot number",required=True)
    parser.add_argument('--run',type=int,action='store',dest='run',help="IMAS run number",required=True)

    parser.add_argument('--core_profiles',type=str,action='store',default=True,dest='use_core_profiles',help="toggle whether to dump test core_profiles to IMAS",required=False)
    parser.add_argument('--distribution_sources',type=str,action='store',default=True,dest='use_distribution_sources',help="toggle whether to dump test distribution_sources to IMAS",required=False)
    parser.add_argument('--equilibrium',type=str,action='store',default=True,dest='use_equilibrium',help="toggle whether to dump test equilibrium to IMAS",required=False)
    parser.add_argument('--mhd_linear',type=str,action='store',default=False,dest='use_mhd_linear',help="toggle whether to dump test mhd_linear to IMAS",required=False)
    parser.add_argument('--wall',type=str,action='store',default=False,dest='use_wall',help="toggle whether to dump test wall to IMAS",required=False)
    
    args=parser.parse_args()

    #convert 
    generate_test_IDS(shot=args.shot,run=args.run,
        use_core_profiles=args.use_core_profiles,
        use_distribution_sources=args.use_distribution_sources,
        use_equilibrium=args.use_equilibrium,
        use_mhd_linear=args.use_mhd_linear,
        use_wall=args.use_wall)

#################################

##################################################################

###################################################################################################