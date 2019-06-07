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

from classes.input_classes.temperature import Temperature as T
from classes.input_classes.number_density import Number_Density as ND
from classes.input_classes.beam_deposition import Beam_Deposition as BD
from classes.input_classes.equilibrium import Equilibrium as EQ
from classes.input_classes.perturbation import Perturbation as P
from classes.input_classes.wall import Wall as W

##################################################################
#Main

def generate_test_IDS(shot,run,use_temperature=True,use_number_density=True,
    use_beam_deposition=True,use_equilibrium=True,use_perturbation=False,use_wall=False):
    """
    generate test IDS filled with dummy LOCUST input data

    notes:
    """

    if use_temperature:
        DUMMY_profile=np.linspace(1,100,100)
        DUMMY_flux_pol=np.linspace(0,1,100)
        temp_i=T(ID='LOCUST_IO test ion temperature',species='ions')
        temp_e=T(ID='LOCUST_IO test electron temperature',species='electrons')
        temp_i.set(flux_pol=DUMMY_flux_pol)
        temp_i.set(T=DUMMY_profile)
        temp_e.set(flux_pol=DUMMY_flux_pol)
        temp_e.set(T=DUMMY_profile)
        temp_i.dump_data(data_format='IDS',shot=shot,run=run,species='ions')
        temp_e.dump_data(data_format='IDS',shot=shot,run=run,species='electrons')

    if use_number_density:
        DUMMY_profile=np.linspace(1,100,100)
        DUMMY_flux_pol=np.linspace(0,1,100)
        numd_i=ND(ID='LOCUST_IO test ion density',species='ions')
        numd_e=ND(ID='LOCUST_IO test electron density',species='electrons')
        numd_i.set(flux_pol=DUMMY_flux_pol)
        numd_i.set(n=DUMMY_profile)
        numd_e.set(flux_pol=DUMMY_flux_pol)
        numd_e.set(n=DUMMY_profile)
        numd_i.dump_data(data_format='IDS',shot=shot,run=run,species='ions')
        numd_e.dump_data(data_format='IDS',shot=shot,run=run,species='electrons')

    if use_beam_deposition:
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

    if use_perturbation:
        filename_perturbation='BPLASMA_n3'
        pert=P(ID='LOCUST_IO test perturbation')

    if use_wall:
        filename_wall=''
        wall=W(ID='LOCUST_IO test wall')


if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='generate dummy LOCUST input data in IMAS')

    parser.add_argument('--shot',type=int,action='store',dest='shot',help="IMAS shot number",required=True)
    parser.add_argument('--run',type=int,action='store',dest='run',help="IMAS run number",required=True)

    parser.add_argument('--temperature',type=str,action='store',default=True,dest='use_temperature',help="toggle whether to dump test temperature to IMAS",required=False)
    parser.add_argument('--number_density',type=str,action='store',default=True,dest='use_number_density',help="toggle whether to dump test number_density to IMAS",required=False)
    parser.add_argument('--beam_deposition',type=str,action='store',default=True,dest='use_beam_deposition',help="toggle whether to dump test beam_deposition to IMAS",required=False)
    parser.add_argument('--equilibrium',type=str,action='store',default=True,dest='use_equilibrium',help="toggle whether to dump test equilibrium to IMAS",required=False)
    parser.add_argument('--perturbation',type=str,action='store',default=False,dest='use_perturbation',help="toggle whether to dump test perturbation to IMAS",required=False)
    parser.add_argument('--wall',type=str,action='store',default=False,dest='use_wall',help="toggle whether to dump test wall to IMAS",required=False)
    
    args=parser.parse_args()

    #convert 
    generate_test_IDS(shot=args.shot,run=args.run,use_temperature=args.use_temperature,
        use_number_density=args.use_number_density,
        use_beam_deposition=args.use_beam_deposition,
        use_equilibrium=args.use_equilibrium,
        use_perturbation=args.use_perturbation,
        use_wall=args.use_wall)

#################################

##################################################################

###################################################################################################