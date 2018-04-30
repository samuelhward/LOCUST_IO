# LOCUST_IO Package

This "package" is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. There are two class sub-modules - input and output - which contain python objects, each representing the types of data that LOCUST handles.

The code aims to be...

* __simple__ - It is designed to abstract the user away from LOCUST - operations such as converting from one data format to another can be done in two lines. Examples on how to use the API are included below, along with an example project. Transparency is also key as the code is heavily commented.

* __extensible__ - Contributing to the project is simple: to add a new filetype, for example, simply copypaste the read/write functions and use the [data dictionary](#global-data-definitions-and-variable-names) to match up your variables. Adding your own plotting and processing routines means just adding them to the respective files in *processing/*. 

* __portable__ - LOCUST_IO has all the infrastructure needed to encapsulate your automated pre/post processing scripts and can be ran 'out of the box' using the *example_project/* with the sample input/output files. Integration with LOCUST simply means cloning the latest version into *LOCUST/* folder!


Got any burning questions? Want to feedback? Please raise an issue here on github!
















Table of Contents
-----------------

* [Requirements](#requirements)
* [Getting Started](#getting-started)
* [Usage](#usage)
* [Appendix](#appendix)
    * [Global Data Definitions And Variable Names](#global-data-definitions-and-variable-names)
    * [Processing Routines](#processing-routines) 
* [Licence](#licence)
















## Requirements

Tested on:

* IMAS ≥3.16.0
* Python ≥v3.5.1
    * copy, re, itertools, math and time modules in the standard library 
    * numpy ≥v1.14.2
    * matplotlib ≥v2.0.0
    * scipy module ≥v0.17.0


















## Getting Started

* Set up a folder within LOCUST_IO e.g. *LOCUST_IO/my_project* (this is where you will run the package)
* Copy a context.py file from *LOCUST_IO/testing* or *LOCUST_IO/example_project* to *LOCUST_IO/my_project* (this sets up your paths etc)
* Copy the files which you want to manipulate to *LOCUST_IO/input_files* and *LOCUST_IO/output_files* (former/latter for inputs/outputs to/from LOCUST)

* Import context into your python session
* Import some classes from classes.input_classes or classes.output_classes
* You're good to go!





















## Usage


As well as the included *example_project/*, some basic usage is outlined below:

```python
#start by setting up python paths and import classes which hold LOCUST inputs:
import context
from classes.input_classes.equilibrium import Equilibrium

my_equilibrium=Equilibrium('ID_tag_for_this_equilibrium',data_format='GEQDSK',input_filename='some.eqdsk') 
#my_equilibrium now holds all the data in one object
#take a quick look at the equilibrium and its data
my_equilibrium.look()                               


#to initialise empty equilibrium to fill later with the read_data() (must always specify an ID):
my_equilibrium=Equilibrium('some_identification') 
#read data from GEQDSK
my_equilibrium.read_data(data_format='GEQDSK',input_filename='some.eqdsk')
#dump equilibrium to IDS format 
my_equilibrium.dump_data(output_data_format='IDS',shot=1,run=1) 


#you can set individual pieces of data with the .set() method. This will overwrite the default data format, which is numpy array:
#set multiple values simultaneously
my_equilibrium.set(nw=5,fpol=[1,2,3,4])  
#equally
my_equilibrium.set(**some_dict)                     


#your input/output objects can also be copied using the .copy() method:
#copy all data from one object to another
my_equilibrium.copy(some_other_equilibrium)
#copy specific fields                                        
my_equilibrium.copy(some_other_equilibrium,'B_field','some_key','some_other_key')  


#to get a quick glimpse of what you're working with, LOCUST_IO can also plot input/output data: 
from processing import plot_input
plot_input.plot_equilibrium(my_equilibrium)                                         
#(you can also stack plots onto the same axis object with the ax arguement)


#to check what data two objects share, use .compare():
my_equilibrium.compare(another_equilibrium,verbose=True)               


#to check if your object contains enough information for running LOCUST, use .run_check():
my_equilibrium.run_check(verbose=True)                                            


#you can also calculate new pieces of data using methods in the processing folder:
from processing import process_input
#calculate and set the magnetic field
my_equilibrium.set(B_field=process_input.B_calc(myeq))        
```



















## Appendix


### Global Data Definitions And Variable Names


Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, I've included all the different variable names used by different codes and file formats for similar quantities. '-' means they are not directly written to/read from and may be derived in other ways. If there are duplicates, that may be because the same class may contain data that can be written out in different ways (e.g. the temperature IDS may hold electron or ion temperature - which will dictate how the data is read/written).


#### Equilibrium:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[GEQDSK](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)/[Equilibrium IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/equilibrium.html)/LOCUST)

    0D data
        nR_1D/nw/-/nEQ_R                                                #number of points in R (x or width)
        nZ_1D/nh/-/nEQ_Z                                                #number of points in Z (y or height)
        idum/idum/-/IDUM                                                #dummy variable
        rdim/rdim/-/RDIM                                                #size of the R dimension in m
        zdim/zdim/-/ZDIM                                                #size of the Z dimension in m
        rcentr/rcentr/...vacuum_toroidal_field.r0/RCENTR                #reference value of R at magnetic axis
        zcentr/-/-/-                                                    #reference value of Z at magnetic axis
        bcentr/bcentr/...vacuum_toroidal_field.b0/BCENTR                #vacuum toroidal magnetic field at rcentr
        rleft/rleft/-/RLEFT                                             #R at left (inner) boundary
        zmid/zmid/-/ZMID                                                #Z at middle of domain (from origin)
        rmaxis/rmaxis/...global_quantities.magnetic_axis.r/Rmagh        #R at magnetic axis (O-point)
        zmaxis/zmaxis/...global_quantities.magnetic_axis.z/Zmagh        #Z at magnetic axis (O-point)
        simag/simag/...global_quantities.psi_axis/PSI_magh              #poloidal flux psi at magnetic axis (Weber / rad)
        sibry/sibry/...global_quantities.psi_boundary/SIBRY             #poloidal flux psi at plasma boundary (Weber / rad)
        current/current/...global_quantities.ip/CURRENT                 #plasma current [Amps]   
        xdum/xdum/-/XDUM                                                #dummy variable - just contains int(zero)
        lcfs_n/nbbbs/-/nb                                               #number of points in the plasma boundary
        limitr/limitr/-/IDUM                                            #number of points in the wall boundary
    1D data
        fpol/fpol/...profiles_1d.f/RBphih                               #poloidal flux function on uniform flux grid (1D array of f(psi)=R*B_toroidal [meter-Tesla]) (negative for positive plasma current in GEQDSK)
        pres/pres/...profiles_1d.pressure/PRES                          #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        ffprime/ffprime/...profiles_1d.f_df_dpsi/FFPRIM                 #F*d(F)/d(psi) where psi is poloidal flux per radian and  F=diamagnetic function=R*B_Phi 
        pprime/pprime/...profiles_1d.dpressure_dpsi/PPRIM               #plasma pressure * d(plasma pressure)/d(psi) (check papers)
        qpsi/qpsi/...profiles_1d.q/QPSI                                 #q profile (function of flux_pol)
        rlim/rlim/...boundary.outline.r/                                #r coordinates` of wall boundary
        zlim/zlim/...boundary.outline.z/                                #z coordinates of wall boundary
        lcfs_r/rbbbs/...boundary.lcfs.r/Rp                              #r coordinates of plasma boundary
        lcfs_z/zbbbs/...boundary.lcfs.z/Zp                              #z coordinates of plasma boundary
        R_1D/-/...profiles_2d[0].grid.dim1/                             #R dimension (m)
        Z_1D/-/...profiles_2d[0].grid.dim2/                             #Z dimension (m)
        flux_pol/-/...profiles_1d.psi/                                  #poloidal flux from magnetic axis up to the plasma boundary (Weber / rad)
        flux_tor/-/...profiles_1d.phi/                                  #toroidal flux from magnetic axis up to the plasma boundary (Weber / rad)
    2D data
        psirz/psirz/...profiles_2d[0].psi/psi_equil_h                   #poloidal flux at coordinate [r,z] in (Weber / rad) 

LOCUST reads an equilibrium in GEQDSK format.

#### Beam Deposition:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[Distribution_Sources IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/distribution_sources.html)/LOCUST)

    1D data
        R/...source[0].markers[0].positions[0][p]/                      #R coordinates of particle p
        phi/...source[0].markers[0].positions[1][p]/                    #toroidal angle of particle p
        Z/...source[0].markers[0].positions[2][p]/                      #Z coordinates of particle p
        V_R/...source[0].markers[0].positions[3][p]/                    #R component of v of particle p
        V_tor/...source[0].markers[0].positions[4][p]/                  #toroidal component of v of particle p
        V_Z/...source[0].markers[0].positions[5][p]/                    #Z component of v of particle p

LOCUST reads a birth profile in ASCII format as (R | phi | Z | V_R | V_tor | V_Z) .

#### Temperature:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[Core_Profiles IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/core_profiles.html)/LOCUST)

    1D data
        flux_pol/...profiles_1d[0].grid.psi/                           #poloidal flux (Weber / rad)
        flux_pol_norm/-/-                                              #normalised poloidal flux
        T/...profiles_1d[0].ion[0].temperature/                        #ion temperature
        T/...profiles_1d[0].electrons.temperature/                     #electron temperature
        flux_tor_coord/...profiles_1d[0].grid.rho_tor/                 #toroidal flux coordinate
        flux_tor/-/                                                    #toroidal flux (Weber / rad)
        q/...profiles_1d[0].q/                                         #safety factor

LOCUST reads a temperature profile in ASCII format as (normalised poloidal flux | temperature) 

#### Number Density:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[Core_Profiles IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/core_profiles.html)/LOCUST)

    1D data
        flux_pol/...profiles_1d[0].grid.psi/                            #poloidal flux (Weber / rad)
        flux_pol_norm/-/-                                               #normalised poloidal flux
        n/...profiles_1d[0].ion[0].density/                             #ion number density
        n/...profiles_1d[0].electrons.density/                          #electron number density
        flux_tor_coord/-/...profiles_1d[0].grid.rho_tor/                #toroidal flux coordinate
        flux_tor/-/-                                                    #toroidal flux (Weber / rad)
        q/...profiles_1d[0].q/                                          #safety factor

LOCUST reads a density profile in ASCII format as (normalised poloidal flux | number density) 

#### Orbits:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/LOCUST)

    0D data
        number_particles/                                               #total number of particles
        number_timesteps/                                               #total number of timesteps
    3D data
        orbits[t,i,p]/                                                  #spatial coordinate i for particle p at time step t

LOCUST dumps orbits in ASCII format as (number_particles \n R | phi | Z \n number_timesteps)

#### Final Particle List:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/TRANSP/LOCUST)

    0D data
        n/                                                              #number of particles per GPU (blocks per grid * threads per block)
        ngpu/                                                           #number of GPUs (OMP threads)
        niter/                                                          #
        npt_/                                                           #number of particle info slots
        nphc/                                                           #
        ntri/                                                           #
        number_particles/-/-                                            #total number of particles = n*ngpu
    1D data
        R/R/                                                            #r coordinate of particle
        phi/-/                                                          #phi coordinate of particle
        Z/Z/                                                            #z coordinate of particle
        V_R/-/                                                          #v_r coordinate of particle
        V_tor/-/                                                        #v_tor coordinate of particle
        V_Z/-/                                                          #v_z coordinate of particle
        V_pitch/vpll_v/                                                 #v_parallel/v
        energy/E/                                                       #energy of particle
        t///                                                            #time coordinate of particle
        status_flag//                                                   #status value of particle at this time
        status_flags//                                                  #possible status flags and their associated values

LOCUST dumps final particle lists in ASCII format as (n \n ngpu \n niter \n npt_ \n nphc \n ntri \n array[particle,npt_,nphc] )

#### Distribution Function:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/LOCUST)

    0D data
        IDFTYP/                                                         #dfn structure ID (= IDFTYP LOCUST flag)
        EQBM_MD5/                                                       #checksum ID of the equilibrium
        nE/                                                             #number of points in energy dimension
        nV/                                                             #number of points in velocity dimension
        nR/nF                                                           #number of R points on Dfn grid
        nZ/nF                                                           #number of Z points on Dfn grid
        nL/                                                             #number of points in Vphi/V dimension
        nPP/                                                            #number of points in Pphi dimension
        nMU/                                                            #number of points in Mu dimension
        nPSI/nPSIF                                                      #number of poloidal flux surface contours
        nPOL/nPOLF                                                      #number of poloidal cells for volume calculation
        nP/nPF                                                          #poloidal gyro-phase cell boundaries
        dEh/                                                            #energy bin width
        dVh/                                                            #velocity bin width
        dPPh/                                                           #Pphi bin width
        dMUh/                                                           #Mu bin width
        cpu_time/                                                       #
        0_1/                                                            #zero indicates fast ions only
        nc/                                                             #number of slices that Dfn is written to file in
        nR_1D/nEQ_R                                                     #number of R points on field grid
        nZ_1D/nEQ_Z                                                     #number of Z points on field grid
        Ai_1                                                            #first value of Ab
        Zi_1                                                            #first value of Zb
        Vsclh                                                           #velocity grid upper limit
        Vnrm                                                            #normalising velocity
        icoll                                                           #collisions (1=on 0=off)
        iscat                                                           #scattering (1=on 0=off)
        idiff                                                           #energy diffusion (1=on 0=off)
        iloss                                                           #charge exchange events (1=on 0=off)
        iterm                                                           #terminate if ptcl. leaves plasma
        niter                                                           #number of iterations for isym=1 simulation
        integrator/LEIID                                                #integrator type
        npnt                                                            #points per gyration
        0_2                                                             #the maximum integrator step size
        dt0                                                             #
        threads_per_block/threadsPerBlock                               #gpu threads per block
        blocks_per_grid/blocksPerGrid                                   #gpu blocks per grid
    1D data
        E/EF                                                            #energy dimension of Dfn grid
        V/VF                                                            #velocity dimension of Dfn grid
        PP/                                                             #Pphi dimension of Dfn grid
        MU/                                                             #Mu dimension of Dfn grid
        Jh/                                                             #Jacobian
        Jh_s/                                                           #Jacobian error
        R_1D/R_equil_h                                                  #R dimension of field grid
        Z_1D/Z_equil_h                                                  #Z dimension of field grid
        npolh                                                           #
        R/R                                                             #R dimension of Dfn grid
        Z/Z                                                             #Z dimension of Dfn grid
        L/L                                                             #pitch dimension of Dfn grid
        PP/PP                                                           #Pphi dimension of Dfn grid
        V/VF                                                            #velocity dimension of Dfn grid
        P/PG                                                            #special dimension - simulation specific
        Ab                                                              #fast ion masses
        Zb                                                              #trace particle Z
        Pdep/E0                                                         #normalised injected power
        tau_s                                                           #zeroth order slowing-down time
        E0                                                              #energy (plasma frame)
        EC                                                              #zeroth order critical energy
        rho                                                             #
        siglg                                                           ##r.m.s. width of test src         
    2D data
        psirz/psi_equil_h*PSI_sclh                                      #poloidal flux field [r,z] grid
        dVOL                                                            #
    3D data        
        dfn/Fh                                                          #Dfn grid (for IDFTYP=3)
        dfn_s/Fh_s                                                      #Dfn grid error (for IDFTYP=3)
    4D data
        csb/                                                            #volume element data
    5D data
        dfn/Fh                                                          #Dfn grid (for IDFTYP!=3)
        dfn_s/Fh_s                                                      #Dfn grid error (for IDFTYP!=3)

LOCUST dumps distribution functions in unformatted binary format. Different run-time flag combinations will dictate whether the above fields are written to file. The Dfn grid is in markers m^-3 (one must integrate to get per bin for plotting).
    
### Processing Routines

LOCUST_IO contains a few simple physics routines to process data (please refer to source code for instructions): 


     process_input

        QTP_calc                                 calculates the missing quantity out of Q, toroidal or poloidal flux (given two)
        fpolrz_calc                              calculates the 1D flux function on the 2D computational grid
        B_calc                                   calculates the components of the axisymmetric magnetic field on the 2D computational grid
        transform_marker_velocities              transforms marker phase space velocities to LOCUST r,phi,z format
        interpolate_1D                           returns a 1D interpolator (RBF)
        interpolate_2D                           returns a 2D interpolator (RBF/RBS)

     process_output

        dfn_integrate                            integrate a LOCUST dfn from s^3/m^6 to /bin
        dfn_collapse                             collapse a LOCUST dfn to the specified dimensions by summing
        particle_list_compression                opens and processes >>GB LOCUST particle lists in memory-efficient way

     plot_input

        plot_number_density                      plots the number density
        plot_temperature                         plots the temperature
        plot_beam_deposition                     plot a histogram/scatter of the beam deposition profile
        plot_equilibrium                         plot the equilibrium
        plot_B_field_stream                      plot the 2D poloidal magnetic vector field
        plot_B_field_line                        plot a single field line in 3D

     plot_output

        plot_orbits                              plot particle orbits
        plot_final_particle_list                 plot a histogram/scatter of the final particle list
        plot_distribution_function               plot the final distribution function














## Licence:


    LOCUST_IO
    Copyright (C) 2018  Samuel Ward

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public Licence as published by
    the Free Software Foundation, either version 3 of the Licence, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public Licence for more details.

    You should have received a copy of the GNU General Public Licence
    along with this program. If not, see [http://www.gnu.org/licences/](http://www.gnu.org/licences/).










## Useful Stuff:

* [Fortran format guide](https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/format.html)
* [GEQDSK format](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)