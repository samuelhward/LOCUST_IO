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

* IMAS ≥3.16.0 (with imasdb set)
* Python ≥v2.7
* numpy module ≥v1.14
* matplotlib module
* scipy module ≥v1.0
* copy, re, itertools, math and time modules in the standard library 


















## Getting Started

* Set up a folder within LOCUST_IO e.g. *LOCUST_IO/my_project* (this is where you will run the package)
* Copy a context.py file from *LOCUST_IO/testing* or *LOCUST_IO/example_project* to *LOCUST_IO/my_project* (this sets up your paths etc)
* Copy the files which you want to manipulate to *LOCUST_IO/input_files* and *LOCUST_IO/output_files* (former/latter for inputs/outputs to/from LOCUST)

* Import context into your python session
* Import input_classes or output_classes into your python session
* You're good to go!





















## Usage


As well as the included *example_project/*, some basic usage is outlined below:

```python
import context
import input_classes

my_equilibrium=input_classes.Equilibrium('ID_tag_for_this_equilibrium',data_format='GEQDSK',input_filename='some.eqdsk') #my_equilibrium now holds all the data in one object
my_equilibrium.look()                               #quick look at the equilibrium and its data
```

```python
my_equilibrium=input_classes.Equilibrium('some_identification') #you can initialise an empty equilibrium this way which you can then fill later with the read_data() method as below (must always specify an ID)

my_equilibrium.read_data(data_format='GEQDSK',input_filename='some.eqdsk')

my_equilibrium.dump_data(output_data_format='IDS',shot=1,run=1) #you've now written out your equilibrium to IDS format
```

* You can set individual pieces of data with the .set() method. This will overwrite the default data format, which is numpy array:

```python
my_equilibrium.set(nw=5,fpol=[1,2,3,4])             #to set multiple values simultaneously
my_equilibrium.set(**some_dict)                     #equally

```

* Your input/output objects can also be copied using the .copy() method:

```python
my_equilibrium.copy(some_other_equilibrium)                                        #copy all data from one object to another
my_equilibrium.copy(some_other_equilibrium,'B_field','some_key','some_other_key')  #copy specific fields
```

* To get a quick glimpse of what you're working with, LOCUST_IO can also plot input/output data using the methods in the *processing/* directory:

```python
import plot_input

plot_input.plot_equilibrium(my_equilibrium)                                         #default psirz data plot
```

* To check if two objects share the same data, use .compare()

```python
my_equilibrium.compare(another_equilibrium,verbose=True)                            #returns true if my_equilibrium contains all data in another_equilibrium and that data is the same, verbose prints a summary
```

* To check if your object contains enough information to input into a LOCUST run, use .run_check()

```python
my_equilibrium.run_check(verbose=True)                                              #returns true if data contains correct fields
```

* You can also calculate new pieces of data using methods in the processing folder

```python
import process_input
my_equilibrium.set(B_field=process_input.B_calc(myeq))        #calculate and set the magnetic field
```



















## Appendix


### Global Data Definitions And Variable Names


Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, I've included all the different variable names used by different codes and file formats for similar quantities. '-' means they are not directly written to/read from and may be derived in other ways. If there are duplicates, that may be because the same class may contain data that can be written out in different ways (e.g. the temperature IDS may hold electron or ion temperature - which will dictate how the data is read/written).


#### Equilibrium:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[GEQDSK](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)/[Equilibrium IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/equilibrium.html)/LOCUST)

    0D data
        nw/nw/-/nEQ_R                                                   #number of points in R (x or width)
        nh/nh/-/nEQ_Z                                                   #number of points in Z (y or height)
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
        nbbbs/nbbbs/-/nb                                                #number of points in the plasma boundary
        limitr/limitr/-/IDUM                                            #number of points in the wall boundary
    1D data
        fpol/fpol/...profiles_1d.f/RBphih                               #poloidal flux function on uniform flux grid (1D array of f(psi)=R*B_toroidal [meter-Tesla]) (negative for positive plasma current in GEQDSK)
        pres/pres/...profiles_1d.pressure/PRES                          #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        ffprime/ffprime/...profiles_1d.f_df_dpsi/FFPRIM                 #F*d(F)/d(psi) where psi is poloidal flux per radian and  F=diamagnetic function=R*B_Phi 
        pprime/pprime/...profiles_1d.dpressure_dpsi/PPRIM               #plasma pressure * d(plasma pressure)/d(psi) (check papers)
        qpsi/qpsi/...profiles_1d.q/QPSI                                 #q profile (function of flux_pol)
        rlim/rlim/...boundary.outline.r/                                #r coordinates of wall boundary
        zlim/zlim/...boundary.outline.z/                                #z coordinates of wall boundary
        rbbbs/rbbbs/...boundary.lcfs.r/Rp                               #r coordinates of plasma boundary
        zbbbs/zbbbs/...boundary.lcfs.z/Zp                               #z coordinates of plasma boundary
        R_1D/-/...profiles_2d[0].grid.dim1/                             #R dimension (m)
        Z_1D/-/...profiles_2d[0].grid.dim2/                             #Z dimension (m)
        flux_pol/-/...profiles_1d.psi/                                  #poloidal flux from magnetic axis up to the plasma boundary (Weber / rad)
        flux_tor/-/...profiles_1d.phi/                                  #toroidal flux from magnetic axis up to the plasma boundary (Weber / rad)
    2D data
        psirz[r,z]/psirz/...profiles_2d[0].psi/psi_equil_h              #poloidal flux at coordinate r,z in (Weber / rad) 

LOCUST reads an equilibrium in GEQDSK format.

#### Beam Deposition:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[Distribution_Sources IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/distribution_sources.html)/LOCUST)

    1D data
        r/...source[0].markers[0].positions[0][p]/                      #r coordinates of particle p
        phi/...source[0].markers[0].positions[1][p]/                    #phi coordinates of particle p
        z/...source[0].markers[0].positions[2][p]/                      #z coordinates of particle p
        v_r/...source[0].markers[0].positions[3][p]/                    #r component of v of particle p
        v_phi/...source[0].markers[0].positions[4][p]/                  #phi component of v of particle p
        v_z/...source[0].markers[0].positions[5][p]/                    #z component of v of particle p

LOCUST reads a birth profile in ASCII format as (r | phi | z | v_r | v_phi | v_z) .

#### Temperature:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[Core_Profiles IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/core_profiles.html)/LOCUST)

    1D data
        flux_pol/...profiles_1d[0].grid.psi/                           #poloidal flux (Weber / rad)
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

LOCUST dumps orbits in ASCII format as (number_particles \n r | phi | z \n number_timesteps)

#### Final Particle List:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/LOCUST)

    0D data
        n/                                                              #number of particles per GPU (blocks per grid * threads per block)
        ngpu/                                                           #number of GPUs (OMP threads)
        niter/                                                          #
        npt_/                                                           #number of particle info slots
        nphc/                                                           #
        ntri/                                                           #
        number_particles/                                               #total number of particles = n*ngpu
    1D data
        r/                                                              #r coordinate of particle
        phi/                                                            #phi coordinate of particle
        z/                                                              #z coordinate of particle
        v_r/                                                            #v_r coordinate of particle
        v_phi/                                                          #v_phi coordinate of particle
        v_z/                                                            #v_z coordinate of particle
        t/                                                              #time coordinate of particle
        status_flag/                                                    #status of particle at this time
        status_flags/                                                   #possible status flags and their associated values

LOCUST dumps final particle lists in ASCII format as (n \n ngpu \n niter \n npt_ \n nphc \n ntri \n array[particle,npt_,nphc] )

    
### Processing Routines

LOCUST_IO contains a few simple physics routines to process data (please refer to source code for instructions): 

```python
    #process_input

        QTP_calc(Q=None,T=None,P=None)          #calculates the missing quantity out of Q, toroidal or poloidal flux (given two)
        fpolrz_calc                             #calculates the 1D flux function on the 2D computational grid
        B_calc                                  #calculates the components of the axisymmetric magnetic field on the 2D computational grid
        transform_marker_velocities             #transforms marker phase space velocities to LOCUST r,phi,z format
        interpolate_1D                          #returns a 1D interpolator (RBF)
        interpolate_2D                          #returns a 2D interpolator (RBF/RBS)

    #process_output

    #plot_input

        plot_number_density                     #plots the number density
        plot_temperature                        #plots the temperature
        plot_beam_deposition                    #plot a histogram/scatter of the beam deposition profile
        plot_equilibrium                        #plot the equilibrium
        plot_B_field_stream                     #plot the 2D poloidal magnetic vector field
        plot_B_field_line                       #plot a single field line in 3D

    #plot_output

        plot_orbits                             #plot particle orbit trails
        plot_final_particle_list                #plot a histogram/scatter of the final particle list
```







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
