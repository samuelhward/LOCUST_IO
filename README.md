# LOCUST_IO Package

This "package" is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. There are two sub-modules - input and output - which contain python objects, each representing LOCUST data.

It is designed such that individual file read/write methods can easily be swapped out if new methods/filetypes come and go - all that needs to be changed is the respective function (which would be the same amount of work as writing the function bespoke).


Got any burning questions? Want to feedback? Please email the author, Sam, at shw515@york.ac.uk. Alternatively, raise an issue on the CCFE github!


Table of Contents
-----------------

* [Requirements](#requirements)
* [Getting Started](#getting-started)
* [Usage](#usage)
* [Overview](#overview)
	* [Classes Subpackage](#classes-subpackage)
		* [Input Module](#input-module)
		* [Output Module](#output-module)
	  	* [Support Module](#support-module)
	* [Testing Subpackage](#testing-subpackage)
	* [Plotting Subpackage](#plotting-subpackage)
* [Further Ideas And Contributing](#further-ideas-and-contributing)
* [Appendix](#appendix)
	* [Global Data Definitions And Variable Names](#global-data-definitions-and-variable-names)
		* [Equilibrium](#equilibrium)



## Requirements

* IMAS module with imasdb environment variable set
* Numpy/Matplotlib for plotting
* Python ___ or higher
*






## Getting Started

* Set up a folder within LOCUST_IO e.g. LOCUST_IO/my_project (this is where you will run the package)
* Copy a context.py file from LOCUST_IO/testing to LOCUST_IO/my_project
* Copy the file which you want to manipulate to LOCUST_IO/input_files and LOCUST_IO/output_files
* Import context
* Import input_classes or output_classes
* You're good to go!






## Usage

Some basic usage is outlined below. This includes importing/exporting an equilibrium for LOCUST between two different file formats.

```python
import context
import input_classes

my_equilibrium=input_classes.Equilibrium('some_identification',data_format='GEQDSK',input_filename='some.eqdsk') #my_equilibrium now holds all the data in one object
```

```python
my_equilibrium=input_classes.Equilibrium('some_identification') #you can initialise an empty equilibrium this way! (must always specify an ID) fill later with the read_dat() method

my_equilibrium.read_data(data_format='GEQDSK',input_filename='some.eqdsk')

my_equilibrium.dump_data(output_data_format='IDS',shot=1,run=1)
```

* You can also set individual pieces of data with the .set() method:

```python
my_equilibrium.copy(some_other_equilibrium) to copy all data
my_equilibrium.copy(some_other_equilibrium,'nh','nw','some_other_arg') to copy specific fields
```

* LOCUST inputs/outputs can be copied using the .copy() method:

```python
my_equilibrium.set(nw=5,fpol=[1,2,3,4]) or myeq.set(**{'nh':100,'fpol':some_external_array}) #to set multiple values simultaneously
```

* LOCUST_IO will also plot input/output data using the methods in the plotting/ directory:









## Overview

The overall layout consists of the main subpackage "classes" which takes data from the "input_files" and "output_files" folders, as well as individual folders for projects which draw from these folders. Every project folder, such as "testing" and "examples" must contain a context.py to ensure the package can see its folders and files. 

### Classes Subpackage

#### Input Module

Import this module if you want to handle inputs for LOCUST. File formats currently supported are:

* Equilibrium - GEQDSK/IDS
* NBI - ASCII
* Number Density - ASCII
* Temperature - ASCII

#### Output Module

#### Support Module


### Testing Subpackage

### Plotting Subpackage




## Further Ideas And Contributing

Please feel free to raise pull requests with this repo, it's designed to be easily contributed to. For example:

* Want to add a new file format type for an input equilibrium? Just standardise your data according to the .data[] member dictionary (see Appendix), copypaste the generic read/write functions (which you may have already written!) and add to the Equilibrium.read_data()/Equilibrium.dump_data() methods as below:

```python
elif data_format=='Baby': #say I want to import files that are encoded in the well-known Baby format, which requires a input_lemon and a pirate_tag

            if not none_check(self.ID,self.LOCUST_input_type,'cannot read_data from Baby format - lemon and pirate required\n',input_lemon,pirate_tag) #include a safety check 

                self.data_format=data_format #now just add the appropriate member data 
                self.lemon=lemon
                self.pirate=pirate #add these arguments to the read_data() argument list too
                self.data=read_Baby(self.lemon,self.pirate) #call my external read_Baby function
```





## Appendix

### Global Data Definitions & Variable Names

Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, here are all the different variable names used by different codes and file formats for similar quantities. '-' means they are not directly written to/read from and may be derived in other ways.

#### Equilibrium:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[GEQDSK](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)/[Equilibrium IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/equilibrium.html)/LOCUST)

    0D data
        /nw/-/nEQ_R                 			                            #number of points in R (x or width)
        /nh/-/nEQ_Z                                                			#number of points in Z (y or height)
        /idum/-/IDUM          			                                    #number of spatial dimensions?
        /rdim/-/RDIM                                             			#size of the R dimension in m
        /zdim/-/ZDIM         		                                        #size of the Z dimension in m
        r0/rcentr/...vacuum_toroidal_field.r0/RCENTR                 		#reference value of R
        b0/bcentr/...vacuum_toroidal_field.b0/BCENTR     		            #vacuum toroidal magnetic field at rcentr
        /rleft/-/RLEFT     			                                        #R at left (inner) boundary
        /zmid/-/ZMID    			                                        #Z at middle of domain (from origin)
        /rmaxis/...global_quantities.magnetic_axis.r/Rmagh      			#R at magnetic axis (O-point)
        /zmaxis/...global_quantities.magnetic_axis.z/Zmagh      			#Z at magnetic axis (O-point)
        /simag/...global_quantities.psi_axis/PSI_magh      		            #poloidal flux psi at magnetic axis (Weber / rad)
        /sibry/...global_quantities.psi_boundary/SIBRY      			    #poloidal flux psi at plasma boundary (Weber / rad)
        /current/...global_quantities.ip/CURRENT      		                #plasma current [Amps]   
        /xdum/-/XDUM                                              			#dummy variable - just contains zero
        /nbbbs/-/nb                                         				#number of points in the plasma boundary
        /limitr/-/IDUM                                          			#number of points in the wall boundary
    1D data
        /fpol/...profiles_1d.f/RBphih                            			#poloidal current function on uniform flux grid (1D array of f(psi)=R*B_toroidal [meter-Tesla])
        /pres/...profiles_1d.pressure/PRES      			                #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        /ffprime/...profiles_1d.f_df_dpsi/FFPRIM                   			#workk1 (check papers)
        /pprime/...profiles_1d.dpressure_dpsi/PPRIM     		          	#workk1 (check papers)
        /qpsi/...profiles_1d.q/QPSI    			                            #q values on uniform flux grid
        /rlim/...boundary.outline.r/	                           			#r wall boundary
        /zlim/...boundary.outline.z/ 				                        #z wall boundary
        /rbbbs/...boundary.lcfs.r/Rp  			                           	#r plasma boundary
        /zbbbs/...boundary.lcfs.z/Zp   		                           		#z plasma boundary
    2D data
        /psirz/...profiles_2d[0].psi/psi_equil_h     	                    #array (nx,ny) of poloidal flux on rectangular grid points (array of arrays)   


# TODO

## General

* Integrate with JET SAL API for instant access to JET data (will need error handling for use of this module on systems without access to JET SAL)
* Use fabric/paramiko for remote host handling stuff https://dtucker.co.uk/hack/ssh-for-python-in-search-of-api-perfection.html (fab and spur both built on paramiko and simplifies it, although more options obviously with paramiko)
* Make an example project which uses argparse to input command line arguments and then use the rest of the module to do batch operations or something 
* Add plotting functionality
* In the read_data/dump_data functions, could take data_format and then just **kwargs and then have none_check look in those /**kwargs for what the user has supplied? would mean that users can contribute new file formats but would not need to edit the arguement list in the read_data/dump_data functions - they would only need to copypaste the chunk of if logic as outlined above. also then it wouldn't matter what order users supplied their arguements at runtime - as long as they supply all the ones that are needed! this is good for hand holding.

* Reorganise classes into individual equilibrium.py, another_input.py...files if input_files.py gets too long
* Warn if writing to a filetype which holds less data than class instance currently holds - i.e. data will go missing! e.g. class has a "colour" and wants to write to a GEQDSK file (which doesn't have a colour field)
* Add a feature to warn if pre-existing file exists when writing out (to stop unwanted overwriting)
* need to decide how to standardise data in dicts 

## Equilibrium Things

* Need to add functionality to calculate toroidal current density when reading in GEQDSK. equivalent IDS is #time_slice(itime)/profiles_2d(i1)/j_tor 
* Generate the R,Z coordinate arrays when reading in a GEQDSK (currently only generated when writing out to IDS)
* Need to pass grid name and description to equilibrium IDS write out function / DECIDE WHAT TO DO WITH THIS DATA AND HOW TO PASS IT
