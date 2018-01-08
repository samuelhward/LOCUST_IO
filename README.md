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


ID='eq_ID_here' #specify the object data here
input_filename="test.eqdsk"
data_format='GEQDSK'
my_equilibrium=input_classes.Equilibrium(ID,data_format,input_filename) #my_equilibrium now holds all the data in one object
```

```python
ID='eq_ID_here' #specify the object data here
input_filename="test.eqdsk"
data_format='GEQDSK'
my_equilibrium=input_classes.Equilibrium(ID) #you can initialise an empty equilibrium this way! (must always specify an ID) fill later with the read_dat() method

my_equilibrium.read_data(data_format,input_filename)

output_data_format='IDS'
shot=1
run=1
my_equilibrium.dump_data(output_data_format=output_data_format,shot=shot,run=run)
```


* LOCUST inputs/outputs can be copied using the .copy() method:

* You can also set individual pieces of data with the .set() method:

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

Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, here are all the different variable names used by different codes and file formats for similar quantities. '-' means they are not used/are empty/not written to.

#### Equilibrium:

([LOCUST_IO](https://github.com/armoured-moose/LOCUST_IO)/[GEQDSK](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)/[Equilibrium IDS](https://portal.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/imas-3.7.3/equilibrium.html)/LOCUST)

    0D data
        /nw/-/nEQ_R            			#number of points in R (x or width)
        /nh/-/nEQ_Z            			#number of points in Z (y or height)
        /idum/-/IDUM          			#number of spatial dimensions?
        /rdim//RDIM         			#size of the R dimension in m
        /zdim//ZDIM         			#size of the Z dimension in m
        r0/rcentr/.r0/RCENTR      		#reference value of R
        b0/bcentr/.b0/BCENTR     		#vacuum toroidal magnetic field at rcentr
        /rleft/    /RLEFT     			#R at left (inner) boundary
        /zmid/      /ZMID    			#Z at middle of domain (from origin)
        /rmaxis/  /Rmagh      			#R at magnetic axis (O-point)
        /zmaxis/  /Zmagh      			#Z at magnetic axis (O-point)
        /simag /  /PSI_magh      		#poloidal flux psi at magnetic axis (Weber / rad)
        /sibry /  /SIBRY      			#poloidal flux psi at plasma boundary (Weber / rad)
        /current/ /CURRENT      		#plasma current [Amps]   
        /xdum/    /XDUM      			#dummy variable - just contains zero
        /nbbbs/   /nb      				#plasma boundary
        /limitr/   /IDUM     			#wall boundary
    1D data
        /fpol/    /RBphih      			#poloidal current function on uniform flux grid (1D array of f(psi)=R*B_toroidal [meter-Tesla])
        /pres/    /PRES      			#plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        /ffprime/  /FFPRIM     			#workk1 (check papers)
        /pprime/   /PPRIM     			#workk1 (check papers)
        /qpsi/      /QPSI    			#q values on uniform flux grid
        /rlim/          /				#r wall boundary
        /zlim/         / 				#z wall boundary
        /rbbbs/       /Rp  				#r plasma boundary
        /zbbbs/      /Zp   				#z plasma boundary
    2D data
        /psirz/    /psi_equil_h     	#array (nx,ny) of poloidal flux on rectangular grid points (array of arrays)   


#TODO

* Integrate with JET SAL API for instant access to JET data (will need error handling for use of this module on systems without access to JET SAL)
* Use fabric/paramiko for remote host handling stuff https://dtucker.co.uk/hack/ssh-for-python-in-search-of-api-perfection.html (fab and spur both built on paramiko and simplifies it, although more options obviously with paramiko)
* Use argparse to input command line arguments etc. to the python code which calls from input_classes. e.g. pass file names etc. via command line and argparse to a python script which implements this module
* Currently, this package will remain in the given folder structure. I understand that this may limit flexibility somewhat if the user is looking for something ultra-light, but sticking to this folder format means that the environment stays controlled and limits the variation across user environments to make reading/using/debugging/contributing as easy as possible. This may change in the future (may be I will implement a light and heavy version...)
* Add plotting functionality