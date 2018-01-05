# LOCUST_IO Package

This "package" is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. There are two sub-modules - input and output - which contain python objects, each representing LOCUST data.

It is designed such that individual file read/write methods can easily be swapped out if new methods/filetypes come and go - all that needs to be changed is the respective function (which would be the same amount of work as writing the function bespoke).




Table of Contents
-----------------

* [Requirements](#Requirements)
* [Basic Usage](#Basic-Usage)
* [Advanced Usage](#Advanced-Usage)
* [Overview](#Overview)
	* [Classes Subpackage](#Classes-Subpackage)
		* [Input Module](#Input-Module)
		* [Output Module](#Output-Module)
	  	* [Support Module](#Support-Module)
	* [Testing Subpackage](#Testing-Subpackage)
		* [Input Class Tests](#Input-Class-Tests)
		* [Output Class Tests](#Output-Class-Tests)
	* [Plotting Subpackage](#Plotting-Subpackage)
* [Further Ideas And Contributing](#Further-Ideas-And-Contributing)
* [Appendix](#Appendix)
	* [Global Data Definitions & Variable Names](#Global-Data-Definitions-&-Variable-Names)
		* [Equilibrium](#Equilibrium)
	* [GEQDSK File Format](#GEQDSK-File-Format)




## Requirements

* IMAS module with imasdb environment variable set
* Numpy/Matplotlib for plotting
* Python ___ or higher
*






## Basic Usage

* Set up a folder within LOCUST_IO e.g. LOCUST_IO/my_project
* Copy a context.py file from LOCUST_IO/testing to LOCUST_IO/my_project
* Copy the file which you want to manipulate to LOCUST_IO/input_files and LOCUST_IO/output_files


```python
ID='eq_ID_here' #specify the object data here
input_filename="test.eqdsk"
data_format='GEQDSK'
my_equilibrium=input_classes.Equilibrium(ID,input_filename,data_format) #my_equilibrium now holds all the data in one object
```

```python
ID='eq_ID_here' #specify the object data here
input_filename="test.eqdsk"
data_format='GEQDSK'
my_equilibrium=input_classes.Equilibrium(ID) #initialise a blank my_equilibrium

my_equilibrium.read_data(input_filename, data_format)
my_equilibrium.dump_data()
```

* Everything MUST have an ID









## Advanced Usage

* You can instantiate a blank LOCUST input/output by instantiating by not specifying data_format or input_filename (you must always specify an ID). Fill the input/output later by using the .read_data() method. 
* The .copy() method can be used to...
* The .set() method can be used to...
* You can use this package interactively/non-interactively anywhere on your system as long as you append your python path to find the classes/ folder. However, input/output files will still be written to the input/output folders within the package directory.











## Overview

The overall layout consists of the main subpackage "classes" which takes data from the "input_files" and "output_files" folders, as well as individual folders for projects which draw from these folders. Every project folder, such as "testing" and "examples" must contain a context.py to ensure the package can see its folders and files. 

### Classes Subpackage

#### Input Module

#### Output Module

#### Support Module


### Testing Subpackage

#### Input Class Tests

#### Output Class Tests








## Further Ideas And Contributing


* Integrate with JET SAL API for instant access to JET data (will need error handling for use of this module on systems without access to JET SAL)
* Use fabric/paramiko for remote host handling stuff https://dtucker.co.uk/hack/ssh-for-python-in-search-of-api-perfection.html (fab and spur both built on paramiko and simplifies it, although more options obviously with paramiko)
* Use argparse to input command line arguments etc. to the python code which calls from input_classes. e.g. pass file names etc. via command line and argparse to a python script which implements this module
* Currently, this package will remain in the given folder structure. I understand that this may limit flexibility somewhat if the user is looking for something ultra-light, but sticking to this folder format means that the environment stays controlled and limits the variation across user environments to make reading/using/debugging/contributing as easy as possible. This may change in the future (may be I will implement a light and heavy version...)


Please feel free to raise pull requests with this repo, it's designed to be easily contributed to. For example, adding a new file format type for a LOCUST input equilibrium simply requires you to standardise your data according to the .data[] member dictionary, add the generic read/write functions (which you may have written before!) and add to the data_format "if" statement in the object read/write methods. Easy!

Got any burning questions? Want to feedback? Please email the author, Sam, at shw515@york.ac.uk. Alternatively, raise an issue on the CCFE github!






## Appendix


### Global Data Definitions & Variable Names

Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, here are all the different variable names used by different codes and file formats for similar quantities. '-' means they are not used/are empty/not written to.

#### Equilibrium (LOCUST_IO/GEQDSK/Equilibrium IDS/LOCUST):

    0D data
        /nw/-/nEQ_R            #number of points in R (x or width)
        /nh/-/nEQ_Z            #number of points in Z (y or height)
        /idum/-/IDUM          #number of spatial dimensions?
        /rdim//RDIM         #size of the R dimension in m
        /zdim//ZDIM         #size of the Z dimension in m
        r0/rcentr/.r0/RCENTR      #reference value of R
        b0/bcentr/.b0/BCENTR     #vacuum toroidal magnetic field at rcentr
        /rleft/    /RLEFT     #R at left (inner) boundary
        /zmid/      /ZMID    #Z at middle of domain (from origin)
        /rmaxis/  /Rmagh      #R at magnetic axis (O-point)
        /zmaxis/  /Zmagh      #Z at magnetic axis (O-point)
        /simag /  /PSI_magh      #poloidal flux psi at magnetic axis (Weber / rad)
        /sibry /  /SIBRY      #poloidal flux psi at plasma boundary (Weber / rad)
        /current/ /CURRENT      #plasma current [Amps]   
        /xdum/    /XDUM      #dummy variable - just contains zero
        /nbbbs/   /nb      #plasma boundary
        /limitr/   /IDUM     #wall boundary
    1D data
        /fpol/    /RBphih      #poloidal current function on uniform flux grid (1D array of f(psi)=R*B_tor [meter-Tesla])
        /pres/    /PRES      #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        /ffprime/  /FFPRIM     #workk1
        /pprime/   /PPRIM     #workk1
        /qpsi/      /QPSI    #q values on uniform flux grid
        /rlim/          /#r wall boundary
        /zlim/         / #z wall boundary
        /rbbbs/       /Rp  #r plasma boundary
        /zbbbs/      /Zp   #z plasma boundary
    2D data
        /psirz/    /psi_equil_h     #array (nx,ny) of poloidal flux on rectangular grid points (array of arrays)   




###[GEQDSK File Format](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)

The G EQDSK file contains information on P', FF', the flux on the
rectangular grid used, q, the boundary, and the limiter contour. A
right-handed cylindrical coordinate system (R, phi, Z) is used. The G EQDSK
provides information on the pressure , poloidal current function, q profile
on a uniform flux grid from the magnetic axis to the plasma boundary and the
poloidal flux function on the rectangular computation grid. Information on
the plasma boundary and the surrounding limiter contour is also provided. 

Variables 
---------

CASE: Identification character string 

NW: Number of horizontal R grid points 

NH: Number of vertical Z grid points 


Namelist OUT1:  

BCENTR: Vacuum toroidal magnetic field in Tesla at RCENTR 

CURRENT: Plasma current in Ampere 

FPOL: Poloidal current function in m-T, F = RBT on flux grid 

FFPRIM: FF'(psi) in (mT)2 / (Weber/rad) on uniform flux grid 

LIMITR: Number of limiter points (length of Rlim or Zlim)

NBBBS: Number of boundary points (length of RBBBS or ZBBBS)

NQPSI: 

PPRIME: P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid 

PRES: Plasma pressure in nt / m2 on uniform flux grid 

PSIZR: Poloidal flux in Weber/rad on the rectangular grid points 

QPSI: q values on uniform flux grid from axis to boundary 

RBBBS: R of boundary points in meter 

RCENTR: R in meter of vacuum toroidal magnetic field BCENTR 

RDIM: Horizontal dimension in meter of computational box 

RLEFT: Minimum R in meter of rectangular computational box 

RLIM: R of surrounding limiter contour in meter 

RMAXIS: R of magnetic axis in meter 

SIMAG: poloidal flux at magnetic axis in Weber/rad 

SIBRY: poloidal flux at the plasma boundary in Weber/rad 

ZBBBS: Z of boundary points in meter 

ZDIM: Vertical dimension in meter of computational box 

ZLIM: Z of surrounding limiter contour in meter 

ZMAXIS: Z of magnetic axis in meter 

ZMID: Z of center of computational box in meter (from origin)


Namelist BASIS  

  KPPFNC  KPPKNT  PPKNT  PPTENS  KFFFNC  KFFKNT  FFKNT 
  FFTENS  KWWFNC  KWWKNT  WWKNT  WWTENS  PPBDRY  PP2BDRY  KPPBDRY  KPP2BDRY 
  FFBDRY  FF2BDRY  KFFBDRY  KFF2BDRY  WWBDRY  WW2BDRY  KWWBDRY  KWW2BDRY 

Namelist CHIOUT

  SAISIL  SAIMPI  SAIPR  SAIIP  
  
Toroidal Current Density  The toroidal current JT related to P'(psi) and
FF'(psi) through  JT (Amp/m2) = R P'(psi) + FF'(psi)/R 
