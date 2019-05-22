# LOCUST_IO Package

This "package" is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. 

The code aims to be...

* __simple__ - It is designed to abstract the user away from LOCUST - operations such as converting from one data format to another can be done in two lines. Examples on how to use the API are included below, along with an example project in *LOCUST_IO/docs/example_project* and multiple tutorials in the *LOCUST_IO/docs/HOWTO.md* guide. Transparency is also key and the code is heavily commented.

* __extensible__ - Contributing to the project is simple: to add a new input/cache/output type, for example, simply copypaste the read/write functions and use the data dictionary to match up your variables. Many functions aren't made class methods to maximise flexibility (and those which take LOCUST_IO objects tend to accept generic dictionaries via the __ getitem __ method). Adding your own high level plot scripts, analysis tools and workflows means just adding them to the respective files in *src/*. 

* __portable__ - LOCUST_IO has all the infrastructure needed to encapsulate your automated pre/post processing workflow and can be ran 'out of the box'. Integration with LOCUST simply means cloning the latest version into *LOCUST/* folder!


Got any burning questions? Want to feedback? Please raise an issue here on github! or email me at samuel.ward@york.ac.uk - and for quick help check the docstrings!



Table of Contents
-----------------

* [Requirements](#requirements)
* [Getting Started](#getting-started)
* [Usage](#usage)
* [Licence](#licence)




## Requirements

Tested with:

* IMAS v3.16.0
* Python ≥v3.6
    * numpy ≥v1.14.2
    * matplotlib ≥v2.0.0
    * scipy ≥v0.17.0
    * h5py ≥2.6.0




## Getting Started

* Set up a folder anywhere within LOCUST_IO e.g. *LOCUST_IO/my_project* (this is where you will create all of your analysis scripts)
* Copy a context.py file from *LOCUST_IO/docs/example_project* to *LOCUST_IO/my_project*
* Copy LOCUST input/cache/output files to *LOCUST_IO/data/input_files*, *LOCUST_IO/data/cache_files* or *LOCUST_IO/output_files* respectively
* import context and off you go!






## Usage


As well as the included *LOCUST_IO/docs/example_project/* some basic usage is outlined below:

```python
import context #tell LOCUST_IO where we are
from classes.input_classes.equilibrium import Equilibrium #import classes which encapsulate LOCUST inputs, e.g. an equilibrium

#to read a GEQDSK from LOCUST_IO/data/input_files/locust_run_1/:
my_equilibrium=Equilibrium(ID='ID_tag_describing_this_equilibrium - mandatory!',data_format='GEQDSK',filename='locust_run_1/some.eqdsk') 
#my_equilibrium now holds all the data in one object
#take a quick look at the equilibrium and its data
my_equilibrium.look()                               

#to initialise empty equilibrium to fill later with the read_data() (must always specify an ID):
my_equilibrium=Equilibrium(ID='a blank equilibrium') 
#read data at a later time from GEQDSK from input_files/
my_equilibrium.read_data(data_format='GEQDSK',filename='some.eqdsk',property1='made using EFIT')
#dump equilibrium to IMAS IDS format 
my_equilibrium.dump_data(output_data_format='IDS',shot=1,run=1) 


#you can set individual pieces of data with the .set() method
#this will overwrite the default data format, which is numpy array:
#set multiple values simultaneously
my_equilibrium.set(nw=5,fpol=[1,2,3,4])  
#equally
my_equilibrium.set(**some_dict)                     


#your input/output objects can also be copied using the .copy() method:
#copy all data from one object to another
my_equilibrium.copy(some_other_equilibrium)
#copy specific fields                                        
my_equilibrium.copy(some_other_equilibrium,'B_field_R','some_key','some_other_key')  


#to get a quick glimpse of what you're working with, LOCUST_IO can also plot input/output data: 
my_equilibrium.plot()                                         
#(you can also stack plots onto the same axis object with the ax arguement - see example_project)
my_equilibrium.plot(ax=some_ax,fig=some_fig)                                         


#to check what data two objects share, use .compare():
my_equilibrium.compare(another_equilibrium,verbose=True)               


#to check if your object contains enough information for running LOCUST, use .run_check():
my_equilibrium.run_check(verbose=True)                                            


#you can also calculate new pieces of data using methods or functions in the processing folder
my_equilibrium.B_calc()         
```



## Licence:


    LOCUST_IO
    Copyright (C) 2019  Samuel Ward

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