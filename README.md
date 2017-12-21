# LOCUST_IO Package

This "package" is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. There are two sub-modules - input and output - which contain python objects, each representing LOCUST data.

It is designed such that individual file read/write methods can easily be swapped out if new methods/filetypes come and go - all that needs to be changed is the respective function (which would be the same amount of work as writing the function bespoke).







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











## Structure

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











Template things below


* Bullet point here

`highlight this writing`



```
this code will just be in a little box
```

[this_appears_as_a_url](www.with.this.address)


Table of Contents
-----------------

* [Title Here](#title-here)
  * [Sub Title Here](#sub-title-here)
