# LOCUST_IO Package

This package is designed to process input and output for the LOCUST-GPU code, developed by Rob Akers at CCFE. There are two sub-modules - input and output - which contain python objects, each representing LOCUST data.

It is designed such that individual file read/write methods can easily be swapped out if new methods/filetypes come and go - all that needs to be changed is the respective function (which would be the same amount of work as writing the function bespoke).






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



## Advanced Usage
* Import the classes submodule by...






## Structure

The overall layout consists of the main subpackage "classes" which takes data from the "input_files" and "output_files" folders, as well as individual folders for projects which draw from these folders. Every project folder, such as "testing" and "examples" must contain a context.py to ensure the package can see its folders and files. 

### Classes Subpackage

#### Input Module

#### Output Module

#### Support Module


### Testing Subpackage

#### Input Class Tests

#### Output Class Tests





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
