# HOW TO

This is a collection of tutorials - enjoy!


Table of Contents
-----------------

* [Combine beam depositions](#Combine-input-particle-lists)
* [Create a workflow](#Create-a-workflow)
* [Use mars_builder](#Use-mars_builder)
* [Use LOCUST_edit_var](#Use-LOCUST_edit_var)


## Combine beam depositions

```python
import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd 

bd_combined=bd(ID='total combined beam depo') #create an empty beam deposition

#if we have multiple beam depositions in list 'beam_depos_to_combine'
for beam_depo in beam_depos_to_combine: 
    bd_combined.combine(beam_depo)    
```


## Create a workflow

A simple workflow platform is implemented in LOCUST_IO, some quick demos are below. Essentially all that is happening is a wrapping of a dispatch table. If a command fails, the next command will be tried - this way post-run cleanups can still take place.

```python
import context
from run_scripts.workflow import Workflow as wf #import the workflow template

class my_custom_workflow(wf):
    def __init__(self,message,*args,**kwargs):
        super().__init__() #must call the Workflow base constructor to properly initialise a workflow
        
        #do something with *args and **kwargs here
        self.message=message 

        #now declare workflow stages defined below

        #use .add_command to make functions visible to Workflow base class
        #here we define the name we use to call the command, the corresponding function and the position in the overall workflow 
        self.add_command(command_name='print a message',command_function=self.workflow_stage_1,position=1) 
        self.add_command(command_name='print finish message',command_function=self.workflow_stage_2,position=2) 

    #now define a step in the workflow
    def workflow_stage_1(self,*args,**kwargs): #best practice keep *args,**kwargs in workflow stage definitions
        print(self.message)

    def workflow_stage_2(self,*args,**kwargs): #print completion message!
        print("workflow completed!")

#test time
a_workflow=my_custom_workflow(message="this is a workflow")
a_workflow.run_command('print a message') #use the .run_command() to run an individual stage of a command
a_workflow.run() #use the .run() command to run all stages in the workflow
```

Builds and Environments are also classes I add to make Workflows more powerful. Builds define a particular state of a code piece of code that a Workflow might execute - the commit hash, source code modifications, compile flags. It also adds some methods to automate this code compilation.

Environments contain all the information needed to describe a system environment e.g. modules to load, environment variables to define. Both Builds and Workflows can contain an Environment to describe the build and runtime environments respectively. A selection of defaults for various systems is also included in LOCUST_IO.


```python
import context
from run_scripts.workflow import Workflow as wf
from run_scripts.environment import Environment as ev
from run_scripts.build import Build as bu

class better_workflow(wf):
    def __init__(self,*args,**kwargs):
        super().__init__() #must call the Workflow base constructor to properly initialise a workflow

        #say we want to build and execute some external code called some_code in our workflow
        
        #first we can define the build version of some_code that this workflow wants to execute
        self.build=bu(system_name='CUMULUS') #system_name will define the build environment i.e. which system you are on and which modules to load etc, you can also define your own
        self.build.flags_add(TOKAMAK=1,STDOUT=True) #add some compile flags
        self.build.source_code_mods_add(source_code_filename='some_code.f90',some_variable_in_some_code=5,some_string_in_some_code="'a_string'") #if some_code is in the fortran language, we can modify the default declaration of variables
        
        #now define the runtime environment for this workflow
        self.environment=ev(system_name='CUMULUS') 

        #attach commands as before
        self.add_command(command_name='compile',command_function=self.workflow_stage_1,position=1) 
        self.add_command(command_name='print finish message',command_function=self.workflow_stage_2,position=2) 
        print(self.commands) #print commands in order

    def workflow_stage_1(self,*args,**kwargs):
        self.build.make(directory='/source/code/directory',clean=True) #execute make clean
        self.build.make(directory='/source/code/directory',clean=False) #modifies source code with edits and compiles code with flags we have defined

    def workflow_stage_2(self,*args,**kwargs): #print completion message!
        print("workflow completed with following environment:")
        self.environment.display() #add extra step to display the environment
```

There are a number of advantages to handling chains of functions in this way. For example, you can ask workflows to execute stages based on a 'string', but that string can remain the same whilst the function behaviour changes depending on the workflow that has been constructed. Take the two workflows we just made as an example:

```python
some_workflow_1.run_command('print finish message') #these functions do different things
some_workflow_2.run_command('print finish message') #but are called in the same way
```

Some examples of workflows are included in LOCUST_IO, namely LOCUST_run and MARS_builder_run. The former adds functions to send off a single run of LOCUST and the latter can be used to automatically trigger mars_builder to transform 3D field data.  

## Use mars_builder

mars_builder is a piece of FORTRAN90 written by Rob Akers which aims to process MARS-F output. More specifically, it linearly combines separate MARS-F fields by interpolating them onto a cylindrical RZ grid to arbitrary precision. You can call mars_builder from within Python using the script MARS_builder_run.py. Some commonly used flags and variables are:

| flag    | description                                                                                            |
|---------|--------------------------------------------------------------------------------------------------------|
| TOKAMAK | integer value of desired tokamak (according to prec_mod.f90) e.g. -DTOKAMAK=1 corresponds to ITER      |
| COILROW | integer value to select individual coil row to dump to file (floors contributions of other coil rows)  |
| MATCH   | use to generate optimal RMP settings (e.g. phase, amplitude) to match field defined in variable `mtch` |

| variable | description                                                                                                   |
|----------|---------------------------------------------------------------------------------------------------------------|
| `root`   | prepends this string to all written filepaths                                                                 |
| `file`   | target 3D field filename string                                                                               |
| `TAIL`   | array holding strings to append to `file` in the case of multiple 3D field files - usually denoting coil rows |

## Use LOCUST_edit_var

LOCUST_edit_var is a function which allows you to directly edit FORTRAN90 source code to change the starting value of a declared variable. To change this:

```
character( len=10057 )         :: file = 'some/nice/& !variables declared like this
                                &file/somewhere' !delete this comment
```

to this:

```
character( len=10057 )         :: file = 'some/other/file/' !variables declared like this
```

you could do this:


```python
source_code_mods={}
source_code_mods['file']="'some/other/file'" 
run_scripts.LOCUST_edit_var.LOCUST_edit_var(filepath_in=some_filepath,filepath_out=some_filepath,**source_code_mods)
```
noting the double brackets (since whatever is inside the brackets is interpreted literally). Comments on the first line of the variable declaration are preserved, and the entire declaration is condensed onto a single line.
