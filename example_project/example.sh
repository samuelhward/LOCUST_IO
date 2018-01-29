#!/usr/bin/env bash

imasdb test #set the location of the IDS database

#kinplot -r run_number -s shot_number #make some nice IDS plots
#scenplot -r run_number -s shot_number

./example.py

#point LOCUST to files here using symlink

#make FLAGS='' #compile and run LOCUST
#locust