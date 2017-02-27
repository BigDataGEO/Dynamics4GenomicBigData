#!/bin/bash

# The following line needs to be updated to meet the local installation.
export PATH=$PATH:/usr/local/MATLAB/R2016a/bin/

matlab -nodesktop -r "run compare.m;exit"
