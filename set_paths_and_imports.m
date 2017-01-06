close all;
clc;
clear;

global Dynamics4GenomicBigData_HOME;
Dynamics4GenomicBigData_HOME = strcat(pwd,'/');

global pipeline_version;
pipeline_version = 'Pipeline Version 1.41';

%  cd([Dynamics4GenomicBigData_HOME,'SBEToolbox-1.3.3\'])
%  install

%Add Paths
addpath(Dynamics4GenomicBigData_HOME);
addpath(genpath([Dynamics4GenomicBigData_HOME,'fdaM/']));
addpath(genpath([Dynamics4GenomicBigData_HOME,'SBEToolbox-1.3.3/']));

py.importlib.import_module('DAVIDWS');

warning('off','all');
clc;
clc;
