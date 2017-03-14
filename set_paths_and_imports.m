close all;
clc;
clear;

global Dynamics4GenomicBigData_HOME;
Dynamics4GenomicBigData_HOME = strcat(pwd,'/');

global pipeline_version;
pipeline_version = 'Pipeline Version 1.56';

%Add Paths
addpath(Dynamics4GenomicBigData_HOME);
addpath(genpath([Dynamics4GenomicBigData_HOME,'lib/fdaM/']));
addpath(genpath([Dynamics4GenomicBigData_HOME,'lib/SBEToolbox-1.3.3/']));

addpath(genpath([Dynamics4GenomicBigData_HOME,'lib/FTP/']));
rehash toolboxcache

%  py.importlib.import_module('DAVIDWS');

warning('off','all');
clc;
clc;
