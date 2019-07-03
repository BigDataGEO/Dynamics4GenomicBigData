close all;
%  clc;
%  clear;

global Dynamics4GenomicBigData_HOME;
Dynamics4GenomicBigData_HOME = pwd;

global pipeline_version;
pipeline_version = 'Pipeline Version 1.64';

%Add Paths
addpath(Dynamics4GenomicBigData_HOME);
addpath(genpath([Dynamics4GenomicBigData_HOME filesep 'lib' filesep 'fdaM']));
addpath(genpath([Dynamics4GenomicBigData_HOME filesep 'lib' filesep 'SBEToolbox-1.3.3']));

addpath(genpath([Dynamics4GenomicBigData_HOME filesep 'lib' filesep 'FTP']));
rehash toolboxcache

%  py.importlib.import_module('DAVIDWS');

warning('off','all');
%  clc;
%  clc;
