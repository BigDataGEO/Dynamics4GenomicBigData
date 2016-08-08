close all;
clc;
clear;

path = strcat(pwd,'/');

%if first time running on computer run line 9 and 10 
%cd([path,'SBEToolbox-1.3.3\'])
%install

%Add Paths
addpath(path)
addpath(genpath([path,'fdaM/']))
addpath(genpath([path,'SBEToolbox-1.3.3/']))

[~,GEO_num] = xlsread('GEO_list.xlsx');

py.importlib.import_module('DAVIDWS');

for i = 1:size(GEO_num,1)

  GEO_number              = char(GEO_num(i));
  Preprocessing_technique = 'Default';
  [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);

  [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);
  [Data_2, Subject_2, Pos_2, str_ind_2, pr_ind_2, tb_2, Subject_name_2] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);  
  
  [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
  con = strrep(con,' ','_');
  con = strrep(con,'/','_');
  con = strrep(con,'.','_');
  flder = strcat(path,'Results/',GEO_number,'/',con);
  mkdir(flder)
  cd(flder)
  options = struct('format','html','outputDir',flder,'showCode',true);
  Paper2;
  close all;
  cd(path);
end





