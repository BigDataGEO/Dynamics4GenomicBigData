function main()

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

    % -----------------------------------------------------------------------
    %  Data Preprocessing
    % -----------------------------------------------------------------

    %Retrieve Data From GEO
    [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);

    %Some data sets have more than one condition. By default the Pipeline will
    %analyze the first condtion only. If you would like to specify another
    %condition which you would like to run the pipeline please specify condtion below
    conditions_analyzed = cell(20);
    cont = 1;
    while cont == 1
    
      [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);
    
      [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
      con = strrep(con,' ','_');
      con = strrep(con,'/','_');
      con = strrep(con,'.','_');
      flder = strcat(path,'Results/',GEO_number,'/',con);
      mkdir(flder)
      cd(flder)

      conditions_analyzed{cont} = con;

      %Run pipeline each subject at a time
      options = struct('format','html','outputDir',flder,'showCode',true);
      publish('Paper.m',options);
      web('Paper.html', '-browser');

      prompt = '\nWould you like to start a new analysis with a different set of samples (i.e., a different subject/condition)? ([1 "yes", 0 "no"]) ';
      cont   = input(prompt);

      close all;
      cd(path);
    end %while cont == 1

    %% Create Manuscript
    prompt = '\nWhich subjects would you like the manuscript to include? (format [1,2,3]) ';
    cond   = input(prompt);
    
    % If running on Windows.  
    if ispc()
      options = struct('format','doc','showCode',false,'outputDir',flder,'stylesheet','document.xsl');
      publish('Paper_dum.m',options);
    % elseif isunix() % If running on Unix/Linux.
    else
      options = struct('format','latex','showCode',false,'outputDir',flder,'stylesheet','document.xsl', 'imageFormat', 'pdf');
      publish('Paper_dum.m',options);
      copyfile('latex/bibliography.bib', flder);
      copyfile('latex/plos2015.bst', flder);
      
      % The following line compiles the .tex file into a .pdf.
      % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
      [x, y]=system(['./latex/compile.sh ' flder]);
    end
    google(GEO_number,'scholar');
  end

