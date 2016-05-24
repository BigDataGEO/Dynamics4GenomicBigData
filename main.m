close all;
clc;
clear;

%Modify the paths for your computer
path = '/home/jramirez/Dynamics4GenomicBigData/';
%  path = '/Users/michellecarey/Dropbox/Dynamics4GenomicBigData/';

%if first time running on computerrun line 9 and 10 
%cd([path,'SBEToolbox-1.3.3\'])
%install


%Add Paths
addpath(path)
addpath(genpath([path,'fdaM/']))
addpath(genpath([path,'SBEToolbox-1.3.3/']))

[~,GEO_num] = xlsread('GEO_list.xlsx');

for i = 1:size(GEO_num,1)

%User Inputs/options
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
no_iter = 0;
while cont == 1
    
%Ask which condition to analyze
tb = tabulate(titles);
dis = strcat(cellstr(arrayfun(@num2str, 1:length(tb(:,1)), 'UniformOutput', false))',' : ',tb(:,1));
display(dis);    
    
prompt = 'Which probe sets would you like to analyze? (format common string)';
str_ind = input(prompt);
%  The following two lines were replaced with the one line further below in order to allow to user to enter a subject's name by parts and thus be able to look it up when this name appears split in the dataset's record.
%  index  = strfind(dis,str_ind);
%  pr_ind = find(~cellfun(@isempty,index));
pr_ind = ExtractSubjects(str_ind, dis);
display(dis(pr_ind));

%Get Data for that condtion
Data     = Data_GEO(:,pr_ind);

if not(isempty(Info))
  if(no_iter==0)
  %Display the Characteristics of the GEO series
  display(strcat(cellstr(arrayfun(@num2str, 1:length({Info{:,pr_ind(1)}}), 'UniformOutput', false))',' : ',{Info{:,pr_ind(1)}}'));

  %Find out where the time is
  prompt = 'Which row has the time indicator in it? (format [1,2,3])';
  tm_ind = input(prompt);
  Pos    = {Info{tm_ind,pr_ind}};
  if ~isempty(cell2mat(strfind(Pos,'Baseline')))
      Pos = strrep(Pos, 'Baseline', '0');
  end
  Pos    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), Pos, 'UniformOutput', false));

  display(Pos);
  prompt = 'These are time values are they right?';
  sane_check = input(prompt);
  if(sane_check==0)
      display('email/hipchat carey or michelle with your GSE number')
  end

  %% Find out where the subject is
  prompt = 'Which row has the subject/condition indicator in it? (format [1,2,3] or all)';
  su_ind = input(prompt);
  [Subject_name,~,Subject] = unique({Info{su_ind,pr_ind}});
  end

  display(Subject);
  prompt = 'These are subject values are they right?';
  sane_check = input(prompt);
  if(sane_check==0)
      display('email/hipchat carey or michelle with your GSE number')
  end

  no_iter = no_iter+1;

  %% Set up Directory of Folders
  [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
  con = strrep(con,' ','_');
  flder = strcat(path,'Results/',GEO_number,'/',con);
  mkdir(flder)
  cd(flder)

  conditions_analyzed{cont} = con;

  %% Run Pipeline steps

  %Run pipeline each subject at a time
  options = struct('format','html','outputDir',flder,'showCode',true);
  publish('Paper.m',options)
  web('Paper.html', '-browser')

  save(strcat(GEO_number,con,date))

  prompt = 'Would you like to continue to the next subject/condition? ([1 "yes", 0 "no"])';
  cont   = input(prompt);

  close all;
  
  
  
  
  
else
%    if(no_iter==0)
%    %Display the Characteristics of the GEO series
%    display(strcat(cellstr(arrayfun(@num2str, 1:length({Info{:,pr_ind(1)}}), 'UniformOutput', false))',' : ',{Info{:,pr_ind(1)}}'));

  %Find out where the time is
%    prompt = 'Which row has the time indicator in it? (format [1,2,3])';
%    tm_ind = input(prompt);
%    Pos    = {Info{tm_ind,pr_ind}};
%    if ~isempty(cell2mat(strfind(Pos,'Baseline')))
%        Pos = strrep(Pos, 'Baseline', '0');
%    end
%    Pos    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), Pos, 'UniformOutput', false));

%    Pos = ExtractTimePoints(titles);
  Pos = ExtractTimePoints(dis(pr_ind));
  display(Pos);
  prompt = 'These are time values are they right?';
  sane_check = input(prompt);
  if(sane_check==0)
      display('email/hipchat carey or michelle with your GSE number')
  end

  %% Find out where the subject is
%    prompt = 'Which row has the subject/condition indicator in it? (format [1,2,3] or all)';
%    su_ind = input(prompt);
%    [Subject_name,~,Subject] = unique({Info{su_ind,pr_ind}});
%    end

  Subject_name = str_ind;
  Subject = repmat(1, 1, size(Pos,2));

  display(Subject);
  prompt = 'These are subject values are they right?';
  sane_check = input(prompt);
  if(sane_check==0)
      display('email/hipchat carey or michelle with your GSE number')
  end

  no_iter = no_iter+1;

  %% Set up Directory of Folders
  [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
  con = strrep(con,' ','_');
  con = strrep(con,'/','_');
  flder = strcat(path,'Results/',GEO_number,'/',con);
  mkdir(flder)
  cd(flder)

  conditions_analyzed{cont} = con;

  %% Run Pipeline steps

  %Run pipeline each subject at a time
  options = struct('format','html','outputDir',flder,'showCode',true);
  publish('Paper.m',options)
  web('Paper.html', '-browser')

%    disp(['I am saving the following: ', GEO_number, ', ', con,', ', date]);
  save(strcat(GEO_number,con,date))

  prompt = 'Would you like to continue to the next subject/condition? ([1 "yes", 0 "no"])';
  cont   = input(prompt);

  close all;

end
end

%% Create Manuscript

prompt = 'Which condtions would you like the manuscript to include ? (format [1,2,3])';
cond   = input(prompt);


% ind = pr_ind(1);
% 
% title_string = char(title(ind));
% words = eval(strrep(['{',sprintf('{''%s''},', title_string),'}'],' ',''' '''));
% display([title_string,geoStruct.Header.Samples.organism_ch1(ind)])
% prompt = 'What words would you like to search in pubmed? (format common string in cell {...})';
% wrd_ind = input(prompt);

%% Create latex document
options = struct('format','latex','showCode',false,'outputDir',flder,...
    'stylesheet','document.xsl');
publish('Paper_dum.m',options)

%opens your web browser and looks for papers related to your data set. Read
%these papers.
google(GEO_number,'scholar');

end

