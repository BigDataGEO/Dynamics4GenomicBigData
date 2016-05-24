function [dat,gid,title,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique)

%% Read Data from the GEO website and if cell files are avaiable perform RMA normalisation

%Request that the User input the GSE number 
%prompt = 'Enter a string specifying a unique identifier for a GEO Sample (GSM), Data Set (GDS), Platform (GPL),or Series (GSE) record in the GEO database.  ';
%GEO_number = input(prompt);

%Obtain information from the GEO website 
geoStruct = getgeodata(GEO_number);

%Extract the required information on the Study
gid    = rownames(geoStruct.Data);  %Gene Id
sam_id = colnames(geoStruct.Data);  %Sample Id

%Supplementary_file = geoStruct.Header.Series.supplementary_file; %http for supplementary_files
%cell_file          = strfind(Supplementary_file,'.tar'); %Are the supplementary_files cell files

Data_type = unique(geoStruct.Header.Samples.type);  %Data Type 

Normalisation_technique = [];
if(isfield(geoStruct.Header.Samples,'data_processing'))
Normalisation_technique = unique(geoStruct.Header.Samples.data_processing); %Normalisation technique 
end

Organism = [];
if(isfield(geoStruct.Header.Samples,'organism_ch1'))
Organism = unique(geoStruct.Header.Samples.organism_ch1); %Organism
end

PInfo = {Data_type,Normalisation_technique,Organism};

%  2016-05-16 (JCRI): Added a conditional block to check first if the information obtained from the GEO database includes the data characteristics of the queried dataset.
if isfield(geoStruct.Header.Samples,'characteristics_ch1')
  Info = geoStruct.Header.Samples.characteristics_ch1; %Obtain the Charateristics of the Data
else
  Info = []; %Obtain the Charateristics of the Data
end

if(isfield(geoStruct.Header.Samples,'title'))
title = geoStruct.Header.Samples.title;
end
% for i = 1:size(Info,1)
%     %if(length(unique(Info(i,:)))==1)
%     %    Info_DC(k,1) = unique(Info(i,:));
%     %    k=k+1;
%     if(mstrfind(char(Info{i}),{'time','day'})) %Obtain times if its in Charateristics
%             for j = 1:length(Info(i,:))
%                 tmp = regexp(char(Info(i,j)),'\d*\.?\d*','match');
%                 if(isempty(tmp))
%                     time(j) = -1;
%                 else
%                     time(j) = str2num(char(tmp));
%                 end
%                 if(strcmp(char(Info(i,j)),'time point: Baseline'))
%                     time(j) = 0;
%                 end
%             end
%     elseif (~isempty(strfind(char(Info{i}),'subject'))) %Obtain subject if its in Charateristics
%             for j = 1:length(Info(i,:))
%                 tmp = regexp(char(Info(i,j)),'\d+','match');
%                 subject(j)  = str2num(char(tmp(end)));
%             end
%     elseif (~isempty(strfind(char(Info{i}),'virus')) | ~isempty(strfind(char(Info{i}),'vaccine'))) %Obtain Virus/Vacine if its in Charateristics 
%            condition = unique(lower(Info(i,:)));
%            for j = 1:length(condition)
%            ind = find(strcmpi(lower(Info(i,:)), condition(j)));
%            condition_ind(ind) = j;
%            end
%     end
% %    else
% %    Info_C{i} = unique(Info(i,:));
% %    for j = 1:length(Info_C{i})
% %        ind = find(strcmpi(Info(i,:), Info_C{i}(j)));
% %        Chara(l,ind) = j;
% %        l=l+1;
% %    end
% %    end
% end

% if(isempty(subject))
%     des = geoStruct.Header.Samples.description;
%     for j = 1:length(des)
%         tmp = regexp(des{j},'subject (\d+)','tokens','once');
%         subject(j) = str2num(char(tmp));
%     end
% end
% 
% 
% if(isempty(time))
%     des = geoStruct.Header.Samples.description;
%     for j = 1:length(des)
%         tmp = regexp(des{j},'time (\d+)','tokens','once');
%         subject(j) = str2num(char(tmp));
%     end
% end
% 
% 
% if(isempty(condition_ind))
%     condition_ind(1:size(Info,2)) = 1;
%     condition = geoStruct.Header.Series.title;
% end

%Extract the processed data if the cell files are not avaiable

%% Read the Cell Files (Read from Cell Files if they exist) 

%if(isempty(cell_file) || Preprocessing_technique=='Default')
     dat    = double(geoStruct.Data);
%else 
    
% FTPobj = ftp('ftp.ncbi.nlm.nih.gov');    
% t=strsplit(Supplementary_file,'gov');
% filepath=t{2};
% mget(FTPobj,filepath);
% 
% % Performs RMA background adjustment, quantile normalization, and summarization procedures 
% % on the PM probe intensity values, and returns a DataMatrix object,
% % containing the metadata and processed data.
% 
% wd = pwd;
% 
% foldername = [wd,'/pub/geo/DATA/supplementary/series/',GEO_number,'/'];
% 
% filename = [GEO_number,'_RAW.tar'];
% 
% untar([foldername,filename])
% 
% gunzip('*.gz');

%if (Preprocessing_technique=='RMA')
% Expression = affyrma('*', 'GPL9188_Hs133Av2_Hs_ENTREZG.cdf');
%elseif(Preprocessing_technique=='GCRMA')
% Expression = affygcrma('*', 'GPL9188_Hs133Av2_Hs_ENTREZG.cdf');
%elseif(Preprocessing_technique=='GCRMA') 
% dmwrite(Expression, 'Data.txt')
% 
% dat    = double(Expression);
% gid    = rownames(Expression);
% sam_id = colnames(Expression);
% 
% end

end