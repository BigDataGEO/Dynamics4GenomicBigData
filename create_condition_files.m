function create_condition_files()

set_paths_and_imports;

fprintf('\n');
display(['The pipeline analysis requires that you enter the accession number of the dataset you want to analyze.']);

fprintf('\n');
display(['The accession number is the alphanumeric code starting with ''GSE'' that identifies the dataset on the GEO.']);

fprintf('\n');
GEO_number = input(['Please enter the accession number of your dataset enclosed in single quotes (e.g., ''GSE52428''): ']);

try
  fprintf('\n');
  display(['Loading dataset. This can take some time, please wait...']);
  [geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number);
catch
  fprintf('\n');
  display(['Could not retrieve dataset ' GEO_number ' from the Gene Expression Omnibus.']);
  fprintf('\n');
  display(['This is possibly because the GEO refused the FTP connection or because the dataset does not exist.']);
  fprintf('\n');
  display(['Please download manually ' GEO_number '''s matrix to ' pwd '/GEO_cache/' GEO_number '.txt and try again.']);
  return;
end

Preprocessing_technique = 'Default';

index_of_analysis = 1;

prompt = ['\nThe pipeline analysis requires that you specify the samples from ' GEO_number '\nthat refer to your desired subject/condition and the time points.\n\nThis interactive interface will allow you to enter this information.\n\nPress Enter to proceed.  '];
input(prompt);

mkdir('Input');

while true
  
  [raw_gene_expression_array{index_of_analysis}, raw_time_points_array{index_of_analysis}, subject_name_array{index_of_analysis}, condition_array{index_of_analysis}, number_of_top_DRGs_considered_array{index_of_analysis}, sample_codes{index_of_analysis}] = prepare_data(geoStruct);

  raw_time_points_array{index_of_analysis} = strtrim(strcat(cellstr(num2str(raw_time_points_array{index_of_analysis})),' hours'));
  
  cd('Input');
  
  writetable(cell2table([num2cell(raw_time_points_array{index_of_analysis}) sample_codes{index_of_analysis}]), [GEO_number '_-_' strrep(condition_array{index_of_analysis}, ' ', '_') '_-_' num2str(number_of_top_DRGs_considered_array{index_of_analysis}) '.csv'], 'WriteRowNames', false, 'WriteVariableNames', false, 'Delimiter', ',');
  
  cd('..');
  
  fprintf('\n\n'); 
  display(['The information for the analysis of subject/condition "' condition_array{index_of_analysis} '" has been loaded successfully.']);
  
  prompt = '\nWould you like to also run another analysis with a different subject/condition? ([1 "yes", 0 "no"]) ';
  continue_analysis   = input(prompt);
  
  if(continue_analysis)
    prompt = '\nNow you will be required to enter the information (samples and time points) of the new subject/condition.\n\npress Enter to proceed. ';
    input(prompt);
  else
    fprintf('\n');
    display(['The input files for all the subjects/conditions you entered have been created and output to ' Dynamics4GenomicBigData_HOME 'Input/.']);
    break;
  end
  
  index_of_analysis = index_of_analysis + 1;
end

end


function [raw_gene_expression, raw_time_points, subject_name, condition, number_of_top_DRGs_considered, sample_codes] = prepare_data(geoStruct)
  
  [raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo] = get_data_from_geo_struct(geoStruct);
  
  [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs_considered, sample_codes] = capture_data(geoStruct, raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo);
  
  raw_time_points = raw_time_points';
  
end

function [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs, sample_codes] = capture_data(geoStruct, raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo)

  %Ask which condition to analyze
  tb = tabulate(titles);
  list_of_sample_record_titles = strcat(cellstr(arrayfun(@num2str, 1:length(tb(:,1)), 'UniformOutput', false))',' : ',tb(:,1));
  display(list_of_sample_record_titles);    
	
  condition_name = input(['Enter a name for the experimental condition: ']);
  display(sprintf(['\n\nAll the samples associated to %s are listed above. You must now specify the samples that correspond to condition "' condition_name '".\n'], geoStruct.Header.Series.geo_accession));
  display(sprintf('You must do this by entering a term or string that is common ONLY to the desired samples from the list above.\n'));
  prompt = ['Which samples would you like to include as part of condition "' condition_name '"? (enter the common string) '];
  common_string_identifying_subject_samples = input(prompt);
  indices_of_selected_sample_record_titles = ExtractSubjects(common_string_identifying_subject_samples, list_of_sample_record_titles);
  display(sprintf('\nThe samples that match your search terms are listed below.\n'));
  display(list_of_sample_record_titles(indices_of_selected_sample_record_titles));

  raw_gene_expression = raw_gene_expression_data_in_all_samples(:,indices_of_selected_sample_record_titles);
  
  sample_codes = geoStruct.Header.Samples.geo_accession(indices_of_selected_sample_record_titles)';
  
  display(sprintf('\nThe list below shows the information available for each one of the selected samples.\n\n'));
    
  %Display the Characteristics of the GEO series
  display(strcat(cellstr(arrayfun(@num2str, 1:length({metadata_for_all_samples{:,indices_of_selected_sample_record_titles(1)}}), 'UniformOutput', false))',' : ',{metadata_for_all_samples{:,indices_of_selected_sample_record_titles(1)}}'));

  %Find out where the time is
  prompt = 'Which item in the list corresponds to the time point? (Enter item number or -1 if none): ';
  index_of_time_row_in_sample_metadata = input(prompt);
  
  if index_of_time_row_in_sample_metadata > 0
    % Traditional case: time points can be read from characteristics matrix.
    format bank;
    raw_time_points    = {metadata_for_all_samples{index_of_time_row_in_sample_metadata,indices_of_selected_sample_record_titles}};
    if ~isempty(cell2mat(strfind(raw_time_points,'Baseline')))
      raw_time_points = strrep(raw_time_points, 'Baseline', '0');
    end
    raw_time_points    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), raw_time_points, 'UniformOutput', false));

    sane_check = 0;
    while sane_check == 0
      display(raw_time_points');
      prompt = '\nThese are all the time values measured in hours. Are they correct? (Enter 1 for "Yes" or 0 for "No") ';
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      raw_time_points = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	raw_time_points = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    %% Find out where the subject is
%      prompt = '\nWhich item in the list corresponds to the subject/condition? (format [1,2,3] or all) ';
%      index_of_subject_in_sample_metadata = input(prompt);
%      [subject_name,~,Subject] = unique({metadata_for_all_samples{index_of_subject_in_sample_metadata,indices_of_selected_sample_record_titles}});
    
    if iscellstr(common_string_identifying_subject_samples)
      subject_name = strjoin(common_string_identifying_subject_samples, '_');
    else
      subject_name = common_string_identifying_subject_samples;
    end
    Subject = repmat(1, 1, size(raw_time_points,2));
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(geoStruct.Header.Series.geo_accession, size(raw_gene_expression, 1));
    
  else
    % New case where time points must be read from title field or somewhere else.      
    timePointsMatrix = ExtractTimePoints(list_of_sample_record_titles(indices_of_selected_sample_record_titles));
      
    raw_time_points = [];
    timeUnit = 'N/A';
      
    if not(isempty(timePointsMatrix))
      raw_time_points = cell2mat(timePointsMatrix(:,1))';
      timeUnit = timePointsMatrix(1,2);
    end
	
    sane_check = 0;
    while sane_check == 0
      display(raw_time_points');
      prompt = ['These are all the time values measured in ' timeUnit{1} '. Are they correct? (Enter 1 for "Yes" or 0 for "No") '];
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      raw_time_points = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	raw_time_points = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    if iscellstr(common_string_identifying_subject_samples)
      subject_name = strjoin(common_string_identifying_subject_samples, '_');
    else
      subject_name = common_string_identifying_subject_samples;
    end
    Subject = repmat(1, 1, size(raw_time_points,2));
    
    [number_of_top_DRGs] = capture_top_number_of_DRGs(geoStruct.Header.Series.geo_accession, size(raw_gene_expression, 1));
  end

  condition = condition_name;
end

function [number_of_top_DRGs] = capture_top_number_of_DRGs(GEO_number, total_number_of_genes)
    prompt = ['Dataset ' GEO_number ' contains a total of ' num2str(total_number_of_genes) ' genes.\n\nEnter the number of top DRGs you want to consider in the analysis (or -1 to include them all): '];
    number_of_top_DRGs = input(['\n' prompt]);
    if(~isnumeric(number_of_top_DRGs) | number_of_top_DRGs < 1)
      number_of_top_DRGs = 3000;
    end
    number_of_top_DRGs = ceil(number_of_top_DRGs);
end

% Input argument subjectIDParts is either a string or an array of strings.
% Input argument Matrix is a matrix of strings.

% The function returns the row numbers of Matrix whose cells contain all of the strings in subjectIDParts that do not start with '~' and excluding those who contain at least one of the strings in subjectIDParts that start with '~'.


% subjectIDParts = {'BBB_CCC', '~min'}
% Matrix={'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_3435.3min_DDD_XXXX', 'AAAAA_BBB_CCC_65.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_656.3hrs_DDD_XXXX','AAAAA_BBB_CCC_965.3hrs_DDD_XXXX'}'

% Example output:
%  ans =
%  
%       1
%       3
%       4
%       5

function rowsWhereSubjectWasFound = ExtractSubjects(subjectIDParts, Matrix)

  subjectNamePartsToInclude = [];
  subjectNamePartsToExclude = [];
  
  % Function strjoin, used further below, needs to receive an array of strings. Therefore, if the function received only one string then it must be converted to an array of strings containing only the original one.
  if ischar(subjectIDParts)
    subjectIDParts = {subjectIDParts};
  end
  
  for i=1:length(subjectIDParts)
    % If this condition is true, then the string must be excluded.
    if(strcmp(subjectIDParts{i}(1:1),'~'))
      subjectNamePartsToExclude = [subjectNamePartsToExclude, {subjectIDParts{i}(2:length(subjectIDParts{i}))}];
    % Otherwise the string must be included.
    else
      subjectNamePartsToInclude = [subjectNamePartsToInclude, {subjectIDParts{i}}];
    end
  end
  
  rowsWhereSubjectWasFound = ExtractSubjectsWithInclusionsAndExclusions(subjectNamePartsToInclude, subjectNamePartsToExclude, Matrix);
end


function rowsWhereSubjectWasFound = ExtractSubjectsWithInclusionsAndExclusions(subjectNamePartsToInclude, subjectNamePartsToExclude, Matrix)
  rowsWhereSubjectWasFound = [];
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      regularExpressionToSearchFor = strjoin(subjectNamePartsToInclude,'.*');     
      startIndex = regexp(stringToSearch, regularExpressionToSearchFor);      
      % If the following condition is true, it means that stringToSearch at least contains the substrings that in subjectNamePartsToInclude.
      if (startIndex~=0)
	if(length(subjectNamePartsToExclude) == 0)
	  rowsWhereSubjectWasFound = cat(2, rowsWhereSubjectWasFound, i);
	else
	  stringToSearchShouldBeIncluded = 1;
	  for k=1:length(subjectNamePartsToExclude)
	    % if 1, then the string should be included. If [] then the string should not be included.
	    found = regexp(stringToSearch, subjectNamePartsToExclude{k});
	    
	    % If this condition is true, then the string does not contain the substring to be excluded. Thus the string should be included in the output list.
	    if isempty(found)
	      % Nothing to do here.
	    else % Otherwise the string contains the substring to be excluded. Thus the string should not be included in the output list.
	      stringToSearchShouldBeIncluded = 0;
	      break;
	    end
	  end
	  if(stringToSearchShouldBeIncluded)
	    rowsWhereSubjectWasFound = cat(2, rowsWhereSubjectWasFound, i);
	  end
	end
      end
    end
  end
  rowsWhereSubjectWasFound = rowsWhereSubjectWasFound.';
end

function [dat, title, Info, PInfo] = get_data_from_geo_struct(geoStruct)

%% Read Data from the GEO website and if cell files are avaiable perform RMA normalisation

%Request that the User input the GSE number 
%prompt = 'Enter a string specifying a unique identifier for a GEO Sample (GSM), Data Set (GDS), Platform (GPL),or Series (GSE) record in the GEO database.  ';
%GEO_number = input(prompt);  
    
%Extract the required information on the Study
row_identifiers_of_gse_matrix = rownames(geoStruct.Data);
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

PInfo = {Data_type, Normalisation_technique, Organism};

Info = [];
names = fieldnames(geoStruct.Header.Samples);
index = find(~cellfun(@isempty,strfind(names,'ch')));
index2 = find(~cellfun(@isempty,strfind(names,'description')));

if(isfield(geoStruct.Header.Samples,'characteristics_ch1'))
  Info = geoStruct.Header.Samples.characteristics_ch1; %Obtain the Charateristics of the Data
elseif(~isempty(index))
  tmp   = struct2cell(geoStruct.Header.Samples);
  tmp   = tmp(index);
  Info  = vertcat(tmp{:});
  if(~isempty(index2))
      tmp   = struct2cell(geoStruct.Header.Samples);
      tmp   = tmp(index2);
      Info  = [Info;vertcat(tmp{:})];
  end
end

title =[];
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
% row_identifiers_of_gse_matrix    = rownames(Expression);
% sam_id = colnames(Expression);
% 
% end

end


function timePoints = InputTimePointsManually()  
  timePoints = [];
  prompt = 'Please enter the time points as a list of text strings enclosed in curved brackets (e.g. {''60 min'', ''2 hr'',''3 hr''}):\n';
  timePoints = ExtractTimePoints(input(prompt).');
end
