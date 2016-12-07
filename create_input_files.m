function create_input_files()

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
  display(['Please download manually ' GEO_number '''s matrix to ' pwd '/GEO_cache/' GEO_number '.soft and try again.']);
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
    display('The analysis will begin now for all the subjects/conditions you entered.');
    break;
  end
  
  index_of_analysis = index_of_analysis + 1;
end

end