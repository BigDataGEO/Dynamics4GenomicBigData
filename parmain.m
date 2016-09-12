set_paths_and_imports;
  
[~,GEO_num] = xlsread('GEO_list.xlsx');

GEO_number = char(GEO_num(1));
Preprocessing_technique = 'Default';

index_of_analysis = 1;

prompt = ['\nThe pipeline analysis requires that you specify the samples from ' GEO_number '\nthat refer to your desired subject/condition and the time points.\n\nThis interactive interface will allow you to enter this information.\n\nPress Enter to proceed.  '];
input(prompt);

while true
  
  [list_of_genes_array{index_of_analysis}, raw_gene_expression_array{index_of_analysis}, raw_time_points_array{index_of_analysis}, subject_name_array{index_of_analysis}, condition_array{index_of_analysis}, gene_ID_type_array{index_of_analysis}, number_of_top_DRGs_considered_array{index_of_analysis}] = step_1(GEO_number);
  
  fprintf('\n\n'); 
  display(['The information for the analysis of subject/condition "' condition_array{index_of_analysis} '" has been entered successfully.']);
  
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

total_number_of_analyses_to_run = index_of_analysis;

parfor index_of_analysis=1:total_number_of_analyses_to_run

  list_of_genes = list_of_genes_array{index_of_analysis};
  raw_gene_expression = raw_gene_expression_array{index_of_analysis};
  raw_time_points = raw_time_points_array{index_of_analysis};
  subject_name = subject_name_array{index_of_analysis};
  condition = condition_array{index_of_analysis};
  gene_ID_type = gene_ID_type_array{index_of_analysis};
  number_of_top_DRGs_considered = number_of_top_DRGs_considered_array{index_of_analysis};

  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
  
  pworkspace = 'Paper.mat';
  
  if ispc() % If running on Windows.
    options = struct('format','doc','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl');
    

  elseif isunix() % If running on Unix/Linux.
    options = struct('format','latex','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl', 'imageFormat', 'png');
    
    parsave('Paper.mat', GEO_number, list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered, Dynamics4GenomicBigData_HOME, path);
    
%      publish('Paper.m',options);
    [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, true);

    copyfile([Dynamics4GenomicBigData_HOME, '/latex/bibliography.bib'], output_folder);
    copyfile([Dynamics4GenomicBigData_HOME, '/latex/plos2015.bst'], output_folder);
    
    delete('Paper.mat');

  end
  
  cd(Dynamics4GenomicBigData_HOME);
end

fprintf('\n');
display('The analysis is complete for all the subjects/conditions.');

fprintf('\n');
display(['All results have been output to folder ' Dynamics4GenomicBigData_HOME 'Results/' GEO_number '/']);
