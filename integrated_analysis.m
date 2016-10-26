function integrated_analysis()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  output_dir = parts{1};

  inputData = readtable(input_file_name, 'Delimiter', ',');
  inputData = table2cell(inputData);
  
  general_comparison_folder = [Dynamics4GenomicBigData_HOME, 'Results/Comparison/', output_dir];
  mkdir(general_comparison_folder);
    
  statistics_of_analyses = {'Series', 'Condition', '# of time points', '# of DRGs', '# of GRMs'};
  
  DRGs = {};
  GRMs = {};

  for i = 1:size(inputData,1)  
    [gene_expression_1, time_points_1, list_of_DRGs_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, coefficients_1, adjacency_matrix_of_gene_regulatory_network_1, network_graph_1, graph_statistics_1, node_statistics_1, subject_name_1, gene_ID_type_1, indices_of_DRGs_1] = load_analysis(inputData{i,1}, inputData{i,2});
    
    statistics_of_current_analysis = {inputData{i,1}, inputData{i,2}, num2str(size(time_points_1,1)), num2str(size(list_of_DRGs_1,1)), num2str(size(list_of_gene_clusters_1,2))};    
    statistics_of_analyses = [statistics_of_analyses; statistics_of_current_analysis];
    
    DRGs{i} = list_of_DRGs_1;    
    GRMs{i} = list_of_gene_clusters_1;    
  end
  
  common_genes =  get_common_genes_across_conditions(DRGs);
  list_of_conditions = inputData(:,2)';
  
  cd(general_comparison_folder);
  
  output_common_genes_reports(common_genes, list_of_conditions);
  writetable(cell2table(statistics_of_analyses), 'Comparison.csv', 'WriteVariableNames', false);
  
  cd(Dynamics4GenomicBigData_HOME);
 
end

function output_common_genes_reports(common_genes, list_of_conditions)
  
  for size_of_output_combinations = 1:length(list_of_conditions)
    output_common_genes_report_for_one_combination(common_genes, size_of_output_combinations, list_of_conditions)
  end
  
end

function output_common_genes_report_for_one_combination(common_genes, size_of_output_combinations, list_of_conditions)
  output_matrix_by_combination_size = [repmat({'Condition'}, 1, size_of_output_combinations) {'#of DRGs in Common'} {'DRG file'}];
  
  common_drgs_folder = 'common_DRGs';
  mkdir(common_drgs_folder);
  cd(common_drgs_folder);
  
  for combination_of_conditions_index = 1:size(common_genes,1)
    combination_of_conditions = common_genes(combination_of_conditions_index,:);
    
    if(length(combination_of_conditions{1}) == size_of_output_combinations)
      conditions = combination_of_conditions{1};
      common_drgs = combination_of_conditions{2};
      common_drgs_file = [num2str(combination_of_conditions_index) '.csv'];
      output_matrix_by_combination_size = [output_matrix_by_combination_size; list_of_conditions(conditions) {num2str(length(common_drgs))} {common_drgs_file}];
      writetable(cell2table(common_drgs'), common_drgs_file, 'WriteVariableNames', false)
    end
  end
  
  cd('..');
  writetable(cell2table(output_matrix_by_combination_size), ['Common_genes_' num2str(size_of_output_combinations) '.csv'], 'WriteVariableNames', false);
end

function common_genes =  get_common_genes_across_conditions(DRGs)

  n = size(DRGs,2);

  % all_combinations_of_condition_indices{1} is the list of size-1 combinations of condition indices.
  all_combinations_of_condition_indices = [];
  for k=1:n
    all_combinations_of_condition_indices = [all_combinations_of_condition_indices; {combnk(1:n,k)}];
  end
  
  common_genes = [];
  
  % n is the number of conditions
  % k is the  combination size.
  % C(n,k) is the number of k-sized combinations from n n-sized set.
  for k=1:n
    k_sized_combinations = all_combinations_of_condition_indices{k};
    
    for combination_index = 1:size(k_sized_combinations,1)
      combination = k_sized_combinations(combination_index,:);      
      sets = [];
      
      for combination_member_index = 1:length(combination)
	combination_member = combination(combination_member_index);	
	genes_in_current_condition = unique(DRGs{combination_member})';	
	sets = [sets; {genes_in_current_condition}];	
      end
      
      intersection = get_intersection(sets);      
      common_genes = [common_genes; {combination} {intersection}];      
    end    
  end
end


% This function returns the intersection of all the rows in sets. That is to say, each row is a set. Each row is a cell array.

% Example input.

% sets = {[1, 2, 3]; [5 9 3 6 96 2]; [2 96 3]}
function intersection = get_intersection(sets)

  intersection = sets{1,:};
  
  for row_index = 2:size(sets,1)
    intersection = intersect(intersection, sets{row_index,:});
  end

end
