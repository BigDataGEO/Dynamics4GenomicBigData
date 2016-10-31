function kompare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  output_dir = parts{1};

  inputData = readtable(input_file_name, 'Delimiter', ',');
  inputData = table2cell(inputData);
  
  for i = 1:size(inputData,1)
    condition{i} = inputData{i,2};
    GEO_number{i} = inputData{i,1};
    [gene_expression{i}, time_points{i}, list_of_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_DRGs{i}, number_of_statistically_significant_DRGs{i}] = load_analysis(GEO_number{i}, condition{i});
  end
  
  
  for i = 1:size(inputData,1)
    for j = 1:size(inputData,1)
      if(i ~= j)
	output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number{i},'/',condition{i});
	mkdir(output_folder);
	cd(output_folder);
	
	output_comparison_plots(subject_name{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, time_points{i}, subject_name{j}, zscore(gene_expression{j}')', indices_of_DRGs{i});
	
	plot_cluster_matches(subject_name{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, time_points{i}, subject_name{j}, gene_expression_by_cluster{j}, list_of_cluster_means{j}, time_points{j});
	
	cd(Dynamics4GenomicBigData_HOME)
      end
    end
  end 
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
