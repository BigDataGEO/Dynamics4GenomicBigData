function integrated_analysis()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  macro_condition = parts{2};

  inputData = readtable(['Input/' input_file_name], 'Delimiter', ',', 'ReadVariableNames', false);
  inputData = table2cell(inputData);
  
  GEO_number = parts{1};
  
  geoStruct = get_geo_data(GEO_number);
  
  GEO_number_folder_path = [Dynamics4GenomicBigData_HOME, 'Output/', GEO_number];
  conditions_folder_path = [GEO_number_folder_path, '/', 'Conditions'];
  output_folder_path = [GEO_number_folder_path];
  mkdir(output_folder_path);
  
  conditions = inputData;
  
  gene_expression = {};
  time_points = {};
  list_of_top_DRGs = {};
  list_of_gene_clusters = {};
  gene_expression_by_cluster = {};
  list_of_cluster_means = {};
  coefficients = {};
  adjacency_matrix_of_gene_regulatory_network = {};
  network_graph = {};
  graph_statistics = {};
  node_statistics = {};
  subject_name = {};
  gene_ID_type = {};
  indices_of_top_DRGs = {};
  number_of_statistically_significant_DRGs = {};
  list_of_genes = {};
  gene_expression_sorted_by_F_value = {};
  list_of_probe_ids = {};
  
  list_of_statistically_significant_DRGs = {};
  
  for i = 1:size(conditions,1)  
    [gene_expression{i}, time_points{i}, list_of_top_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_top_DRGs{i}, number_of_statistically_significant_DRGs{i}, list_of_genes{i}, gene_expression_sorted_by_F_value{i}, list_of_probe_ids{i}, indices_of_genes_sorted_by_F_value{i}] = load_analysis(GEO_number, conditions{i});
    
    list_of_statistically_significant_DRGs{i} = gene_expression_sorted_by_F_value{i}(1:number_of_statistically_significant_DRGs{i},1:2);
    
    list_of_statistically_significant_DRGs{i} = cellfun(@num2str, list_of_statistically_significant_DRGs{i}, 'UniformOutput', false);
    
  end
  
  [frequency_of_DRGs, common_probes] =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs);
  
  cd(output_folder_path);
  
  output_folder = pwd;
 
  
  statistics_of_analyses = {'Series', 'Condition', '# of time points', '# of DRGs', '# of Top DRGs for comparison', '# of GRMs'};
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    
    statistics_of_current_analysis = {GEO_number, condition, num2str(size(time_points{condition_iter_index},1)), num2str(number_of_statistically_significant_DRGs{condition_iter_index}), num2str(size(list_of_top_DRGs{condition_iter_index},1)), num2str(size(list_of_gene_clusters{condition_iter_index},2))};    
    statistics_of_analyses = [statistics_of_analyses; statistics_of_current_analysis];

  end
    
  writetable(cell2table(statistics_of_analyses), [macro_condition '_Summary.csv'], 'WriteVariableNames', false);
  writetable(cell2table(frequency_of_DRGs), [macro_condition '_Frequency_of_DRGs.csv'], 'WriteVariableNames', false);

  cd(Dynamics4GenomicBigData_HOME);
end
