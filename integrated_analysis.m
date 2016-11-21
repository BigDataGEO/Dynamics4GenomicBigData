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
    
  statistics_of_analyses = {'Series', 'Condition', '# of time points', '# of DRGs', '# of Top DRGs for comparison', '# of GRMs'};
  
  DRGs = {};
  GRMs = {};

  for i = 1:size(inputData,1)  
    [gene_expression, time_points, list_of_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_DRGs, number_of_statistically_significant_DRGs, list_of_genes, list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value] = load_analysis(inputData{i,1}, inputData{i,2});
    
    statistics_of_current_analysis = {inputData{i,1}, inputData{i,2}, num2str(size(time_points,1)), num2str(number_of_statistically_significant_DRGs), num2str(size(list_of_DRGs,1)), num2str(size(list_of_gene_clusters,2))};    
    statistics_of_analyses = [statistics_of_analyses; statistics_of_current_analysis];
    
    DRGs{i} = list_of_DRGs;    
    GRMs{i} = list_of_gene_clusters;    
  end
  
  frequency_of_DRGs =  get_frequency_of_DRGs(DRGs);
  
  cd(general_comparison_folder);
  
  writetable(cell2table(frequency_of_DRGs), 'Frequency_of_DRGs.csv', 'WriteVariableNames', false);
  writetable(cell2table(statistics_of_analyses), 'Comparison.csv', 'WriteVariableNames', false);
  
  cd(Dynamics4GenomicBigData_HOME);
 
end

