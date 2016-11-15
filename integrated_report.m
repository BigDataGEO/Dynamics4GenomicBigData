function integrated_report()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  GEO_number = input('Enter the GEO number: ');
  
  GEO_number_folder_path = [Dynamics4GenomicBigData_HOME, 'Results/', GEO_number];
  output_folder_name = 'Integrated';
  output_folder_path = [GEO_number_folder_path, '/', 'Integrated'];
  mkdir(output_folder_path);
  
  conditions = get_subdirs(GEO_number_folder_path);  
  idx = find(strcmp(conditions, output_folder_name));
  conditions(idx) = [];
  
  gene_expression = {};
  time_points = {};
  list_of_DRGs = {};
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
  indices_of_DRGs = {};
  number_of_statistically_significant_DRGs = {};
  list_of_genes = {};
  list_of_genes_sorted_by_F_value = {};
  gene_expression_sorted_by_F_value = {};

  % This loop simply reads all the results from all the conditions provided.
  for i = 1:size(conditions,1)  
    [gene_expression{i}, time_points{i}, list_of_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_DRGs{i}, number_of_statistically_significant_DRGs{i}, list_of_genes{i}, list_of_genes_sorted_by_F_value{i}, gene_expression_sorted_by_F_value{i}] = load_analysis(GEO_number, conditions{i}); 
  end
  
  cd(output_folder_path);
  
  write_study_report(GEO_number, conditions, gene_expression, time_points, list_of_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_DRGs, number_of_statistically_significant_DRGs, list_of_genes, list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value);
  
  cd(Dynamics4GenomicBigData_HOME);
 
end
