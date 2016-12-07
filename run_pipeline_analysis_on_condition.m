function run_pipeline_analysis_on_condition(GEO_number, list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered, list_of_probe_ids, geoStruct)

  global Dynamics4GenomicBigData_HOME;
  
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/Conditions/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
    
  [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, true);

  [gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_top_DRGs, list_of_top_DRGs, indices_of_genes_sorted_by_F_value] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs_considered, list_of_probe_ids, true);

  [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(list_of_probe_ids, list_of_genes, gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, true);

  [coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_top_DRGs, fd_smooth_coefficients, true);

  [network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, true);

  [chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type);
  
  path_to_results_file = ['Results.mat'];
  
  save(path_to_results_file, 'gene_expression', 'time_points', 'list_of_top_DRGs', 'list_of_gene_clusters', 'gene_expression_by_cluster', 'list_of_cluster_means', 'coefficients', 'adjacency_matrix_of_gene_regulatory_network', 'network_graph', 'graph_statistics', 'node_statistics', 'subject_name', 'gene_ID_type', 'indices_of_top_DRGs', 'number_of_statistically_significant_DRGs', 'list_of_genes', 'gene_expression_sorted_by_F_value', 'list_of_probe_ids', 'indices_of_genes_sorted_by_F_value');
  
  close all;
  
  cd(Dynamics4GenomicBigData_HOME);
end