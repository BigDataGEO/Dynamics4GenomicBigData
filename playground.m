
set_paths_and_imports;

GEO_number = 'GSE52428';

condition = 'H1N1_001';

[gene_expression, time_points, list_of_top_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids, indices_of_genes_sorted_by_F_value, standardized_gene_expression] = load_analysis(GEO_number, condition);

simulated_gene_expression = run_simulation(GEO_number, condition);



[gene_expression, time_points, smooth_gene_trajectories, standardized_gene_expression] = step_2(simulated_gene_expression, time_points, false);
  
[gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_top_DRGs, list_of_top_DRGs, indices_of_genes_sorted_by_F_value] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, 3000, list_of_probe_ids, standardized_gene_expression, false);

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(list_of_probe_ids, list_of_genes, standardized_gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, false);

[coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_top_DRGs, fd_smooth_coefficients, false);

[network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, false);

[chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type);