
set_paths_and_imports;

GEO_number = 'GSE59015';

condition = 'D10';

[gene_expression, time_points, list_of_top_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids, indices_of_genes_sorted_by_F_value, standardized_gene_expression] = load_analysis(GEO_number, condition);

number_of_non_DRGs = size(gene_expression,1) - length(indices_of_top_DRGs);

simulated_gene_expression = generate_simulated_curves(list_of_gene_clusters, gene_expression, number_of_non_DRGs, indices_of_top_DRGs);

mkdir(['Output/' GEO_number '_simulated/' condition]);
cd(['Output/' GEO_number '_simulated/' condition]);

[gene_expression, time_points, smooth_gene_trajectories, standardized_gene_expression] = step_2(simulated_gene_expression, time_points, true);

[gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_top_DRGs, list_of_top_DRGs, indices_of_genes_sorted_by_F_value] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, length(indices_of_top_DRGs), list_of_probe_ids, standardized_gene_expression, true);

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(list_of_probe_ids, list_of_genes, standardized_gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, true);

[coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_top_DRGs, fd_smooth_coefficients, true);

[network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, true);

close all;

cd('../../..')
