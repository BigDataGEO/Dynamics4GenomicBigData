set_paths_and_imports;

Condition_file = 'GSE59015_-_Wildtype_-_3000.csv';

cd('Input');
[GEO_number, condition, samples, time_points, number_of_top_DRGs] = read_input([Condition_file]);
cd('..');

% The following function call may take some time to complete.
[geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number);

[raw_gene_expression, raw_time_points] = step_1(geoStruct, samples, time_points);

[gene_expression, time_points, smooth_gene_trajectories, standardized_gene_expression] = step_2(raw_gene_expression, raw_time_points, true);

[gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_top_DRGs, list_of_top_DRGs, indices_of_genes_sorted_by_F_value] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs, list_of_probe_ids, standardized_gene_expression, true);

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(list_of_probe_ids, list_of_genes, standardized_gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, true);

[coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_top_DRGs, fd_smooth_coefficients, true);

[network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, true);

[chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type);
