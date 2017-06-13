set_paths_and_imports;

Condition_file = 'GSE19392_-_HBEs_infected_with_delNS1_post_trypsin_delNS1_1_-_3000.csv';

cd('Input');
[GEO_number, condition, samples, preprocessed_time_points, number_of_top_DRGs_considered] = read_input([Condition_file]);
cd('..');

[gene_expression, standardized_gene_expression, time_points, gene_ID_type, smooth_gene_trajectories] = step_2(GEO_number, samples, preprocessed_time_points, false);

[gene_expression_sorted_by_F_value, standardized_gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients] = step_3(gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs_considered, false);

list_of_grms = step_4(standardized_gene_expression_sorted_by_F_value, time_points, number_of_top_DRGs_considered, false);

[coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_grms, time_points, fd_smooth_coefficients, false);

[network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, false);

[chartReport, tableReport] = step_7(list_of_grms, gene_ID_type);
