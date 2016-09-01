set_paths_and_imports;

[list_of_genes, raw_gene_expression, raw_time_points, number_of_top_DRGs_considered, gene_ID_type] = step_1('GSE47962');

[gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, false);

[fd_smooth_coefficients, smooth_gene_expression, indices_of_DRGs, list_of_DRGs, indices_of_genes_sorted_by_F_value] = step_3(time_points, smooth_gene_trajectories, gene_expression, list_of_genes, number_of_top_DRGs_considered, false);

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(smooth_gene_expression, indices_of_DRGs, time_points, number_of_top_DRGs_considered, gene_expression, indices_of_genes_sorted_by_F_value, list_of_DRGs, false);

step_5(list_of_gene_clusters, list_of_genes, indices_of_DRGs, gene_ID_type);

adjacency_matrix_of_gene_regulatory_network = step_7(list_of_gene_clusters, indices_of_DRGs, fd_smooth_coefficients, time_points, false);

step_8(adjacency_matrix_of_gene_regulatory_network);



