set_paths_and_imports;

Condition_file = 'GSE41067_-_Replicate_1_-_3000.csv';

cd('Input');
[GEO_number, condition, samples, time_points, number_of_top_DRGs] = read_input([Condition_file]);
cd('..');

% The following function call may take some time to complete.
[geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number);

[raw_gene_expression, raw_time_points] = step_1(geoStruct, samples, time_points);

[gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, false);

[gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_top_DRGs, list_of_top_DRGs] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs, list_of_probe_ids, false);

[list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, false);

[coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_top_DRGs, fd_smooth_coefficients, false);

[network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, false);

[chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type);


% OPTIONAL
% Further data manipulation that could be carried out after the pipeline steps.

% 1. Plotting the expression of a particular cluster of genes.
cluster_index = 10;
gene_expression_plot(gene_expression(list_of_gene_clusters{cluster_index},:), time_points, 'Gene expression', 'Time', 'Genes', 'Expression');

% 2. To plot the expression of list of genes read from a file.
filename = 'Genes_in_M1.txt';
cell_array_of_strings_to_find = table2cell(readtable(filename, 'Delimiter', ','));
indices = find_strings_in_cell_array(list_of_genes, cell_array_of_strings_to_find);
gene_expression_plot(gene_expression(indices,:), time_points, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['Gene_expression.pdf']);
