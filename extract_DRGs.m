

GEO_number = 'GSE59015';
condition = 'D10';

[gene_expression, time_points, list_of_top_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids] = load_analysis(GEO_number, condition);

DRGs_probes = gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,1);
DRGs_gene_names = gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,2);


writetable(cell2table(DRGs_gene_names), '', 'WriteRowNames', false, 'WriteVariableNames', false, 'Delimiter', ',');