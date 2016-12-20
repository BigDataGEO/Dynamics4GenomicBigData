% OPTIONAL
% Further data manipulation that could be carried out after the pipeline steps.

% 1. Plotting the expression of a particular cluster of genes.
load Results.mat
cluster_index = 1;
gene_expression = zscore(gene_expression(indices_of_top_DRGs,:)')';
gene_expression_plot(gene_expression(list_of_gene_clusters{cluster_index},:), time_points, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['M3_R4.pdf']);
close all;


load Results.mat
gene_expression_1 = zscore(gene_expression(indices_of_top_DRGs,:)')';
time_points_1 = time_points;
list_of_gene_clusters_1 = list_of_gene_clusters;

load Results.mat
gene_expression_2 = zscore(gene_expression(indices_of_top_DRGs,:)')';
time_points_2 = time_points;
list_of_gene_clusters_2 = list_of_gene_clusters;

load Results.mat
gene_expression_3 = zscore(gene_expression(indices_of_top_DRGs,:)')';
time_points_3 = time_points;
list_of_gene_clusters_3 = list_of_gene_clusters;

cluster_index = 2;
m1 = intersect(list_of_gene_clusters_3{cluster_index}, intersect(list_of_gene_clusters_1{cluster_index}, list_of_gene_clusters_2{cluster_index}));

gene_expression_plot(gene_expression_1(m1,:), time_points_1, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['M1_R1_common.pdf']);

gene_expression_plot(gene_expression_2(m1,:), time_points_2, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['M1_R2_common.pdf']);

gene_expression_plot(gene_expression_3(m1,:), time_points_3, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['M1_R3_common.pdf']);

close all;


% 2. To plot the expression of list of genes read from a file.
filename = 'Genes_in_M1.txt';
cell_array_of_strings_to_find = table2cell(readtable(filename, 'Delimiter', ','));
indices = find_strings_in_cell_array(list_of_genes, cell_array_of_strings_to_find);
gene_expression_plot(gene_expression(indices,:), time_points, 'Gene expression', 'Time', 'Genes', 'Expression');
print(gcf,'-dpdf', ['Gene_expression.pdf']);