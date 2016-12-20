set_paths_and_imports;

global Dynamics4GenomicBigData_HOME;

cd('Input');
s = dir('*.csv');
file_list = {s.name}';

GEO_number={};
condition = {};
samples = {};
time_points = {};
number_of_top_DRGs = {};

for i=1:length(file_list)
  [GEO_number{i}, condition{i}, samples{i}, time_points{i}, number_of_top_DRGs{i}] = read_input(file_list{i});
end

cd('..');

fprintf('\n');
display('Gene annotation is starting. This can take some time, please wait.');

for i=1:length(file_list)
  
  [gene_expression, time_points, list_of_top_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids, indices_of_genes_sorted_by_F_value, standardized_gene_expression] = load_analysis(GEO_number{i}, condition{i});
  
  
  path_to_condition_folder = ['Output/' GEO_number{i} '/Conditions/' condition{i} '/'];
  
  cd(path_to_condition_folder);
  
  fprintf('\n');
  display(['Annotating genes of condition ''' condition{i} ''' from GEO series ' GEO_number{i} '.']);
  
  [chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_top_DRGs, gene_ID_type);
  
  fprintf('\n');
  display(['Finished the annotation of genes of condition ''' condition{i} ''' from GEO series ' GEO_number{i} '.']);
  
  cd(Dynamics4GenomicBigData_HOME);
  
end

fprintf('\n');
display('The annotation is complete for all the subjects/conditions.');
