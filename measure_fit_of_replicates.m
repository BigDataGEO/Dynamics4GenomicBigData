function measure_fit_of_replicates()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  macro_condition = parts{2};

  inputData = readtable(['Input/' input_file_name], 'Delimiter', ',', 'ReadVariableNames', false);
  inputData = table2cell(inputData);
  
  GEO_number = parts{1};
  
  geoStruct = get_geo_data(GEO_number);
  
  general_comparison_folder = [Dynamics4GenomicBigData_HOME, 'Results/' parts{1} '/Comparison_of_replicates/' macro_condition];
  mkdir(general_comparison_folder);
  
  conditions = inputData;
  
  gene_expression = {};
  time_points = {};
  list_of_top_DRGs = {};
  list_of_gene_clusters = {};
  gene_expression_by_cluster = {};
  list_of_cluster_means = {};
  coefficients = {};
  adjacency_matrix_of_gene_regulatory_network = {};
  network_graph = {};
  graph_statistics = {};
  node_statistics = {};
  subject_name = {};
  gene_ID_type = {};
  indices_of_top_DRGs = {};
  number_of_statistically_significant_DRGs = {};
  list_of_genes = {};
  gene_expression_sorted_by_F_value = {};
  list_of_probe_ids = {};
  
  list_of_statistically_significant_DRGs = {};
  
  for i = 1:size(conditions,1)  
    [gene_expression{i}, time_points{i}, list_of_top_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_top_DRGs{i}, number_of_statistically_significant_DRGs{i}, list_of_genes{i}, gene_expression_sorted_by_F_value{i}, list_of_probe_ids{i}, indices_of_genes_sorted_by_F_value{i}] = load_analysis(GEO_number, conditions{i});
    
    list_of_statistically_significant_DRGs{i} = gene_expression_sorted_by_F_value{i}(1:number_of_statistically_significant_DRGs{i},1:2);
    
    list_of_statistically_significant_DRGs{i} = cellfun(@num2str, list_of_statistically_significant_DRGs{i}, 'UniformOutput', false);
    
  end
  
  % The following two lines measure the noise between the replicates.
  % The results are stored in variable list_of_probes_genes_noise, which is a cell array where the first column is the probe ids, the second column is the gene names and the third column is the noise measurements across the replicates.
  % The order of the probe ids/gene names is the same as in the original GEO matrix.
  noise_per_gene = measure_noise_between_replicates(gene_expression);  
  list_of_probes_genes_noise = [list_of_probe_ids{1} list_of_genes{1} num2cell(noise_per_gene)];
  
  % The following lines perform the same function, but the resulting cell array lists probe ids/gene names sorted by noise, from lower to higher noise.
  [B,I] = sort(noise_per_gene);  
  list_of_probes_genes_noise_sorted_by_noise = [list_of_probe_ids{1}(I) list_of_genes{1}(I) num2cell(noise_per_gene(I))];

  
  cd(general_comparison_folder);
  
  % The two cell arrays constructed earlier are exported as .csv files.
  writetable(cell2table(list_of_probes_genes_noise), [parts{2} '_noise_per_gene_ALL_GENES.csv'], 'WriteVariableNames',false);  
  writetable(cell2table(list_of_probes_genes_noise_sorted_by_noise), [parts{2} '_noise_per_gene_ALL_GENES_SORTED_BY_NOISE.csv'], 'WriteVariableNames',false);
  
  
  
  
  for probe_id_index_seq = 1:min(10, length(I))
    probe_id_index = I(probe_id_index_seq);
  
    gene_expression_to_plot = [];
    for condition_index = 1:size(conditions,1)
      gene_expression_to_plot = [gene_expression_to_plot; gene_expression{condition_index}(probe_id_index,:)];
    end
    gene_expression_plot(gene_expression_to_plot, time_points{1}, ['Probe ' list_of_probe_ids{1}(probe_id_index)], 'Time', 'Genes', 'Expression');
    
    print(gcf,'-dpdf', [num2str(probe_id_index_seq) '_Probe_' list_of_probe_ids{1}{probe_id_index}]);
    close all;
  end
  
  
  for probe_id_index_seq = max(0,length(I)-10):length(I)
    probe_id_index = I(probe_id_index_seq);
  
    gene_expression_to_plot = [];
    for condition_index = 1:size(conditions,1)
      gene_expression_to_plot = [gene_expression_to_plot; gene_expression{condition_index}(probe_id_index,:)];
    end
    gene_expression_plot(gene_expression_to_plot, time_points{1}, ['Probe ' list_of_probe_ids{1}(probe_id_index)], 'Time', 'Genes', 'Expression');
  
    print(gcf,'-dpdf', [num2str(probe_id_index_seq) '_Probe_' list_of_probe_ids{1}{probe_id_index}]);
    close all;
  end
  
  
  cd(Dynamics4GenomicBigData_HOME);

end