function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  output_dir = parts{1};

  inputData = readtable(input_file_name, 'Delimiter', ',');
  inputData = table2cell(inputData);
  
  general_comparison_folder = [Dynamics4GenomicBigData_HOME, 'Results/Comparison/'];
  mkdir(general_comparison_folder);
  
  for i = 1:size(inputData,1) 
     
    % One to one analyses, that are perfomed only between conditions belonging to the same series.
    for j = 1:size(inputData,1)
      if (i~=j)
      
	% Get the condition's analysis data.
	[gene_expression_1, time_points_1, list_of_DRGs_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, coefficients_1, adjacency_matrix_of_gene_regulatory_network_1, network_graph_1, graph_statistics_1, node_statistics_1, subject_name_1, gene_ID_type_1, indices_of_DRGs_1, number_of_statistically_significant_DRGs, list_of_genes, list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value] = load_analysis(inputData{i,1}, inputData{i,2});
    
	% Get the second condition's analysis data.
	[gene_expression_2, time_points_2, list_of_DRGs_2, list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2, coefficients_2, adjacency_matrix_of_gene_regulatory_network_2, network_graph_2, graph_statistics_2, node_statistics_2, subject_name_2, gene_ID_type_2, indices_of_DRGs_2, list_of_genes_2, list_of_genes_sorted_by_F_value_2, gene_expression_sorted_by_F_value_2] = load_analysis(inputData{j,1}, inputData{j,2});
	
	one_to_one_comparison_folder = [general_comparison_folder, '/', inputData{i,1}, '_vs_', inputData{j,1} '/', inputData{i,2}, '_vs_', inputData{j,2}];
	
	if(strcmp(inputData{i,1}, inputData{j,1}))
	  one_to_one_comparison_folder = [general_comparison_folder, '/', inputData{i,1}, '/', inputData{i,2}, '_vs_', inputData{j,2}];
	  
	  mkdir(one_to_one_comparison_folder);
	  cd(one_to_one_comparison_folder);
	  
	  output_comparison_plots(subject_name_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, zscore(gene_expression_2')', indices_of_DRGs_1, list_of_genes);

	  plot_cluster_matches(subject_name_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
	  
	  cd(Dynamics4GenomicBigData_HOME);
	end
      end
    end
  end
end
