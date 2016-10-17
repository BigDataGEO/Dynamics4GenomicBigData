function kompare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;

  input = readtable('input.csv', 'Delimiter', ',');
  input = table2cell(input);

  for i = 1:size(input,1)  

    for j = 1:size(input,1)
      if (i~=j)
	% Do the comparison between i and j.
	[gene_expression_1, time_points_1, list_of_DRGs_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, coefficients_1, adjacency_matrix_of_gene_regulatory_network_1, network_graph_1, graph_statistics_1, node_statistics_1, subject_name_1, gene_ID_type_1, indices_of_DRGs_1] = load_analysis(input{i,1}, input{i,2});
	
	[gene_expression_2, time_points_2, list_of_DRGs_2, list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2, coefficients_2, adjacency_matrix_of_gene_regulatory_network_2, network_graph_2, graph_statistics_2, node_statistics_2, subject_name_2, gene_ID_type_2, indices_of_DRGs_2] = load_analysis(input{j,1}, input{j,2});
	
	output_folder = [Dynamics4GenomicBigData_HOME, 'Results/Comparison/', input{i,1}, '_vs_', input{j,1} '/', input{i,2}, '_vs_', input{j,2}];
	
	if(input{i,1} == input{j,1})
	  output_folder = [Dynamics4GenomicBigData_HOME, 'Results/Comparison/', input{i,1}, '/', input{i,2}, '_vs_', input{j,2}];
	end
    
	mkdir(output_folder);
	cd(output_folder);
	
	
	output_comparison_plots(subject_name_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, zscore(gene_expression_2')', indices_of_DRGs_1);

	plot_cluster_matches(subject_name_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
	
	cd(Dynamics4GenomicBigData_HOME);
	
      end
    end
  end
end