function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  GEO_number = parts{1};

  inputData = readtable(['Input/' input_file_name], 'Delimiter', ',', 'ReadVariableNames',false);
  inputData = table2cell(inputData);
  
  if(size(inputData,1) > 1)
  
    general_comparison_folder = [Dynamics4GenomicBigData_HOME, 'Results/' parts{1} '/Comparison/'];
    mkdir(general_comparison_folder);
    
    for i = 1:size(inputData,1) 
      
      % One to one analyses, that are perfomed only between conditions belonging to the same series.
      for j = 1:size(inputData,1)
	if (i~=j)
	
	  % Get the condition's analysis data.
	  [gene_expression_1, time_points_1, list_of_top_DRGs_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, coefficients_1, adjacency_matrix_of_gene_regulatory_network_1, network_graph_1, graph_statistics_1, node_statistics_1, subject_name_1, gene_ID_type_1, indices_of_top_DRGs_1, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids, indices_of_genes_sorted_by_F_value] = load_analysis(GEO_number, inputData{i,1});
      
	  % Get the second condition's analysis data.
	  [gene_expression_2, time_points_2, list_of_top_DRGs_2, list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2, coefficients_2, adjacency_matrix_of_gene_regulatory_network_2, network_graph_2, graph_statistics_2, node_statistics_2, subject_name_2, gene_ID_type_2, indices_of_top_DRGs_2, number_of_statistically_significant_DRGs_2, list_of_genes_2, gene_expression_sorted_by_F_value_2, list_of_probe_ids, indices_of_genes_sorted_by_F_value] = load_analysis(GEO_number, inputData{j,1});
 
%  	  one_to_one_comparison_folder = [general_comparison_folder, '/', GEO_number, '_vs_', GEO_number '/', inputData{i,1}, '_vs_', inputData{j,1}];
	  
	  if(strcmp(GEO_number, GEO_number))
	  
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    
	    % This section finds and outputs the list of common statistically significant DRGs.
	    
	    cd(general_comparison_folder);
	 
	    probes1 = strtrim(gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,1));
	    probes2 = strtrim(gene_expression_sorted_by_F_value_2(1:number_of_statistically_significant_DRGs_2,1));
	    common_DRG_probes = intersect(probes1,probes2);
	    
	    [indices_of_common_DRGs, not_found] = find_strings_in_cell_array(probes1, common_DRG_probes);
	    
	    common_DRG_genes = gene_expression_sorted_by_F_value(indices_of_common_DRGs,2);
	    
	    create_exel_file('Common_statistically_significant_DRGs.xls',[[{'Probe IDs'} {'Gene names'}]; [common_DRG_probes common_DRG_genes]],1,[],Dynamics4GenomicBigData_HOME);
	    cd(Dynamics4GenomicBigData_HOME);
	    
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
%  	    one_to_one_comparison_folder = [general_comparison_folder, '/', GEO_number, '/', inputData{i,1}, '_vs_', inputData{j,1}];
	    one_to_one_comparison_folder = [general_comparison_folder, '/', inputData{i,1}, '_vs_', inputData{j,1}];
	    
	    mkdir(one_to_one_comparison_folder);
	    cd(one_to_one_comparison_folder);
	    
	    output_comparison_plots(subject_name_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, zscore(gene_expression_2')', indices_of_top_DRGs_1, list_of_genes);

	    plot_cluster_matches(subject_name_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
	    
	    cd(Dynamics4GenomicBigData_HOME);
	  end
	end
      end
    end
  end
end
