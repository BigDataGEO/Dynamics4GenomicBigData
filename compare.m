function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  output_dir = parts{1};

  inputData = readtable(input_file_name, 'Delimiter', ',');
  inputData = table2cell(inputData);
  
  for i = 1:size(inputData,1)
    condition{i} = inputData{i,2};
    GEO_number{i} = inputData{i,1};
    [gene_expression{i}, time_points{i}, list_of_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_DRGs{i}, number_of_statistically_significant_DRGs{i}] = load_analysis(GEO_number{i}, condition{i});
  end
  
  
  for i = 1:size(inputData,1)
    for j = 1:size(inputData,1)
      if(i ~= j)
	output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number{i},'/',condition{i});
	mkdir(output_folder);
	cd(output_folder);
	
	output_comparison_plots(subject_name{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, time_points{i}, subject_name{j}, zscore(gene_expression{j}')', indices_of_DRGs{i});
	
	plot_cluster_matches(subject_name{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, time_points{i}, subject_name{j}, gene_expression_by_cluster{j}, list_of_cluster_means{j}, time_points{j});
	
	cd(Dynamics4GenomicBigData_HOME)
      end
    end
  end 
end
