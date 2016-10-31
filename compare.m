function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;

  fprintf('\n');
  GEO_number = input(['Please enter the accession number of your dataset enclosed in single quotes (e.g., ''GSE52428''): ']);
  
  try
    [geoStruct, list_of_genes, gene_ID_type] = get_geo_data(GEO_number);
  catch
    display(['Could not retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.']);
    return;
  end

  [raw_gene_expression, raw_time_points, name_of_first_subject, condition, number_of_top_DRGs_considered] = step_1(geoStruct);
    
  display(sprintf('\nYou have successfully entered the data for the first subject.'));
    
  prompt = '\nNow you will be required to enter the same information for the second subject (press enter to continue).';
  
  uselessVariable = input(prompt);
  
  [raw_gene_expression_2, raw_time_points_2, name_of_second_subject, condition_2, number_of_top_DRGs_considered_2] = step_1(geoStruct);
    
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
    
  mkdir(output_folder);
  cd(output_folder);
    
  % Steps 2, 3, and 4 of the pipeline are run for the first subject.
  [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, false);
  [list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_DRGs, list_of_DRGs] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs_considered, false);
  [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_DRGs, indices_of_DRGs, smooth_gene_expression, false);
    
  % Steps 2, 3, and 4 of the pipeline are run for the second subject.    
  [gene_expression_2, time_points_2, smooth_gene_trajectories_2] = step_2(raw_gene_expression_2, raw_time_points_2, false);
  [list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression_2, fd_smooth_coefficients_2, indices_of_DRGs_2, list_of_DRGs_2] = step_3(list_of_genes, gene_expression_2, time_points_2, smooth_gene_trajectories_2, number_of_top_DRGs_considered, false);  
  [list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2] = step_4(gene_expression_2, time_points_2, list_of_DRGs_2, indices_of_DRGs_2, smooth_gene_expression_2, false);
  
  
  output_comparison_plots(name_of_first_subject, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, zscore(gene_expression_2')', indices_of_DRGs);
    
  plot_cluster_matches(name_of_first_subject, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
    
  close all;
  cd(Dynamics4GenomicBigData_HOME);
end
