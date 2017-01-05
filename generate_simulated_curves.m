function simulated_gene_expression = generate_simulated_curves(list_of_gene_clusters, gene_expression, number_of_non_DRGs, indices_of_top_DRGs)
  simulated_gene_expression = [];
  
  list_of_cluster_means = [];
  gene_expression_by_cluster = [];
  
  for i=1:length(list_of_gene_clusters)
    
    if(size(gene_expression(indices_of_top_DRGs(list_of_gene_clusters{i}),:),1) == 1)
      
      list_of_cluster_means = [list_of_cluster_means; mean(gene_expression(indices_of_top_DRGs(list_of_gene_clusters{i}),:),1)];    
      gene_expression_by_cluster = [gene_expression_by_cluster; {gene_expression(indices_of_top_DRGs(list_of_gene_clusters{i}),:)}];
    
    else
      list_of_cluster_means = [list_of_cluster_means; mean(gene_expression(indices_of_top_DRGs(list_of_gene_clusters{i}),:))];    
      gene_expression_by_cluster = [gene_expression_by_cluster; {gene_expression(indices_of_top_DRGs(list_of_gene_clusters{i}),:)}];
    end

  end
  
  simulated_gene_expression_DRGs = generate_simulated_curves_from_clusters(list_of_cluster_means, gene_expression_by_cluster); % i.e., the curves.
  
  mean_of_all_DRGs = mean(simulated_gene_expression_DRGs,1);
  
  simulated_gene_expression_non_DRGs = generate_simulated_curves_from_mean_curve(number_of_non_DRGs, mean_of_all_DRGs, 1);
  
  simulated_gene_expression = [simulated_gene_expression_DRGs; simulated_gene_expression_non_DRGs];
  
  % The following two lines simply shuffle the simulated curves.
  ordering = randperm(size(simulated_gene_expression,1));
  simulated_gene_expression = simulated_gene_expression(ordering, :);
  
end

