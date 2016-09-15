function node_statistics = calculate_node_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network)

  binary_adjacency_matrix = (adjacency_matrix_of_gene_regulatory_network~=0);

  NS = [bridging_centrality(binary_adjacency_matrix), closeness_centrality(sparse(binary_adjacency_matrix)), eccentricity_centrality(sparse(binary_adjacency_matrix))'];
  
  labels = {'Bridging centrality' 'Closeness centrality' 'Eccentricity centrality'};
  
  values = num2cell(NS);
  
  node_statistics = vertcat(labels, values);
  
  preffix = 'M';
  row_labels = [{''}];
  for row_index = 1:size(adjacency_matrix_of_gene_regulatory_network,1)
    row_labels = [row_labels; {['M' num2str(row_index)]}];
  end
  
  node_statistics = horzcat(row_labels, node_statistics);
  
end