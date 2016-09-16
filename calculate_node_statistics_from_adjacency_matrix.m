function node_statistics = calculate_node_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network)

  binary_adjacency_matrix = (adjacency_matrix_of_gene_regulatory_network~=0);

  NS = [calculate_node_indegree(binary_adjacency_matrix), calculate_node_outdegree(binary_adjacency_matrix), bridging_centrality(binary_adjacency_matrix), closeness_centrality(sparse(binary_adjacency_matrix)), eccentricity_centrality(sparse(binary_adjacency_matrix))', betweenness(binary_adjacency_matrix), clusteringcoeff(binary_adjacency_matrix,1:size(binary_adjacency_matrix,1))];
  
  labels = {'In-degree', 'Out-degree', 'Bridging centrality' 'Closeness centrality' 'Eccentricity centrality', 'Betweenness', 'Clustering coefficient'};
  
  values = num2cell(NS);
  
  node_statistics = vertcat(labels, values);
  
  preffix = 'M';
  row_labels = [{''}];
  for row_index = 1:size(adjacency_matrix_of_gene_regulatory_network,1)
    row_labels = [row_labels; {['M' num2str(row_index)]}];
  end
  
  node_statistics = horzcat(row_labels, node_statistics);
  
end

function node_indegree = calculate_node_indegree(binary_adjacency_matrix)

  node_indegree = [];

  for row_index=1:size(binary_adjacency_matrix,1)
    node_indegree(row_index) = nnz(binary_adjacency_matrix(row_index,:));
  end
  node_indegree = node_indegree';
end

function node_outdegree = calculate_node_outdegree(binary_adjacency_matrix)

  node_outdegree = [];
  
  for column_index=1:size(binary_adjacency_matrix,2)
    node_outdegree(column_index) = nnz(binary_adjacency_matrix(:, column_index));
  end
  node_outdegree = node_outdegree';

end
