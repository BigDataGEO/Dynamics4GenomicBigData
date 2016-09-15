function graph_statistics = calculate_graph_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network)

  binary_adjacency_matrix = (adjacency_matrix_of_gene_regulatory_network~=0);

  GS = [graph_clustercoeff(sparse(binary_adjacency_matrix)), graph_diameter(binary_adjacency_matrix), graph_meandist(binary_adjacency_matrix), graph_density(binary_adjacency_matrix)]';
  
  labels = {'Clustering coefficient' 'Diameter' 'Mean distance' 'Density'}';
  
  values = num2cell(GS);
  
  graph_statistics = horzcat(labels, values);  
  
end