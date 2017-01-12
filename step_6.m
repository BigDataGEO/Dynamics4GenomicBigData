function [network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, output)

  global Dynamics4GenomicBigData_HOME;
  flder = pwd;

  number_of_clusters = size(adjacency_matrix_of_gene_regulatory_network,1);

  dependency_matrix = repmat({''}, [number_of_clusters 2]);

  module_name_preffix = 'M';
  withinFieldSeparator = ';';
  betweenFieldSeparator = ',';

  row_labels = repmat({''}, [number_of_clusters 1]);
  column_labels = {'' 'In (Influenced by) (Inward edges)' 'Out (Influences) (Outward edges)'};

  for u=1:number_of_clusters
    row_labels{u} = strcat(module_name_preffix, num2str(u));
    for v=1:number_of_clusters
      if adjacency_matrix_of_gene_regulatory_network(u,v) ~= 0
	if strcmp(dependency_matrix{u, 2}, '')
	  dependency_matrix{u, 2} = strcat(module_name_preffix, num2str(v));
	else
	  dependency_matrix{u, 2} = strcat(dependency_matrix{u, 2}, withinFieldSeparator, module_name_preffix, num2str(v));
	end
	if strcmp(dependency_matrix{v, 1}, '')
	  dependency_matrix{v, 1} = strcat(module_name_preffix, num2str(u));
	else
	  dependency_matrix{v, 1} = strcat(dependency_matrix{v, 1}, withinFieldSeparator, module_name_preffix, num2str(u));
	end
      end
    end
  end
  
  network_graph=digraph(adjacency_matrix_of_gene_regulatory_network, row_labels);
  graph_statistics = calculate_graph_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network);
  node_statistics = calculate_node_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network);
  
  if(output)
    outputFolder = 'Step_6';
    mkdir(outputFolder);
    
    newtwork_figure=figure('units', 'centimeters', 'position', [0, 0, 35, 40]);
    g_plot=plot(network_graph, 'Layout','force');    
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 60 50]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [60 50]);   

    print('Network_plot_MATLAB.pdf','-dpdf');
    movefile('Network_plot_MATLAB.pdf', outputFolder);
    
    graphStatsFileName = 'Graph_Statistics.xls';
    
    writetable(cell2table(graph_statistics), graphStatsFileName, 'WriteVariableNames', false);
    movefile(graphStatsFileName, outputFolder);
    
    nodeStatsFileName = 'Node_Statistics.xls';
    writetable(cell2table(node_statistics), nodeStatsFileName, 'WriteVariableNames', false);
    movefile(nodeStatsFileName, outputFolder);

    matrix_to_save = [column_labels; [row_labels dependency_matrix]];
    matrix_to_save = cell2table(matrix_to_save);
    
    depMatrixFilename = 'Dependency_matrix.csv';
    writetable(matrix_to_save, depMatrixFilename, 'WriteVariableNames', false);
    movefile(depMatrixFilename, outputFolder);
    
    matrix_to_save = [[{''} row_labels']; [row_labels num2cell(adjacency_matrix_of_gene_regulatory_network)]];
    matrix_to_save = cell2table(matrix_to_save);
    
    adjacencyMatrixFilename='Adjacency_matrix.csv';
    writetable(matrix_to_save, adjacencyMatrixFilename, 'WriteVariableNames', false);
    movefile(adjacencyMatrixFilename, outputFolder);

    networkTGF='Network.tgf';
    tgfFile = fopen(networkTGF,'w');
    
    networkSIF='Network.sif';
    sifFile = fopen(networkSIF,'w');

    for u = 1:length(adjacency_matrix_of_gene_regulatory_network)
      for v = 1:length(adjacency_matrix_of_gene_regulatory_network)
	if(adjacency_matrix_of_gene_regulatory_network(u,v) ~= 0)
	  fprintf(tgfFile,'"M%1.1i"	"M%1.1i"	%2.5f\n',u,v,adjacency_matrix_of_gene_regulatory_network(u,v));
	  fprintf(sifFile,'"M%1.1i"	pp	"M%1.1i"\n',u,v);
	end
      end
    end

    fclose(tgfFile);
    fclose(sifFile);
    
    movefile(networkTGF, outputFolder);
    movefile(networkSIF, outputFolder);

%      if isunix()
%        adjacencyMatrixCSVFilename = 'Network_matrix.csv';
%        csvwrite(adjacencyMatrixCSVFilename, adjacency_matrix_of_gene_regulatory_network');
%      
%        networkPlotFilename = 'Network_plot_R.pdf';
%        [x, y]=system(['Rscript ', Dynamics4GenomicBigData_HOME, 'plot_network.R', ' ', adjacencyMatrixCSVFilename, ' ', networkPlotFilename]);    
%        movefile(networkPlotFilename, outputFolder);
%        
%        delete(adjacencyMatrixCSVFilename);
%      end
    
    matrix_of_files_descs = [{'File name'} {'Description'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{adjacencyMatrixFilename} {'Adjacency matrix of the gene regulatory network (GRN) in Excel format.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{graphStatsFileName} {'Graph metrics of the gene regulatory network (GRN).'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{nodeStatsFileName} {'Node metrics of the gene regulatory network (GRN).'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{depMatrixFilename} {'Matrix of dependencies between the gene response modules (GRM) in the gene regulatory network (GRN).'}]];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{networkSIF} {'Gene regulatory network in .sif format for import into Cytoscape.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{networkTGF} {'Gene regulatory network in .tgf format for import into BioLayoutExpress3D.'}]];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Network_plot_MATLAB.pdf'} {'Plot of the gene regulatory network (GRN).'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Network_plot_R.pdf'} {'Plot of the gene regulatory network (GRN).'}]];
    
    writetable(cell2table(matrix_of_files_descs), 'List_and_description_of_output.csv', 'WriteVariableNames', false);

    movefile('List_and_description_of_output.csv', outputFolder);
  end
end

function graph_statistics = calculate_graph_statistics_from_adjacency_matrix(adjacency_matrix_of_gene_regulatory_network)

  binary_adjacency_matrix = (adjacency_matrix_of_gene_regulatory_network~=0);

  GS = [graph_clustercoeff(sparse(binary_adjacency_matrix)), graph_diameter(binary_adjacency_matrix), graph_meandist(binary_adjacency_matrix), graph_density(binary_adjacency_matrix)]';
  
  labels = {'Clustering coefficient' 'Diameter' 'Mean distance' 'Density'}';
  
  values = num2cell(GS);
  
  graph_statistics = horzcat(labels, values);  
  
end

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
