function [network_graph, graph_statistics, node_statistics] = step_7(adjacency_matrix_of_gene_regulatory_network, output)

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
    outputFolder = 'Step_7';
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
    create_exel_file(graphStatsFileName,graph_statistics,1,[],Dynamics4GenomicBigData_HOME);
    movefile(graphStatsFileName, outputFolder);
    
    nodeStatsFileName = 'Node_Statistics.xls';
    create_exel_file(nodeStatsFileName,node_statistics,1,[],Dynamics4GenomicBigData_HOME);
    movefile(nodeStatsFileName, outputFolder);

    matrix_to_save = [column_labels; [row_labels dependency_matrix]];
    depMatrixFilename = 'Dependency_matrix.xls';
    create_exel_file(depMatrixFilename,matrix_to_save,1,[],Dynamics4GenomicBigData_HOME);
    movefile(depMatrixFilename, outputFolder);
    
    matrix_to_save = [[{''} row_labels']; [row_labels num2cell(adjacency_matrix_of_gene_regulatory_network)]];
    adjacencyMatrixFilename='Adjacency_matrix.xls';
    create_exel_file(adjacencyMatrixFilename,matrix_to_save,1,[],Dynamics4GenomicBigData_HOME);
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
    matrix_of_files_descs = [matrix_of_files_descs; [{depMatrixFilename} {'Matrix of dependencies between the gene response modules (GRM) in the gene regulatory network (GRN).'}]];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{networkSIF} {'Gene regulatory network in .sif format for import into Cytoscape.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{networkTGF} {'Gene regulatory network in .tgf format for import into BioLayoutExpress3D.'}]];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Network_plot_MATLAB.pdf'} {'Plot of the gene regulatory network (GRN).'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Network_plot_R.pdf'} {'Plot of the gene regulatory network (GRN).'}]];
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
end