function step_8(EAS)

  global Dynamics4GenomicBigData_HOME;
  flder = pwd;
  
  outputFolder = 'Step_8';
  mkdir(outputFolder);
  
  G = (EAS~=0);

  GS(1,:) = [graph_clustercoeff(sparse(G)),graph_diameter(G),graph_meandist(G),graph_density(G)];
  NS{1} = [bridging_centrality(G),closeness_centrality(sparse(G)),eccentricity_centrality(sparse(G))'];
  SWI(1) = smallworldindex(G);

  graphStatsFileName = 'Graph_Statistics.xls';
  create_exel_file(graphStatsFileName,GS,1,[],Dynamics4GenomicBigData_HOME);
  movefile('Graph_Statistics.xls', outputFolder);

  adjacency_matrix = EAS;

  number_of_clusters = size(adjacency_matrix,1);

  dependency_matrix = repmat({''}, [number_of_clusters 2]);

  module_name_preffix = 'M';
  withinFieldSeparator = ';';
  betweenFieldSeparator = ',';

  row_labels = repmat({''}, [number_of_clusters 1]);
  column_labels = {'' 'Out (Influences) (Outward edges)' 'In (Influenced by) (Inward edges)'};

  for u=1:number_of_clusters
    row_labels{u} = strcat(module_name_preffix, num2str(u));
    for v=1:number_of_clusters
      if adjacency_matrix(u,v) ~= 0
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

  matrix_to_save = [column_labels; [row_labels dependency_matrix]];
  depMatrixFilename = 'Dependency_matrix.xls';
  create_exel_file(depMatrixFilename,matrix_to_save,1,[],Dynamics4GenomicBigData_HOME);
  movefile(depMatrixFilename, outputFolder);
  
  matrix_to_save = [[{''} row_labels']; [row_labels num2cell(adjacency_matrix)]];
  adjacencyMatrixFilename='Adjacency_matrix.xls';
  create_exel_file(adjacencyMatrixFilename,matrix_to_save,1,[],Dynamics4GenomicBigData_HOME);
  movefile(adjacencyMatrixFilename, outputFolder);

  networkTGF='Network.tgf';
  tgfFile = fopen(networkTGF,'w');
  
  networkSIF='Network.sif';
  sifFile = fopen(networkSIF,'w');

  for u = 1:length(EAS)
    for v = 1:length(EAS)
      if(EAS(u,v) ~= 0)
	fprintf(tgfFile,'"M%1.1i"	"M%1.1i"\n',v,u);
	fprintf(sifFile,'"M%1.1i"	pp	"M%1.1i"\n',v,u);
      end
    end
  end

  fclose(tgfFile);
  fclose(sifFile);
  
  movefile(networkTGF, outputFolder);
  movefile(networkSIF, outputFolder);

  adjacencyMatrixCSVFilename = 'Network_matrix.csv';
  csvwrite(adjacencyMatrixCSVFilename,EAS);
  
  if isunix()
    [x, y]=system(['Rscript ', Dynamics4GenomicBigData_HOME, 'plot_network.R', ' ', adjacencyMatrixCSVFilename]);    
    movefile('Network_plot.pdf', outputFolder);
  end
  
  movefile(adjacencyMatrixCSVFilename, outputFolder);
  
  matrix_of_files_descs = [{'File name'} {'Description'}];
  
  matrix_of_files_descs = [matrix_of_files_descs; [{adjacencyMatrixFilename} {'Adjacency matrix of the gene regulatory network (GRN) in Excel format.'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{adjacencyMatrixCSVFilename} {'Adjacency matrix of the gene regulatory network (GRN) in CSV format.'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{graphStatsFileName} {'Graph metrics of the gene regulatory network (GRN).'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{depMatrixFilename} {'Matrix of dependencies between the gene response modules (GRM) in the gene regulatory network (GRN).'}]];
  
  matrix_of_files_descs = [matrix_of_files_descs; [{networkSIF} {'Gene regulatory network in .sif format for import into Cytoscape.'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{networkTGF} {'Gene regulatory network in .tgf format for import into BioLayoutExpress3D.'}]];
  
  matrix_of_files_descs = [matrix_of_files_descs; [{'Network_plot.pdf'} {'Plot of the gene regulatory network (GRN).'}]];
  
  create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

  movefile('List_and_description_of_output.xls', outputFolder);
  
end