function [GS, NS, SWI] = step_8(N, G, path, EAS, flder)

  for i = 1:N
    GS(i,:) = [graph_clustercoeff(sparse(G{i})),graph_diameter(G{i}),graph_meandist(G{i}),graph_density(G{i})];
  end

  for i = 1:N
    create_exel_file('Graph Statistics.xls',GS,i,[],path);
  end
  
  for i = 1:N
    NS{i} = [bridging_centrality(G{i}),closeness_centrality(sparse(G{i})),eccentricity_centrality(sparse(G{i}))'];
    SWI(i) = smallworldindex(G{i});
  end

  for i = 1:N
    create_exel_file('Node Statistics.xls',GS,i,[],path);
  end
  
  % % % %   Visualization

  adjacency_matrix = EAS{i};

  %  adjacency_matrix = [1 1 1;  0 1 1; 1 0 0]; % example

  number_of_clusters = size(adjacency_matrix,1);

  dependency_matrix = repmat({''}, [number_of_clusters 2]);

  module_name_preffix = 'M';
  withinFieldSeparator = ';';
  betweenFieldSeparator = ',';

  row_labels = repmat({''}, [number_of_clusters 1]);
  column_labels = {'' 'Out' 'In'};

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

  create_exel_file('Dependency_matrix.xls',matrix_to_save,i,[],path);
  disp(strcat('This is the matrix listing the <a href="',flder,'/Dependency_matrix.xls">dependencies between GRMs</a>.'));


  matrix_to_save = [[{''} row_labels']; [row_labels num2cell(adjacency_matrix)]];

  create_exel_file('Adjacency_matrix.xls',matrix_to_save,i,[],path);
  disp(strcat('This is the matrix listing the <a href="',flder,'/Adjacency_matrix.xls">dependencies between GRMs</a>.'));

  % View network in a plot
  tgfFile = fopen('Network.tgf','w');
  sifFile = fopen('Network.sif','w');

  i = 1;
  for u = 1:length(EAS{i})
    for v = 1:length(EAS{i})
      if(EAS{i}(u,v) ~= 0)
	fprintf(tgfFile,'"Cluster %1.1i"	"Cluster %1.1i"\n',u,v);
	fprintf(sifFile,'"Cluster %1.1i"	pp	"Cluster %1.1i"\n',u,v);
      end
    end
  end

  fclose(tgfFile);
  fclose(sifFile);

  disp(strcat('This is <a href="',flder,'/Network.sif">the network file</a> that can be imported into Cytoscape.'));
  disp(strcat('This is <a href="',flder,'/Network.sif">the network file</a> that can be imported into BioLayout Express 3D.'));

  if isunix()
    [x, y]=system(['Rscript ', path, 'plot_network.R', ' Network_matrix.csv']);
  end