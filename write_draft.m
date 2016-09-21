function write_draft(GEO_number, list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered)

  global Dynamics4GenomicBigData_HOME;
  
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
 
  options = struct('format','latex','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl', 'imageFormat', 'png');

  save('paper.mat');
    
  [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, true);
    
  [list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, fd_smooth_coefficients] = step_3(list_of_genes, gene_expression, time_points, number_of_top_DRGs_considered, smooth_gene_trajectories, true);
    
  [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, number_of_top_DRGs_considered, true);
    
  step_5(list_of_genes, list_of_gene_clusters, indices_of_DRGs, gene_ID_type);

  [coefficients, adjacency_matrix_of_gene_regulatory_network] = step_7(list_of_gene_clusters, time_points, indices_of_DRGs, fd_smooth_coefficients, true);

  [network_graph, graph_statistics, node_statistics] = step_8(adjacency_matrix_of_gene_regulatory_network, true);

  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part1.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part2.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part3.tex'], output_folder);
    
  draft = fopen('paper.tex', 'wt');
    
  fid = fopen('Part1.tex');
  F = fread(fid, '*char')';
  fclose(fid);
    
  fprintf(draft,'%-50s\t',F);
    
  text = ['data series ' GEO_number ' for the selected condition.   '];
  fprintf(draft, '%s', text);
    
  fid = fopen('Part2.tex');
  F = fread(fid, '*char')';
  fclose(fid);
    
  fprintf(draft,'%-50s\t',F);
    
  % Graph statistics of the GRN
  text = ['\begin{center} \begin{table} \centering \begin{tabular}{ | l | l | p{5cm} |} \hline Metric & Value  \\ \hline'];
  fprintf(draft, '%s\n', text);

  text = [graph_statistics{1,1} ' & ' num2str(graph_statistics{1,2}) '  \\ \hline'];
  fprintf(draft, '%s\n', text);
    
  text = [graph_statistics{2,1} ' & ' num2str(graph_statistics{2,2}) '  \\ \hline'];
  fprintf(draft, '%s\n', text);
    
  text = [graph_statistics{3,1} ' & ' num2str(graph_statistics{3,2}) '  \\ \hline'];
  fprintf(draft, '%s\n', text);
    
  text = [graph_statistics{4,1} ' & ' num2str(graph_statistics{4,2}) '  \\ \hline'];
  fprintf(draft, '%s\n', text);
    
  text = ['\end{tabular} \caption{Graph metrics of the gene regulatory network.} \label{table:graphstats} \end{table} \end{center}'];
  fprintf(draft, '%s\n\n', text);
    
  % Node statistics of the GRN
  text = ['\begin{center} \begin{table} \centering \begin{tabular}{ | l | c | c |} \hline Metric & Top ranking module & Bottom ranking module \\ \hline'];
  fprintf(draft, '%s\n', text);
    
  for metric_index = 2:size(node_statistics,2)
    
    [max_value, index_of_max] = max([node_statistics{2:size(node_statistics,1),metric_index}]);
    [min_value, index_of_min] = min([node_statistics{2:size(node_statistics,1),metric_index}]);
      
    text = [node_statistics{1,metric_index} ' & ' node_statistics{index_of_max+1,1} ' & ' node_statistics{index_of_min+1,1} '  \\ \hline'];
    fprintf(draft, '%s\n', text);
    
  end
    
  text = ['\end{tabular} \caption{Node metrics of the gene regulatory network.} \label{table:nodestats} \end{table} \end{center}'];
  fprintf(draft, '%s\n\n', text);
    
  fid = fopen('Part3.tex');
  F = fread(fid, '*char')';
  fclose(fid);
    
  fprintf(draft,'%-50s\t',F);
    
  fclose(draft);
    
  close all;

  copyfile([Dynamics4GenomicBigData_HOME, '/latex/bibliography.bib'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/plos2015.bst'], output_folder);

  delete('paper.mat');
  delete('Part1.tex');
  delete('Part2.tex');
  delete('Part3.tex');
  
  if isunix()
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/compile.sh ' output_folder]);
  end
  
  cd(Dynamics4GenomicBigData_HOME);
end