function write_condition_report(GEO_number, list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered, geoStruct)

  global Dynamics4GenomicBigData_HOME;
  
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
      
  mkdir(output_folder);
  cd(output_folder);
 
  options = struct('format','latex','showCode',false,'outputDir',output_folder,'stylesheet','document.xsl', 'imageFormat', 'png');

  save('paper.mat');
    
  [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, true);
    
  [list_of_genes_sorted_by_F_value, gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients, indices_of_DRGs, list_of_DRGs] = step_3(list_of_genes, gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs_considered, true);
    
  [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_DRGs, indices_of_DRGs, smooth_gene_expression, true);

  [coefficients, adjacency_matrix_of_gene_regulatory_network] = step_5(list_of_gene_clusters, time_points, indices_of_DRGs, fd_smooth_coefficients, true);

  [network_graph, graph_statistics, node_statistics] = step_6(adjacency_matrix_of_gene_regulatory_network, true);
  
  [chartReport, tableReport] = step_7(list_of_genes, list_of_gene_clusters, indices_of_DRGs, gene_ID_type);
  
  path_to_results_file = ['Results.mat'];
  
  save(path_to_results_file, 'gene_expression', 'time_points', 'list_of_DRGs', 'list_of_gene_clusters', 'gene_expression_by_cluster', 'list_of_cluster_means', 'coefficients', 'adjacency_matrix_of_gene_regulatory_network', 'network_graph', 'graph_statistics', 'node_statistics', 'subject_name', 'gene_ID_type', 'indices_of_DRGs', 'number_of_statistically_significant_DRGs', 'list_of_genes', 'list_of_genes_sorted_by_F_value', 'gene_expression_sorted_by_F_value');

  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part1.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part2.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part3.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Part4.tex'], output_folder);
  
  draft = fopen('paper.tex', 'wt');
  
  fid = fopen('Part1.tex');
  F = fread(fid, '*char')';
  fclose(fid);
    
  fprintf(draft,'%-50s\n',F);
  
  fprintf(draft,'%s', ['This manuscript applies the pipeline analysis proposed by Carey et al. (2016) in order to identify gene regulatory networks from the time course data available in data series ' GEO_number ' of the \textit{Gene Expression Omnibus (GEO)}. The analysis is focused on condition ``' strrep(condition, '_', '\_') ''''' of the original study associated to data series ' GEO_number '. ']);
  
  fprintf(draft,'%s\n', ['This condition has a total of ' num2str(length(time_points)) ' time points.']);
  
  if(isfield(geoStruct.Header.Series, 'title'))
    fprintf(draft,'%s', ['The original study associated to dataset ' GEO_number ' is titled: \textit{``' geoStruct.Header.Series.title '''''}. ']);
  end
  
  if(isfield(geoStruct.Header.Series, 'summary'))
    fprintf(draft,'%s\n\n', ['The authors summarize this study as follows.']);
    fprintf(draft,'%s\n\n', ['\textit{' geoStruct.Header.Series.summary '}']);
  end
  
  fprintf(draft,'%s', ['The pipeline analysis used in this article (Carey et al., 2016) is composed of a sequence of steps where the data is obtained, preprocessed and analyzed for the identification of dynamic response genes (\textit{i.e.}, genes that exhibit significant changes across time), the clustering of these and the discovery of a gene regulatory network between these clusters. ']);
  
  fprintf(draft,'%s\n\n', ['A broad description of the pipeline steps is provided in the following subsections. The results obtained from application of these to the time course data of condition ``' strrep(condition, '_', '\_') ''''' in series ' GEO_number ' can be found in Section~\ref{section:results}.']);
  
  fid = fopen('Part2.tex');
  F = fread(fid, '*char')';
  fclose(fid);
    
  fprintf(draft,'%-50s\t',F);
    
  text = ['data series ' GEO_number ' for the selected condition.   '];
  fprintf(draft, '%s', text);
    
  fid = fopen('Part3.tex');
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
  
  text = ['In order to annotate all the dynamic response genes, the full list of genes must be submitted to the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. This full list of genes can be found in supplementary file \href{Step_5/All\_DRGs.txt}{All\_DRGs.txt}.'];
  fprintf(draft, '%s', text);
  
  text = ['The list of genes belonging to module $M1$ can be found in supplementary file \href{Step_5/M1/Genes_in_M1.txt}{Genes\_in\_M1.txt}. This list can be used to annotate only the genes in $M1$ using the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. An analogous method should be used to annotate the genes in other gene response modules.'];
  fprintf(draft, '%s', text);
  
  if(not(isempty(tableReport)))
    text = ['The table and chart reports of $M1$ can be found in supplementary files \href{Step_5/M1/Table\_report\_of\_M1.csv}{Table\_report\_of\_M1.csv} and \href{Step_5/M1/Chart\_report\_of\_M1.csv}{Chart\_report\_of\_M1.csv}, respectively. '];
    fprintf(draft, '%s', text);
    
    text = ['The reports of the other modules can be found analogously.'];
    fprintf(draft, '%s\n\n', text);
  end
    
  fid = fopen('Part4.tex');
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
  delete('Part4.tex');
  
  if isunix()
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/compile.sh ' output_folder]);
  end
  
  cd(Dynamics4GenomicBigData_HOME);
end