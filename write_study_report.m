function write_study_report(GEO_number)

  global Dynamics4GenomicBigData_HOME;
  
  geoStruct = get_geo_data(GEO_number);
  
  GEO_number_folder_path = [Dynamics4GenomicBigData_HOME, 'Results/', GEO_number];
  conditions_folder_path = [GEO_number_folder_path, '/', 'Conditions'];
  output_folder_path = [GEO_number_folder_path];
  mkdir(output_folder_path);
  
  conditions = get_subdirs(conditions_folder_path);
  
  gene_expression = {};
  time_points = {};
  list_of_top_DRGs = {};
  list_of_gene_clusters = {};
  gene_expression_by_cluster = {};
  list_of_cluster_means = {};
  coefficients = {};
  adjacency_matrix_of_gene_regulatory_network = {};
  network_graph = {};
  graph_statistics = {};
  node_statistics = {};
  subject_name = {};
  gene_ID_type = {};
  indices_of_top_DRGs = {};
  number_of_statistically_significant_DRGs = {};
  list_of_genes = {};
  gene_expression_sorted_by_F_value = {};
  list_of_probe_ids = {};
  
  list_of_statistically_significant_DRGs = {};
  
  for i = 1:size(conditions,1)  
    [gene_expression{i}, time_points{i}, list_of_top_DRGs{i}, list_of_gene_clusters{i}, gene_expression_by_cluster{i}, list_of_cluster_means{i}, coefficients{i}, adjacency_matrix_of_gene_regulatory_network{i}, network_graph{i}, graph_statistics{i}, node_statistics{i}, subject_name{i}, gene_ID_type{i}, indices_of_top_DRGs{i}, number_of_statistically_significant_DRGs{i}, list_of_genes{i}, gene_expression_sorted_by_F_value{i}, list_of_probe_ids{i}] = load_analysis(GEO_number, conditions{i}); 
    
    list_of_statistically_significant_DRGs{i} = gene_expression_sorted_by_F_value{i}(1:number_of_statistically_significant_DRGs{i},2);
    
    list_of_statistically_significant_DRGs{i} = cellfun(@num2str, list_of_statistically_significant_DRGs{i}, 'UniformOutput', false);
    
  end

  frequency_of_DRGs =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs);
  
  cd(output_folder_path);
  
  output_folder = pwd;
  
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/Part1.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/Part2.tex'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/Part3.tex'], output_folder);
  
  draft = fopen('Report.tex', 'wt');
  
  fid = fopen('Part1.tex');
  F = fread(fid, '*char')';
  fclose(fid);    
  fprintf(draft,'%-50s\n', F);
  
  
  fprintf(draft,'%s\n', ['This manuscript applies the pipeline analysis proposed by Carey et al. (2016) in order to identify gene regulatory networks from the time course data available in data series ' GEO_number ' of the \textit{Gene Expression Omnibus (GEO)}. The analysis is focused on the ' num2str(length(conditions)) ' conditions listed below. ']);
  
  fprintf(draft,'%s\n', ['\begin{enumerate}']);
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    fprintf(draft,'%s\n', ['\item \texttt{' strrep(condition, '_', '\_') '}']);
  end
  
  fprintf(draft,'%s\n', ['\end{enumerate}']);

  fprintf(draft,'%s', ['\par This condition has a total of ' num2str(length(time_points)) ' time points.']);
  
  if(isfield(geoStruct.Header.Series, 'title'))
    fprintf(draft,'%s', ['The original study associated to dataset ' GEO_number ' is titled: \textit{``' geoStruct.Header.Series.title '''''}. ']);
  end
  
  if(isfield(geoStruct.Header.Series, 'summary'))
    fprintf(draft,'%s\n\n', ['The authors summarize this study as follows.']);
    fprintf(draft,'%s\n\n', ['\textit{' geoStruct.Header.Series.summary '}']);
  end
  
  fprintf(draft,'%s', ['The pipeline analysis used in this article (Carey et al., 2016) is composed of a sequence of steps where the data is obtained, preprocessed and analyzed for the identification of dynamic response genes (\textit{i.e.}, genes that exhibit significant changes across time), the clustering of these and the discovery of a gene regulatory network between these clusters. ']);
  
  fprintf(draft,'%s\n\n', ['A broad description of the pipeline steps is provided in the following subsections. The results obtained from application of these to the time course data of the conditions from series ' GEO_number ' listed earlier can be found in Section~\ref{section:results}.']);
  
  fid = fopen('Part2.tex');
  F = fread(fid, '*char')';
  fclose(fid);    
  fprintf(draft,'%-50s\n', F);
  
  % Table with summary of all analysed conditions.  
  fprintf(draft,'%s\n', ['\begin{table}']);
  
  fprintf(draft,'%s\n', ['\centering']);

  fprintf(draft,'%s\n', ['\begin{center}']);  
  
  fprintf(draft,'%s\n', ['\begin{tabular}{|c|c|c|c|} \hline']);
  
  fprintf(draft,'%s\n', ['Condition & \# of time points & \# of DRGs & \# of GRMs \\ \hline']);
  
  statistics_of_analyses = {'Series', 'Condition', '# of time points', '# of DRGs', '# of Top DRGs for comparison', '# of GRMs'};
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    
    fprintf(draft,'\n\t%s\n', ['\texttt{' strrep(condition, '_', '\_') '}' ' & ' num2str(length(time_points{condition_iter_index})) ' & ' num2str(number_of_statistically_significant_DRGs{condition_iter_index}) ' & ' num2str(length(list_of_gene_clusters{condition_iter_index})) ' \\ \hline']);
    
    statistics_of_current_analysis = {GEO_number, condition, num2str(size(time_points{condition_iter_index},1)), num2str(number_of_statistically_significant_DRGs{condition_iter_index}), num2str(size(list_of_top_DRGs{condition_iter_index},1)), num2str(size(list_of_gene_clusters{condition_iter_index},2))};    
    statistics_of_analyses = [statistics_of_analyses; statistics_of_current_analysis];

  end
  
  fprintf(draft,'%s\n', ['\end{tabular}']);
  
  fprintf(draft,'%s\n', ['\end{center}']);
  
  fprintf(draft,'%s\n', ['\caption{Result statistics from the conditions analyzed in series ' GEO_number '.}']);
  
  fprintf(draft,'%s\n', ['\label{table:summary}']);
  
  fprintf(draft,'%s\n', ['\end{table}']);
  
  
  % Table with top 25 frequent DRGs.
  top_25_frequent_DRGs = frequency_of_DRGs(1:min([25 size(frequency_of_DRGs, 1)]), :);
  
  fprintf(draft,'%s\n', ['\begin{table}']);
  
  fprintf(draft,'%s\n', ['\centering']);

  fprintf(draft,'%s\n', ['\begin{center}']);  
  
  fprintf(draft,'%s\n', ['\begin{tabular}{|c|c|} \hline']);
  
  fprintf(draft,'%s\n', ['Gene & Frequency \\ \hline']);
  
  for frequent_DRG_index = 1:size(top_25_frequent_DRGs,1)
    frequent_DRG = top_25_frequent_DRGs(frequent_DRG_index, :);
    
    fprintf(draft,'\n\t%s\n', ['\texttt{' strrep(frequent_DRG{1}, '_', '\_') '} & ' strrep(num2str(frequent_DRG{2}), '_', '\_') ' \\ \hline']);

  end
  
  fprintf(draft,'%s\n', ['\end{tabular}']);
  
  fprintf(draft,'%s\n', ['\end{center}']);
  
  fprintf(draft,'%s\n', ['\caption{Frequent DRGs across the conditions analyzed from series ' GEO_number '. The first column indicates gene name and the second column indicates the number of times the corresponding probe''s expression was measured as a statistically significant DRG. Full list in supplementary file \href{Frequency_of_DRGs.csv}{Frequency\_of\_DRGs.csv}.}']);
  
  fprintf(draft,'%s\n', ['\label{table:frequent_DRGs}']);
  
  fprintf(draft,'%s\n', ['\end{table}']);
  
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    write_condition_section(draft, GEO_number, condition, gene_expression{condition_iter_index}, time_points{condition_iter_index}, list_of_top_DRGs{condition_iter_index}, list_of_gene_clusters{condition_iter_index}, gene_expression_by_cluster{condition_iter_index}, list_of_cluster_means{condition_iter_index}, coefficients{condition_iter_index}, adjacency_matrix_of_gene_regulatory_network{condition_iter_index}, network_graph{condition_iter_index}, graph_statistics{condition_iter_index}, node_statistics{condition_iter_index}, subject_name, gene_ID_type{condition_iter_index}, indices_of_top_DRGs{condition_iter_index}, number_of_statistically_significant_DRGs{condition_iter_index}, gene_expression_sorted_by_F_value{condition_iter_index});
  end

  fid = fopen('Part3.tex');
  F = fread(fid, '*char')';
  fclose(fid);    
  fprintf(draft,'%-50s\n', F);
  
  writetable(cell2table(statistics_of_analyses), 'Summary.csv', 'WriteVariableNames', false);
  writetable(cell2table(frequency_of_DRGs), 'Frequency_of_DRGs.csv', 'WriteVariableNames', false);
  
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/bibliography.bib'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/plos2015.bst'], output_folder);
  
  if isunix()
    % The following line compiles the .tex file into a .pdf.
    % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
    [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/Study/compile.sh ' output_folder]);
  end
  
  delete('Part1.tex');
  delete('Part2.tex');
  delete('Part3.tex');
  
  cd(Dynamics4GenomicBigData_HOME);

end


function write_condition_section(draft, GEO_number, condition, gene_expression, time_points, list_of_top_DRGs, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, subject_name, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, gene_expression_sorted_by_F_value)

  global Dynamics4GenomicBigData_HOME;
  
  fprintf(draft,'%s\n', ['']);
  
  fprintf(draft,'%s\n', ['\subsection{Condition \texttt{' strrep(condition, '_', '\_') '}}']);
  
  fprintf(draft,'%s\n', ['\par Figure~\ref{fig:smoothexp_' condition '} shows the smooth expression curves.']);
  
  fprintf(draft,'%s\n', ['\par Figure~\ref{figure:drgs_' condition '} shows the expression of the top ranking dynamic response genes identified using the method described in Section~\ref{section:identification_of_drgs}.']);
  
  fprintf(draft,'%s\n', ['\par Figure~\ref{figure:grms_' condition '} shows the expression of gene response modules identified using the method described in Section~\ref{section:identification_of_grms}.']);
  
  fprintf(draft,'%s\n', ['\par Figure~\ref{figure:grmstype_' condition '} shows the GRMs'' mean expression curves grouped into four categories by cluster size. The four categories are single-gene modules (SGM) with only one gene in each cluster, small-size modules (SSM) that contain between 2-10 genes in each cluster, medium-size modules (MSM) that consist of 11-99 genes in each of the clusters and large-size modules (LSM) which contain over 100 genes in each cluster.']);
  
  fprintf(draft,'%s\n', ['\par Figure~\ref{fig:generegnet_' condition '} shows the gene regulatory network discovered using the method described in Section~\ref{section:identification_of_grn}.']);
  
  fprintf(draft,'%s\n', ['\par Graph theorists and network analysts have developed a number of metrics to characterize biological networks \cite{huber2007graphs, lee2004coexpression}. These metrics facilitate drug target identification and insight on potential strategies for treating various diseases.']);

  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_2/Paper_01.png}']);
  fprintf(draft,'%s\n', ['\caption{All genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:allgenes_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_3/Paper_03.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Smooth expression of all genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:smoothexp_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_3/Paper_04.png}']);
  fprintf(draft,'%s\n', ['\caption{Dynamic response genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:drgs_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_3/Paper_05.png}']);
  fprintf(draft,'%s\n', ['\caption{Top ranking DRGs in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:top_drgs_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);

  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_4/GRMs_1.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Expression of the gene response modules in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:grms_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);

  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_4/GRMs.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Mean curves of gene response modules grouped by cluster size in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:grmstype_' condition '}}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Results/' GEO_number '/Conditions/' condition '/Step_6/Network_plot_MATLAB.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:generegnet_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
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
    
  text = ['\end{tabular} \caption{Graph metrics of the gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.} \label{table:graphstats_' condition '} \end{table} \end{center}'];
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
  
  text = ['\end{tabular} \caption{Node metrics of the gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.} \label{table:nodestats_' condition '} \end{table} \end{center}'];
  fprintf(draft, '%s\n\n', text);
  
  % Gene annotation
  text = ['In order to annotate all the dynamic response genes in condition \texttt{' strrep(condition, '_', '\_') '}, the full list of genes must be submitted to the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. This full list of genes can be found in supplementary file \href{Conditions/' condition '/Step_7/All\_DRGs.txt}{All\_DRGs.txt}.'];
  fprintf(draft, '%s', text);
  
  text = ['The list of genes belonging to module $M1$ from condition \texttt{' strrep(condition, '_', '\_') '} can be found in supplementary file \href{Conditions/' condition '/Step_7/M1/Genes_in_M1.txt}{Genes\_in\_M1.txt}. This list can be used to annotate only the genes in $M1$ using the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. An analogous method should be used to annotate the genes in other gene response modules.'];
  fprintf(draft, '%s', text);
  
%    cd(Dynamics4GenomicBigData_HOME);
  
end
