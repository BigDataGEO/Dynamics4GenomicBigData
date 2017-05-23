function write_study_report(GEO_number)

  global Dynamics4GenomicBigData_HOME;
  
%    geoStruct = get_geo_data(GEO_number);
  
  GEO_number_folder_path = [Dynamics4GenomicBigData_HOME, 'Output/', GEO_number];
  conditions_folder_path = [GEO_number_folder_path, '/', 'Conditions'];
  output_folder_path = [GEO_number_folder_path];
  mkdir(output_folder_path);
  
  conditions = get_subdirs(conditions_folder_path);
  
  gene_expression = {};
  time_points = {};
  list_of_top_DRGs = {};
  list_of_grms = {};
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
  indices_of_genes_sorted_by_F_value = {};
  standardized_gene_expression = {};
  
  list_of_statistically_significant_DRGs = {};
  
  save('AArgh.mat')
  
  for i = 1:size(conditions,1)

    results{i} = load_analysis(GEO_number, conditions{i});
    
    list_of_statistically_significant_DRGs{i} = results{i}.step_3.gene_expression_sorted_by_F_value(1:results{i}.step_3.number_of_statistically_significant_DRGs,2:3);
    
    list_of_statistically_significant_DRGs{i} = cellfun(@num2str, table2cell(list_of_statistically_significant_DRGs{i}), 'UniformOutput', false);
    
  end
  
  [frequency_of_DRGs, common_probes] =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs);
  
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
  
  
  fprintf(draft,'%s\n', ['This manuscript applies the pipeline analysis proposed by Carey et al. (2016) to analyze the time course data available in data series ' GEO_number ' of the \textit{Gene Expression Omnibus (GEO)} in order to identify differentially expressed genes and the gene regulatory network these comprise. The analysis is focused on the ' num2str(length(conditions)) ' experimental conditions listed below. ']);
  
  fprintf(draft,'%s\n', ['\begin{enumerate}']);
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    fprintf(draft,'%s\n', ['\item \texttt{' strrep(condition, '_', '\_') '}']);
  end
  
  fprintf(draft,'%s\n', ['\end{enumerate}']);

  fprintf(draft,'%s', ['\par The conditions in this study have at least ' num2str(min(cellfun(@length,time_points))) ' time points. ']);
  
%    if(isfield(geoStruct.Header.Series, 'title'))
%      fprintf(draft,'%s', ['The original study associated to dataset ' GEO_number ' is titled: \textit{``' geoStruct.Header.Series.title '''''}. ']);
%    end
  
%    if(isfield(geoStruct.Header.Series, 'summary'))
%      fprintf(draft,'%s\n\n', ['The authors summarize this study as follows.']);
%      fprintf(draft,'%s\n\n', ['\textit{' geoStruct.Header.Series.summary '}']);
%    end
  
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
  
  fprintf(draft,'%s\n', ['\begin{tabular}{|c|c|c|c|c|} \hline']);
  
  fprintf(draft,'%s\n', ['Condition & \# of time points & \# of DRGs & \# of top DRGs clustered & \# of GRMs \\ \hline']);
  
  statistics_of_analyses = {'Series', 'Condition', '# of time points', '# of DRGs', '# of Top DRGs for comparison', '# of GRMs'};
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    
    fprintf(draft,'\n\t%s\n', ['\texttt{' strrep(condition, '_', '\_') '}' ' & ' num2str(length(results{condition_iter_index}.step_2.time_points)) ' & ' num2str(results{condition_iter_index}.step_3.number_of_statistically_significant_DRGs) ' & ' num2str(length(results{condition_iter_index}.step_3.indices_of_top_DRGs)) ' & ' num2str(length(results{condition_iter_index}.step_4.list_of_grms)) ' \\ \hline']);
    
    statistics_of_current_analysis = {GEO_number, condition, num2str(size(results{condition_iter_index}.step_2.time_points,1)), num2str(results{condition_iter_index}.step_3.number_of_statistically_significant_DRGs), num2str(size(results{condition_iter_index}.step_3.list_of_top_DRGs,1)), num2str(size(results{condition_iter_index}.step_4.list_of_grms,2))};    
    statistics_of_analyses = [statistics_of_analyses; statistics_of_current_analysis];

  end
  
  fprintf(draft,'%s\n', ['\end{tabular}']);
  
  fprintf(draft,'%s\n', ['\end{center}']);
  
  fprintf(draft,'%s\n', ['\caption{Result statistics from the ' num2str(length(conditions)) ' conditions analyzed in series ' GEO_number '. Full statistics in can be found in supplementary file \href{Summary.csv}{Summary.csv}.}']);
  
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
  
  fprintf(draft,'%s\n', ['\caption{Frequent DRGs across the conditions analyzed from series ' GEO_number '. The first column indicates gene name and the second column indicates the number of probes where the gene appears as a DRG. Full list in supplementary file \href{Frequency_of_DRGs.csv}{Frequency\_of\_DRGs.csv}.}']);
  
  fprintf(draft,'%s\n', ['\label{table:frequent_DRGs}']);
  
  fprintf(draft,'%s\n', ['\end{table}']);
  
  
  for condition_iter_index = 1:length(conditions)
    condition = conditions{condition_iter_index};
    write_condition_section(draft, GEO_number, condition, results{condition_iter_index}.step_2.gene_expression, results{condition_iter_index}.step_2.time_points, results{condition_iter_index}.step_3.list_of_top_DRGs, results{condition_iter_index}.step_4.list_of_grms, results{condition_iter_index}.step_5.coefficients, results{condition_iter_index}.step_5.adjacency_matrix_of_gene_regulatory_network, results{condition_iter_index}.step_6.network_graph, results{condition_iter_index}.step_6.graph_statistics, results{condition_iter_index}.step_6.node_statistics, results{condition_iter_index}.gene_ID_type, results{condition_iter_index}.step_3.indices_of_top_DRGs, results{condition_iter_index}.step_3.number_of_statistically_significant_DRGs, results{condition_iter_index}.step_3.gene_expression_sorted_by_F_value);
  end

  fid = fopen('Part3.tex');
  F = fread(fid, '*char')';
  fclose(fid);    
  fprintf(draft,'%-50s\n', F);
  
  writetable(cell2table(statistics_of_analyses), 'Summary.csv', 'WriteVariableNames', false);
  writetable(cell2table(frequency_of_DRGs), 'Frequency_of_DRGs.csv', 'WriteVariableNames', false);
  
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/bibliography.bib'], output_folder);
  copyfile([Dynamics4GenomicBigData_HOME, '/latex/Study/plos2015.bst'], output_folder);
  
%    if isunix()
%      % The following line compiles the .tex file into a .pdf.
%      % Two output arguments (x and y) are used simply to prevent the output from being printed onscreen.
%      [x, y]=system([Dynamics4GenomicBigData_HOME 'latex/Study/compile.sh ' output_folder]);
%    end
  
  delete('Part1.tex');
  delete('Part2.tex');
  delete('Part3.tex');
  
  cd(Dynamics4GenomicBigData_HOME);

end

function write_condition_section(draft, GEO_number, condition, gene_expression, time_points, list_of_top_DRGs, list_of_grms, coefficients, adjacency_matrix_of_gene_regulatory_network, network_graph, graph_statistics, node_statistics, gene_ID_type, indices_of_top_DRGs, number_of_statistically_significant_DRGs, gene_expression_sorted_by_F_value)

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
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_2/Paper_01.png}']);
  fprintf(draft,'%s\n', ['\caption{All genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:allgenes_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_3/Smooth_expression_curves.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Smooth expression of all genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:smoothexp_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_3/Smooth_expression_of_DRGs.png}']);
  fprintf(draft,'%s\n', ['\caption{Dynamic response genes in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:drgs_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_3/Smooth_expression_of_top_DRGs.png}']);
  fprintf(draft,'%s\n', ['\caption{Top ranking DRGs in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:top_drgs_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);

  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_4/GRMs_1.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Expression of the gene response modules in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:grms_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);

  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_4/GRMs.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Mean curves of gene response modules grouped by cluster size in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{figure:grmstype_' condition '}}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  fprintf(draft,'%s\n', ['\begin{figure}']);
  fprintf(draft,'%s\n', ['\centering']);
  fprintf(draft,'%s\n', ['\includegraphics[width=\textwidth]{' Dynamics4GenomicBigData_HOME '/Output/' GEO_number '/Conditions/' condition '/Step_6/Network_plot_MATLAB.pdf}']);
  fprintf(draft,'%s\n', ['\caption{Gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.}']);
  fprintf(draft,'%s\n', ['\label{fig:generegnet_' condition '}']);
  fprintf(draft,'%s\n', ['\end{figure}']);
  
  % Graph statistics of the GRN
%    text = ['\begin{center} \begin{table} \centering \begin{tabular}{ | l | l | p{5cm} |} \hline Metric & Value  \\ \hline'];
%    fprintf(draft, '%s\n', text);
%  
%    text = [graph_statistics{1,1} ' & ' num2str(graph_statistics{1,2}) '  \\ \hline'];
%    fprintf(draft, '%s\n', text);
%      
%    text = [graph_statistics{2,1} ' & ' num2str(graph_statistics{2,2}) '  \\ \hline'];
%    fprintf(draft, '%s\n', text);
%      
%    text = [graph_statistics{3,1} ' & ' num2str(graph_statistics{3,2}) '  \\ \hline'];
%    fprintf(draft, '%s\n', text);
%      
%    text = [graph_statistics{4,1} ' & ' num2str(graph_statistics{4,2}) '  \\ \hline'];
%    fprintf(draft, '%s\n', text);
%      
%    text = ['\end{tabular} \caption{Graph metrics of the gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.} \label{table:graphstats_' condition '} \end{table} \end{center}'];
%    fprintf(draft, '%s\n\n', text);
%    
%    % Node statistics of the GRN
%    text = ['\begin{center} \begin{table} \centering \begin{tabular}{ | l | c | c |} \hline Metric & Top ranking module & Bottom ranking module \\ \hline'];
%    fprintf(draft, '%s\n', text);
%      
%    for metric_index = 2:size(node_statistics,2)
%      
%      [max_value, index_of_max] = max([node_statistics{2:size(node_statistics,1),metric_index}]);
%      [min_value, index_of_min] = min([node_statistics{2:size(node_statistics,1),metric_index}]);
%        
%      text = [node_statistics{1,metric_index} ' & ' node_statistics{index_of_max+1,1} ' & ' node_statistics{index_of_min+1,1} '  \\ \hline'];
%      fprintf(draft, '%s\n', text);
%      
%    end
%    
%    text = ['\end{tabular} \caption{Node metrics of the gene regulatory network in condition \texttt{' strrep(condition, '_', '\_') '}.} \label{table:nodestats_' condition '} \end{table} \end{center}'];
%    fprintf(draft, '%s\n\n', text);
  
  % Gene annotation
  text = ['In order to annotate all the dynamic response genes in condition \texttt{' strrep(condition, '_', '\_') '}, the full list of genes must be submitted to the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. This full list of genes can be found in supplementary file \href{Conditions/' condition '/Step_7/All\_DRGs.txt}{All\_DRGs.txt}.'];
  fprintf(draft, '%s', text);
  
  text = ['The list of genes belonging to module $M1$ from condition \texttt{' strrep(condition, '_', '\_') '} can be found in supplementary file \href{Conditions/' condition '/Step_7/M1/Genes_in_M1.txt}{Genes\_in\_M1.txt}. This list can be used to annotate only the genes in $M1$ using the \textit{DAVID} \href{https://david.ncifcrf.gov}{website}. An analogous method should be used to annotate the genes in other gene response modules.'];
  fprintf(draft, '%s', text);
  
end

% This function finds the intersection of probe ids in the list of DRGs provided.

% This intersection of (DRG) probe ids is returned in variable common_probes.

% The list of gene names associated with each one of the common probes is returned in frequency_of_DRGs along with the number of (common) probes where the gene appears (which could be more than one).

% Example: Suppose the input is composed of the following two matrices.

% | Probe | Gene |		| Probe | Gene |
%     A      G1  		   X        G25
%     B      G52 		   B       G52
%     C      G41   		   D       G12
%     H      G15                   C       G41
%     F      G52                   F       G52

% The common probes are B, C, and F.

% The common DRGs are G52 (frequency 1), G41 (frequency 1), and G52 (frequency 2).

% Input:

% list_of_statistically_significant_DRGs is a cell array. Each element is a cell array of size Mx2 where the first column is the (DRG) probe ids and the second column the corresponding gene name, listing the DRGs of a subject/condition. M is the number of (DRG) probes.

% Output:

% frequency_of_DRGs is a Nx2 cell array, where the first column is the gene names (as strings) and the second column is the frequency, as numbers.
% common_probes is a Mx1 cell array with the probe ids of that are DRGs across all the subject/conditions.

function [frequency_of_DRGs, common_probes] =  get_frequency_of_DRGs(list_of_statistically_significant_DRGs)
  
  intersection_of_probes = list_of_statistically_significant_DRGs{1}(:,1);
  
  for k=1:length(list_of_statistically_significant_DRGs)
    intersection_of_probes = intersect(intersection_of_probes, list_of_statistically_significant_DRGs{k}(:,1));
  end
  
  A = list_of_statistically_significant_DRGs{1}(:,1);
  
  B = intersection_of_probes;
  
  common_probes = intersection_of_probes;
  
  intersection_of_probes_and_genes = list_of_statistically_significant_DRGs{1}(find(ismember(A,B)),:);
  
  frequency_of_DRGs = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,2));
  
  frequency_of_probes = get_frequency_of_each_array_element(intersection_of_probes_and_genes(:,1));
  
end

% Receives a cell array of strings of size Nx1 and returns an Nx2 cell array where the first column is the elements in the input array and the second column is the frequency of each element.

%  the_array = [{'a'}; {'s'}; {'a'}; {'e'}; {'q'}];

% Returns 
%  frequency_per_element = 
%  
%      'a'    [2]
%      'e'    [1]
%      'q'    [1]
%      's'    [1]

function frequency_per_element = get_frequency_of_each_array_element(the_array)

  [a b c] = unique(the_array);
  d = hist(c,length(a));
  P = [a num2cell(d')];
  
  [B I] = sort(cell2mat(P(:,2)), 'descend');
  
  frequency_per_element = P(I,:);  
end

function subdirs = get_subdirs(folder_name)
  d = dir(folder_name);
  isub = [d(:).isdir];
  subdirs = {d(isub).name}';
  subdirs(ismember(subdirs,{'.','..'})) = [];
end
