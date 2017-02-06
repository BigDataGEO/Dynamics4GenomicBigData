function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;
  
  input_file_name = input('Enter the input file name: ');
  
  parts = strsplit(input_file_name, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  GEO_number = parts{1};

  inputData = readtable(['Input/' input_file_name], 'Delimiter', ',', 'ReadVariableNames',false);
  inputData = table2cell(inputData);
  
  if(size(inputData,1) > 1)
  
    general_comparison_folder = [Dynamics4GenomicBigData_HOME, 'Output/' parts{1} '/Comparison/'];
    mkdir(general_comparison_folder);
    
    for i = 1:size(inputData,1) 
      
      % One to one analyses, that are perfomed only between conditions belonging to the same series.
      for j = 1:size(inputData,1)
	if (i~=j)
	
	  % Get the condition's analysis data.
	  [gene_expression_1, time_points_1, list_of_top_DRGs_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, coefficients_1, adjacency_matrix_of_gene_regulatory_network_1, network_graph_1, graph_statistics_1, node_statistics_1, subject_name_1, gene_ID_type_1, indices_of_top_DRGs_1, number_of_statistically_significant_DRGs, list_of_genes, gene_expression_sorted_by_F_value, list_of_probe_ids, indices_of_genes_sorted_by_F_value] = load_analysis(GEO_number, inputData{i,1});
      
	  % Get the second condition's analysis data.
	  [gene_expression_2, time_points_2, list_of_top_DRGs_2, list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2, coefficients_2, adjacency_matrix_of_gene_regulatory_network_2, network_graph_2, graph_statistics_2, node_statistics_2, subject_name_2, gene_ID_type_2, indices_of_top_DRGs_2, number_of_statistically_significant_DRGs_2, list_of_genes_2, gene_expression_sorted_by_F_value_2, list_of_probe_ids, indices_of_genes_sorted_by_F_value] = load_analysis(GEO_number, inputData{j,1});
 
%  	  one_to_one_comparison_folder = [general_comparison_folder, '/', GEO_number, '_vs_', GEO_number '/', inputData{i,1}, '_vs_', inputData{j,1}];
	  
	  if(strcmp(GEO_number, GEO_number))
	  
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    
	    % This section finds and outputs the list of common statistically significant DRGs.
	    
	    cd(general_comparison_folder);
	 
	    probes1 = strtrim(gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs,1));
	    probes2 = strtrim(gene_expression_sorted_by_F_value_2(1:number_of_statistically_significant_DRGs_2,1));
	    common_DRG_probes = intersect(probes1,probes2);
	    
	    [indices_of_common_DRGs, not_found] = find_strings_in_cell_array(probes1, common_DRG_probes);
	    
	    common_DRG_genes = gene_expression_sorted_by_F_value(indices_of_common_DRGs,2);
	    
	    writetable(cell2table([[{'Probe IDs'} {'Gene names'}]; [common_DRG_probes common_DRG_genes]]), 'Common_statistically_significant_DRGs.csv', 'WriteVariableNames', false);
	    
	    cd(Dynamics4GenomicBigData_HOME);
	    
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
%  	    one_to_one_comparison_folder = [general_comparison_folder, '/', GEO_number, '/', inputData{i,1}, '_vs_', inputData{j,1}];
	    one_to_one_comparison_folder = [general_comparison_folder, '/', inputData{i,1}, '_vs_', inputData{j,1}];
	    
	    mkdir(one_to_one_comparison_folder);
	    cd(one_to_one_comparison_folder);
	    
	    output_comparison_plots(subject_name_1, list_of_gene_clusters_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, zscore(gene_expression_2')', indices_of_top_DRGs_1, list_of_genes, list_of_probe_ids);

	    plot_cluster_matches(subject_name_1, gene_expression_by_cluster_1, list_of_cluster_means_1, time_points_1, subject_name_2, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
	    
	    cd(Dynamics4GenomicBigData_HOME);
	  end
	end
      end
    end
  end
end


function plot_cluster_matches(name_of_first_subject, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2)

    global Dynamics4GenomicBigData_HOME;

    x = list_of_cluster_means;
    y = list_of_cluster_means_2;
    
    z = corr(x', y', 'type', 'Spearman');
    
    output_folder = 'GRM_matching';
    mkdir(output_folder);
    cd(output_folder);
    
    for current_cluster = 1:size(z,1)
      % The following lines output the correlation of the current module (from the first subject) against all the modules from the second subject.
      [sorted_values sorted_indices] = sort(z(current_cluster,:), 'descend');
      correlations_of_current_cluster = [{['Modules from ' name_of_second_subject]} {['Spearman correlation with ' name_of_first_subject '''s M' num2str(current_cluster)]}];
      correlations_of_current_cluster = [correlations_of_current_cluster; [strcat('M', strread(num2str(sorted_indices),'%s')')' strread(num2str(sorted_values),'%s')]];
      writetable(cell2table(correlations_of_current_cluster), ['M' num2str(current_cluster) '.csv'], 'WriteVariableNames', false);
    
      % The following lines plot the current module against the second individual's modules with the lowest and highest correlation.
      module_with_highest_correlation = find(z(current_cluster,:) == max(z(current_cluster,:)), 1, 'first');
      module_with_lowest_correlation = find(z(current_cluster,:) == min(z(current_cluster,:)), 1, 'first');
      
      expression_of_first_subjects_current_cluster = gene_expression_by_cluster{current_cluster};      
      expression_of_second_subjects_1st_cluster = gene_expression_by_cluster_2{module_with_lowest_correlation};      
      expression_of_second_subjects_2nd_cluster = gene_expression_by_cluster_2{module_with_highest_correlation};
      
      expression_of_first_subjects_cluster = expression_of_first_subjects_current_cluster;
      name_of_first_subjects_cluster = ['M' num2str(current_cluster)];
      mean_of_first_subjects_cluster = list_of_cluster_means(current_cluster,:);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_lowest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_1st_cluster;
      mean_of_second_subjects_cluster = list_of_cluster_means_2(module_with_lowest_correlation,:);
      
      y_limits = [min([min(expression_of_first_subjects_current_cluster), min(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])-.05,max([max(expression_of_first_subjects_current_cluster), max(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])+.05];
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_limits);
      
      print('-dpsc2', '-append', [name_of_first_subjects_cluster '.ps']);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_highest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_2nd_cluster;
      mean_of_second_subjects_cluster = list_of_cluster_means_2(module_with_highest_correlation,:);
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_limits);
      
      print('-dpsc2', '-append', [name_of_first_subjects_cluster '.ps']);
      
      close all;
    end
    cd('..');
end

function plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_axis_limits)
      h8=figure('units', 'centimeters', 'position', [0, 0, 45, 20]);
      axisLabelFontSize = 9;
	    
      set(gcf, 'PaperPositionMode', 'manual');
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf, 'PaperPosition', [0 0 45 20]);
      set(gcf, 'PaperUnits', 'centimeters');
      set(gcf, 'PaperSize', [45 20]);
      
      % First subplot      
      subplot(1,2,1);
      
      plot(expression_of_first_subjects_cluster','-*b');
	  
      set(gca,'XTick', 1:size(time_points));
      set(gca,'XTickLabel', time_points);

      xlabel('Time', 'FontSize', axisLabelFontSize);
      ylabel('Expression', 'FontSize', axisLabelFontSize);
      
      hold on;

      plot(mean_of_first_subjects_cluster,'o-r','LineWidth',1.5);

      ylim(y_axis_limits);

      v = axis;

      handle=title(['Expression of ', name_of_first_subjects_cluster, ' from ', name_of_first_subject, '']);

      set(handle,'Position',[2.6 v(4)*1.03 0]);

      hold off;
      
      % Second subplot
      subplot(1,2,2);
      
      plot(expression_of_second_subjects_cluster','-*b');
	  
      set(gca,'XTick', 1:size(time_points_2));
      set(gca,'XTickLabel', time_points_2);

      xlabel('Time', 'FontSize', axisLabelFontSize);
      ylabel('Expression', 'FontSize', axisLabelFontSize);
      
      hold on;

      plot(mean_of_second_subjects_cluster,'o-r','LineWidth',1.5);

      ylim(y_axis_limits);

      v = axis;
      
      handle=title(['Expression of ', name_of_second_subjects_cluster, ' from ', name_of_second_subject, '']);

      set(handle,'Position',[2.35 v(4)*1.03 0]);

      hold off;
end

function output_comparison_plots(name_of_first_subject, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_2, indices_of_DRGs, list_of_genes, list_of_probe_ids)

      global Dynamics4GenomicBigData_HOME;
  
      output_folder = 'Gene_expression_of_both_subjects';
      mkdir(output_folder);
      cd(output_folder);
      
      alpha_threshold = 0.05;
      
      p_value_output_matrix = {'GRM number', 'p-value (Wilcoxon)', 'p-value (KS)', 'p-value (KW)', 'p-value (correlation with mean)', 'p-value (bootstrap)', 'p-value (permutation)', ['p-value < ' num2str(alpha_threshold)], '', 'Spearman correlation coefficient', 'Coefficient < 0.75?'};

      
      number_of_clusters = size(list_of_cluster_means,1);

      number_of_subplots = 2 * number_of_clusters;
      
      number_of_subplots_per_page = 2;
      number_of_columns_per_page = 2;    
      number_of_rows_per_page = number_of_subplots_per_page / number_of_columns_per_page;
      
      number_of_pages = ceil(number_of_subplots / number_of_subplots_per_page);
      number_of_plots_in_last_page = number_of_subplots_per_page;
      if mod(number_of_subplots, number_of_subplots_per_page)~=0
	number_of_plots_in_last_page = mod(number_of_subplots, number_of_subplots_per_page);
      end
	
      currentClusterIndex = 1;
      
      for b = 1:number_of_pages
	
	number_of_plots_in_current_page = number_of_subplots_per_page;
	
	if(b == number_of_pages) %i.e., if this is the last page
	  number_of_plots_in_current_page = number_of_plots_in_last_page;
	end

	h8=figure('units', 'centimeters', 'position', [0, 0, 45, 20]);
	axisLabelFontSize = 9;
	    
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', [0 0 45 20]);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [45 20]);

	gen = 1;
	    
	while gen <= number_of_plots_in_current_page
	  
	  expression_of_first_subject = gene_expression_by_cluster{currentClusterIndex};
	  expression_of_second_subject = gene_expression_2(indices_of_DRGs(list_of_gene_clusters{currentClusterIndex}),:);
	  
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
	  % This section calculates the Spearman correlation between each gene's expression in the first individual and the expression in the second individual.
	  % The information is output to a spreadsheet.
	  
	  genes_in_current_cluster = list_of_genes(indices_of_DRGs(list_of_gene_clusters{currentClusterIndex}));
	  probes_in_current_cluster = list_of_probe_ids(indices_of_DRGs(list_of_gene_clusters{currentClusterIndex}));
	  
	  correlations_of_current_cluster_header = [{['Probe ID']} {['Gene ID']} {['Spearman correlation between expression in ' name_of_first_subject ' and ' name_of_second_subject]}];
	  correlations_of_current_cluster = [];
	  scores = [];
	  
	  for index_of_genes_in_current_cluster=1:size(expression_of_first_subject,1)
	    exp1 = expression_of_first_subject(index_of_genes_in_current_cluster,:)';
	    exp2 = expression_of_second_subject(index_of_genes_in_current_cluster,:)';
	    
	    z = corr(exp1, exp2, 'type', 'Spearman');
	    
	    correlations_of_current_cluster = [correlations_of_current_cluster; [ strtrim({probes_in_current_cluster{index_of_genes_in_current_cluster}}) strtrim({genes_in_current_cluster{index_of_genes_in_current_cluster}}) {num2str(z)} ]];
	    
	    scores = [scores; {z}];
	  end
	  
	  [values indices] = sort(cell2mat(scores), 'descend');
	  correlations_of_current_cluster = correlations_of_current_cluster(indices,:);
	  correlations_of_current_cluster = [correlations_of_current_cluster_header; correlations_of_current_cluster];
	  
	  writetable(cell2table(correlations_of_current_cluster), ['M' num2str(currentClusterIndex) '.csv'], 'WriteVariableNames', false);

	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  

	  % Plot the expression of the first individual

	  subplot(number_of_rows_per_page,number_of_columns_per_page,gen);

	  plot(expression_of_first_subject','-*b');
	  
	  set(gca,'XTick', 1:size(time_points));
	  set(gca,'XTickLabel', time_points);

	  xlabel('Time', 'FontSize', axisLabelFontSize);
	  ylabel('Expression', 'FontSize', axisLabelFontSize);

	  hold on;

	  plot(list_of_cluster_means(currentClusterIndex,:),'o-r','LineWidth',1.5);

	  y_axis_limits = [min([min(expression_of_first_subject), min(expression_of_second_subject)])-.05,max([max(expression_of_first_subject), max(expression_of_second_subject)])+.05];
	  ylim(y_axis_limits);

	  v = axis;

	  handle=title(['Expression of genes in M', num2str(currentClusterIndex), ' from ', name_of_first_subject, '']);

	  set(handle,'Position',[3 v(4)*1.03 0]);

	  hold off;
		
		
	  % Plot the expression of the second individual
		
	  gen = gen + 1;

	  subplot(number_of_rows_per_page,number_of_columns_per_page,gen)

	  plot(expression_of_second_subject','-*b');
	  
	  hold on;

	  plot(mean(expression_of_second_subject,1),'o-r','LineWidth',1.5);

	  ylim(y_axis_limits);
	  
	  set(gca,'XTick', 1:size(time_points));
	  set(gca,'XTickLabel', time_points);

	  xlabel('Time', 'FontSize', axisLabelFontSize)

	  ylabel('Expression', 'FontSize', axisLabelFontSize)

	  hold on;

	  v = axis;

	  handle=title(['Expression of the same genes in ', name_of_second_subject, '']);

	  set(handle,'Position',[2.815 v(4)*1.03 0]);

	  hold off;
	  
	  [p_value_wilcoxon, p_value_ks, p_value_kruskalwallis, p_value_correlation_with_mean, p_value_bootstrap, p_value_permutation, spearman_correlation_coefficient] = compare_gene_expression(expression_of_first_subject, expression_of_second_subject, time_points);
	  
	  significant = p_value_wilcoxon < alpha_threshold | p_value_ks < alpha_threshold | p_value_kruskalwallis < alpha_threshold | p_value_correlation_with_mean < alpha_threshold | p_value_bootstrap < alpha_threshold | p_value_permutation < alpha_threshold;

	  
%  	  p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(spearman_correlation_coefficient), num2str(spearman_correlation_coefficient > 0.75)}];
	  
	  p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(p_value_wilcoxon), num2str(p_value_ks), num2str(p_value_kruskalwallis), num2str(p_value_correlation_with_mean), num2str(p_value_bootstrap), num2str(p_value_permutation), num2str(significant), '', num2str(spearman_correlation_coefficient), num2str(spearman_correlation_coefficient < 0.75)}];
	  
	  print(h8,'-dpdf', ['M' num2str(currentClusterIndex) '.pdf']);
	  
	  currentClusterIndex = currentClusterIndex + 1;

	  gen = gen + 1;

	end

	close all;
      end
      
      writetable(cell2table(p_value_output_matrix), 'values_control_vs_stimulus.csv', 'WriteVariableNames', false);
      
      cd('..');
end


function [p_value_wilcoxon, p_value_ks, p_value_kruskalwallis, p_value_correlation_with_mean, p_value_bootstrap, p_value_permutation, spearman_correlation_coefficient] = compare_gene_expression(expression_of_first_subject, expression_of_second_subject, time_points)
  
  first_data_set = mean(expression_of_first_subject, 1)';
  
  second_data_set = mean(expression_of_second_subject, 1)';
  
  p_value_wilcoxon = signrank(first_data_set, second_data_set);%requires column vectors
  
%    p_value = ranksum(first_data_set, second_data_set);%requires column vectors

  [h, p_value_ks] = kstest2(first_data_set',second_data_set');%requires row vectors

  [p_value_kruskalwallis,tbl,stats] = kruskalwallis([first_data_set second_data_set], [], 'off');%requires column vectors

  p_value_correlation_with_mean = compare_curves_correlation_with_mean(first_data_set, second_data_set);%requires column vectors
  
  p_value_bootstrap = compare_curves_bootstrap(first_data_set, second_data_set, time_points);%requires column vectors
  
  p_value_permutation = compare_curves_permutation(first_data_set, second_data_set, time_points);%requires column vectors
  
  spearman_correlation_coefficient = corr(first_data_set, second_data_set, 'type', 'Spearman');
  
end


function p_value = compare_curves_correlation_with_mean(first_data_set, second_data_set)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_CorrelationWithMean.R ' filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);

end

function p_value = compare_curves_bootstrap(first_data_set, second_data_set, time_points)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);
  
  time_points_filename = 'time_points.csv';

  csvwrite(time_points_filename, [time_points]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_Bootstrap.R ' filename ' ' time_points_filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);
  
  delete(time_points_filename);

end


function p_value = compare_curves_permutation(first_data_set, second_data_set, time_points)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);
  
  time_points_filename = 'time_points.csv';

  csvwrite(time_points_filename, [time_points]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_Permutation.R ' filename ' ' time_points_filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);
  
  delete(time_points_filename);

end


function [indices, not_found] = find_strings_in_cell_array(cell_array_to_search, cell_array_to_search_for)

  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
%      if(mod(i,1000)==0)
%        display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
%      end
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
%    display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
end
