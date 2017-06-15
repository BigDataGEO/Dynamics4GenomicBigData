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
	  results_1 = load_analysis(GEO_number, inputData{i,1});
	  results_2 = load_analysis(GEO_number, inputData{j,1});
	  
	  % In the following lines we make sure that the two conditions have the same time points
	  % and complete where any are missing.	  
	  % This condition simply checks if the list of time points in both conditions are exactly
	  % the same. That is to say, it verifies that the time points are the same and in the same
	  % order. If the values are the same then they should also appear in the same order
	  % because the pipeline sorts them.
	  if ~isequal(results_1.step_2.time_points, results_2.step_2.time_points)
	    [results_1, results_2] = equalizeConditions(results_1, results_2);
	  end
	  
	  if(strcmp(GEO_number, GEO_number))
	    
	    % The following lines find and output the list of common statistically significant DRGs.
	    cd(general_comparison_folder);
	 
	    probes1 = strtrim(table2cell(results_1.step_3.gene_expression_sorted_by_F_value(1:results_1.step_3.number_of_statistically_significant_DRGs,2)));
	    probes2 = strtrim(table2cell(results_2.step_3.gene_expression_sorted_by_F_value(1:results_2.step_3.number_of_statistically_significant_DRGs,2)));
	    common_DRG_probes = intersect(probes1,probes2);
	    
	    [indices_of_common_DRGs, not_found] = find_strings_in_cell_array(probes1, common_DRG_probes);
	    
	    common_DRGs = results_1.step_3.gene_expression_sorted_by_F_value(indices_of_common_DRGs,2:3);
	    
	    writetable(common_DRGs, 'Common_statistically_significant_DRGs.csv', 'WriteVariableNames', true);
	    
	    cd(Dynamics4GenomicBigData_HOME);

	    % The following lines perform and output the actual comparisons.
	    one_to_one_comparison_folder = [general_comparison_folder, '/', inputData{i,1}, '_vs_', inputData{j,1}];
	    
	    mkdir(one_to_one_comparison_folder);
	    cd(one_to_one_comparison_folder);
	    
	    output_gene_correlations(results_1, results_2);
	    
	    output_comparison_plots(results_1.condition, results_1.step_4.list_of_grms, results_1.step_2.time_points, results_2.condition, cell2mat(table2cell(results_2.step_3.standardized_gene_expression(:,3:size(results_2.step_3.standardized_gene_expression, 2)))), results_1.step_3.indices_of_top_DRGs, results_1.list_of_genes, results_1.list_of_probe_ids);
	    
	    plot_cluster_matches(results_1.condition, results_1.step_4.list_of_grms, results_1.step_2.time_points, results_2.condition, results_2.step_4.list_of_grms, results_2.step_2.time_points);
	    
	    cd(Dynamics4GenomicBigData_HOME);
	  end
	end
      end
    end
  end
end

function [equalized_results_1, equalized_results_2] = equalizeConditions(results_1, results_2)

  old_results_1 = results_1;
  old_results_2 = results_2;
  
  new_time_points = sort(union(results_1.step_2.time_points, results_2.step_2.time_points));
  
  results_1.step_2.original_time_points = results_1.step_2.time_points;
  results_2.step_2.original_time_points = results_2.step_2.time_points;
  
  results_1.step_2.time_points = new_time_points;
  results_2.step_2.time_points = new_time_points;
  
  % The following are the fields where the interpolation will be required.
  % results_1.step_4.gene_expression_by_cluster
  % results_2.step_2.gene_expression  
  % results_1.step_4.list_of_cluster_means
  % results_2.step_4.list_of_cluster_means  
  
  % Interpolating in the clusters of the first condition
  for cluster_index=1:size(results_1.step_4.list_of_grms, 1)
    expression_in_current_cluster = results_1.step_4.list_of_grms{cluster_index};
    
    
    table_names = results_1.step_4.list_of_grms{cluster_index}.Properties.VariableNames(1:3);
    A = 1:size(results_1.step_2.time_points,1);
    A = A';
    A = strtrim(cellstr(num2str(A))');
    A = strcat('T', A);
    table_names = [table_names A];
        
    results_1.step_4.list_of_grms{cluster_index} = cell2table([table2cell(results_1.step_4.list_of_grms{cluster_index}(:,1:3)) num2cell(interpolate(results_1.step_2.original_time_points', cell2mat(table2cell(expression_in_current_cluster(:,4:size(expression_in_current_cluster,2)))), results_1.step_2.time_points'))]);

    results_1.step_4.list_of_grms{cluster_index}.Properties.VariableNames = table_names;
  end
  
  % Interpolating in the clusters of the second condition
  for cluster_index=1:size(results_2.step_4.list_of_grms, 1)
    expression_in_current_cluster = results_2.step_4.list_of_grms{cluster_index};
    
    table_names = results_2.step_4.list_of_grms{cluster_index}.Properties.VariableNames(1:3);
    A = 1:size(results_2.step_2.time_points,1);
    A = A';
    A = strtrim(cellstr(num2str(A))');
    A = strcat('T', A);
    table_names = [table_names A];
    
    results_2.step_4.list_of_grms{cluster_index} = cell2table([table2cell(results_2.step_4.list_of_grms{cluster_index}(:,1:3)) num2cell(interpolate(results_2.step_2.original_time_points', cell2mat(table2cell(expression_in_current_cluster(:,4:size(expression_in_current_cluster,2)))), results_2.step_2.time_points'))]);
    
    results_2.step_4.list_of_grms{cluster_index}.Properties.VariableNames = table_names;
  end
  
  % Interpolating in the expression level matrix of the first condition
  results_1.step_2.gene_expression = interpolate(results_1.step_2.original_time_points', cell2mat(table2cell(results_1.step_2.gene_expression(:,3:size(results_1.step_2.gene_expression,2)))), results_1.step_2.time_points');

  % Interpolating in the expression level matrix of the second condition
  results_2.step_2.gene_expression = interpolate(results_2.step_2.original_time_points', cell2mat(table2cell(results_2.step_2.gene_expression(:,3:size(results_2.step_2.gene_expression,2)))), results_2.step_2.time_points');
  
  % Interpolating in the std expression level of the modules in the first condition
  expression_as_matrix = cell2mat(table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(:,5:size(results_1.step_3.standardized_gene_expression_sorted_by_F_value,2))));
  interpolated_expression_as_matrix = interpolate(results_1.step_2.original_time_points', expression_as_matrix, results_1.step_2.time_points');
  table_names = results_1.step_3.standardized_gene_expression_sorted_by_F_value.Properties.VariableNames(1:4);
  A = 1:size(results_1.step_2.time_points,1);
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T', A);
  table_names = [table_names A];
  results_1.step_3.standardized_gene_expression_sorted_by_F_value = cell2table([table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(:,1:4)) num2cell(interpolated_expression_as_matrix)]);
  results_1.step_3.standardized_gene_expression_sorted_by_F_value.Properties.VariableNames = table_names;
  
  % Interpolating in the std expression level of the modules in the second condition
  expression_as_matrix = cell2mat(table2cell(results_2.step_3.standardized_gene_expression_sorted_by_F_value(:,5:size(results_2.step_3.standardized_gene_expression_sorted_by_F_value,2))));
  interpolated_expression_as_matrix = interpolate(results_2.step_2.original_time_points', expression_as_matrix, results_2.step_2.time_points');
  table_names = results_2.step_3.standardized_gene_expression_sorted_by_F_value.Properties.VariableNames(1:4);
  A = 1:size(results_2.step_2.time_points,1);
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T', A);
  table_names = [table_names A];  
  results_2.step_3.standardized_gene_expression_sorted_by_F_value = cell2table([table2cell(results_2.step_3.standardized_gene_expression_sorted_by_F_value(:,1:4)) num2cell(interpolated_expression_as_matrix)]);
  results_2.step_3.standardized_gene_expression_sorted_by_F_value.Properties.VariableNames = table_names;
  
  % Interpolating the standardized gene expression of subject 1.  
  expression_as_matrix = cell2mat(table2cell(results_1.step_3.standardized_gene_expression(:,3:size(results_1.step_3.standardized_gene_expression,2))));
  interpolated_expression_as_matrix = interpolate(results_1.step_2.original_time_points', expression_as_matrix, results_1.step_2.time_points');
  table_names = results_1.step_3.standardized_gene_expression.Properties.VariableNames(1:2);  
  A = 1:size(results_1.step_2.time_points,1);
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T', A);
  table_names = [table_names A];
  results_1.step_3.standardized_gene_expression = cell2table([table2cell(results_1.step_3.standardized_gene_expression(:,1:2)) num2cell(interpolated_expression_as_matrix)]);
  results_1.step_3.standardized_gene_expression.Properties.VariableNames = table_names;
  
  % Interpolating the standardized gene expression of subject 2.
  expression_as_matrix = cell2mat(table2cell(results_2.step_3.standardized_gene_expression(:,3:size(results_2.step_3.standardized_gene_expression,2))));
  interpolated_expression_as_matrix = interpolate(results_2.step_2.original_time_points', expression_as_matrix, results_2.step_2.time_points');
  table_names = results_2.step_3.standardized_gene_expression.Properties.VariableNames(1:2);  
  A = 1:size(results_2.step_2.time_points,1);
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T', A);
  table_names = [table_names A];
  results_2.step_3.standardized_gene_expression = cell2table([table2cell(results_2.step_3.standardized_gene_expression(:,1:2)) num2cell(interpolated_expression_as_matrix)]);
  results_2.step_3.standardized_gene_expression.Properties.VariableNames = table_names;

  equalized_results_1 = results_1;
  equalized_results_2 = results_2;
  
end

% This function calculates the splines (x,y) for each row of y_matrix and returns the interpolated
% values (yq_matrix) for each element in xq.
% x is a horizontal vector.
% y_matrix is a matrix. Number of columns is the same as the number of elements in x.
% xq is a horizontal vector.
function yq_matrix = interpolate(x, y_matrix, xq)

%    x = [1 2 3 4 5 6];
%    y_matrix = [[1 4 9 16 25 36]; [1 8 27 64 125 216]];
%    xq = [1.5 2.5];

  yq_matrix = [];

  for i=1:size(y_matrix, 1)
    y = y_matrix(i,:);
    
    % The following condition checks that none of the values in y is NaN.
    % If any of the values of y is Nan then the interpolation cannot be performed and a vector of
    % only Nan values will be used instead of the interpolated values.
    % This check is necessary because the expression data in some GEO series inexplicably come with
    % NaN values.
    if ~isempty(find(isnan(y)))
      yq_matrix = [yq_matrix; repmat(nan, 1, length(xq))];
    else
      yq_matrix = [yq_matrix; spline(x, y, xq)];
    end
  end
end

function output_gene_correlations(results_1, results_2)

  name_of_first_subject = results_1.condition;
  name_of_second_subject = results_2.condition;

  matrix_of_correlations_sorted_by_F_value = [];
  
  for drg_index=1:size(results_1.step_3.standardized_gene_expression_sorted_by_F_value, 1)
    probe_id_of_current_drg_in_first_subject = table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(drg_index, 2));
    
    [indices, not_found] = find_strings_in_cell_array(table2cell(results_2.step_3.standardized_gene_expression_sorted_by_F_value(:,2)), probe_id_of_current_drg_in_first_subject);
    
    if isempty(indices)
      matrix_of_correlations_sorted_by_F_value = [matrix_of_correlations_sorted_by_F_value; [table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(drg_index,1:4)) {nan}]];
    else
      expression_of_current_drg_in_first_subject = cell2mat(table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(drg_index,5:size(results_1.step_3.standardized_gene_expression_sorted_by_F_value,2))));
      
      expression_of_current_drg_in_second_subject = cell2mat(table2cell(results_2.step_3.standardized_gene_expression_sorted_by_F_value(indices(1),5:size(results_2.step_3.standardized_gene_expression_sorted_by_F_value,2))));
      
      correlation_index = corr(expression_of_current_drg_in_first_subject', expression_of_current_drg_in_second_subject', 'type', 'Spearman');
      
      matrix_of_correlations_sorted_by_F_value = [matrix_of_correlations_sorted_by_F_value; [table2cell(results_1.step_3.standardized_gene_expression_sorted_by_F_value(drg_index,1:4)) {correlation_index}]];
    end
  end
  
  [useless_variable, indices_of_sorted_correlations] = sort(cell2mat(matrix_of_correlations_sorted_by_F_value(1:size(results_1.step_3.list_of_top_DRGs,1),5)), 'descend');
  
  matrix_of_drgs_correlations_sorted_by_correlations = matrix_of_correlations_sorted_by_F_value(indices_of_sorted_correlations,:);
  
  [useless_variable, indices_of_sorted_correlations] = sort(cell2mat(matrix_of_correlations_sorted_by_F_value(:,5)), 'descend');
  
  matrix_of_correlations_sorted_by_correlations = matrix_of_correlations_sorted_by_F_value(indices_of_sorted_correlations,:);
  
  correlations_matrix_header = [{'Row in GSE matrix'} {'Probe ID'} {'Gene ID'} {'F value'} {['Spearman correlation between expression in ' name_of_first_subject ' and ' name_of_second_subject]}];
  
  matrix_of_correlations_sorted_by_F_value = [correlations_matrix_header; matrix_of_correlations_sorted_by_F_value];
  
  matrix_of_correlations_sorted_by_correlations = [correlations_matrix_header; matrix_of_correlations_sorted_by_correlations];
  
  matrix_of_drgs_correlations_sorted_by_correlations = [correlations_matrix_header; matrix_of_drgs_correlations_sorted_by_correlations];
  
  writetable(cell2table(matrix_of_drgs_correlations_sorted_by_correlations), ['Expression_level_correlation_of_DRGs_(sorted_by_correlation).csv'], 'WriteVariableNames', false);
  writetable(cell2table(matrix_of_correlations_sorted_by_F_value(1:size(results_1.step_3.list_of_top_DRGs,1),:)), ['Expression_level_correlation_of_DRGs(sorted_by_F_value).csv'], 'WriteVariableNames', false);
  writetable(cell2table(matrix_of_correlations_sorted_by_F_value), ['Expression_level_correlation_of_all_genes(sorted_by_F_value).csv'], 'WriteVariableNames', false);
  writetable(cell2table(matrix_of_correlations_sorted_by_correlations), ['Expression_level_correlation__of_all_genes(sorted_by_correlation).csv'], 'WriteVariableNames', false);
  
end

function output_comparison_plots(name_of_first_subject, list_of_grms, time_points, name_of_second_subject, gene_expression_2, indices_of_DRGs, list_of_genes, list_of_probe_ids)

      global Dynamics4GenomicBigData_HOME;
  
      output_folder = 'Gene_expression_of_both_subjects';
      mkdir(output_folder);
      cd(output_folder);
      
      alpha_threshold = 0.05;
      
      p_value_output_matrix = {'GRM number', 'p-value (Wilcoxon)', 'p-value (KS)', 'p-value (KW)', 'p-value (correlation with mean)', 'p-value (bootstrap)', 'p-value (permutation)', ['p-value < ' num2str(alpha_threshold)], '', 'Spearman correlation coefficient', 'Coefficient < 0.75?'};

      
      number_of_clusters = size(list_of_grms,1);

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
	  
	  expression_of_first_subject = cell2mat(table2cell(list_of_grms{currentClusterIndex}(:, 4:size(list_of_grms{currentClusterIndex}, 2))));
	  expression_of_second_subject = gene_expression_2(cell2mat(table2cell(list_of_grms{currentClusterIndex}(:,1))), :);
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
	  % This section calculates the Spearman correlation between the mean expression of the first
	  % subject's module and the mean expression of the same genes from the second subject.
	  
	  mean_exp1 = mean(expression_of_first_subject);
	  mean_exp2 = mean(expression_of_second_subject);
	  
	  correlation_between_mean_expressions = corr(mean_exp1', mean_exp2', 'type', 'Spearman');
	  
	  % The method above only works if there are more than one gene in the module, i.e., if the
	  % expression matrices have more than one row.
	  % Therefore, if there is only one gene, then some extra manipulation is required.
	  if(size(expression_of_first_subject,1) < 2)
	  	correlation_between_mean_expressions = corr(expression_of_first_subject', expression_of_second_subject', 'type', 'Spearman');
	  end
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
	  
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  
	  % This section calculates the Spearman correlation between each gene's expression in the first individual and the expression in the second individual.
	  % The information is output to a spreadsheet.
	  
	  genes_in_current_cluster = table2cell(list_of_grms{currentClusterIndex}(:,3));
	  probes_in_current_cluster = table2cell(list_of_grms{currentClusterIndex}(:,2));
	  
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

	  if size(expression_of_first_subject, 1) == 1
	    plot(mean(cell2mat(table2cell(list_of_grms{currentClusterIndex}(:, 4:size(list_of_grms{currentClusterIndex}, 2)))), 1),'o-r','LineWidth',1.5);
	  else
	    plot(mean(cell2mat(table2cell(list_of_grms{currentClusterIndex}(:, 4:size(list_of_grms{currentClusterIndex}, 2))))),'o-r','LineWidth',1.5);
	  end
	  

	  y_axis_limits = [min([min(expression_of_first_subject), min(expression_of_second_subject)])-.05,max([max(expression_of_first_subject), max(expression_of_second_subject)])+.05];
	  ylim(y_axis_limits);

	  v = axis;

	  handle=title(['Expression of genes in M', num2str(currentClusterIndex), ' from ', strrep(name_of_first_subject,'_','\_'), '']);

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

	  handle=title(['Expression of the same genes in ', strrep(name_of_second_subject,'_','\_'), ' (\rho = ' num2str(correlation_between_mean_expressions) ')']);

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

function plot_cluster_matches(name_of_first_subject, list_of_grms, time_points, name_of_second_subject, list_of_grms_2, time_points_2)

    global Dynamics4GenomicBigData_HOME;

    list_of_cluster_means = [];
    for i = 1:length(list_of_grms)    
      if size(list_of_grms{i}(:, 4:size(list_of_grms{i}, 2)), 1) == 1
	cluster_mean = mean(cell2mat(table2cell(list_of_grms{i}(:, 4:size(list_of_grms{i}, 2)))),1);
	list_of_cluster_means = [list_of_cluster_means; {cluster_mean}];
      else
	cluster_mean = mean(cell2mat(table2cell(list_of_grms{i}(:, 4:size(list_of_grms{i}, 2)))));
	list_of_cluster_means = [list_of_cluster_means; {cluster_mean}];
      end
    end
    
    list_of_cluster_means_2 = [];
    for i = 1:length(list_of_grms_2)    
      if size(list_of_grms_2{i}(:, 4:size(list_of_grms_2{i}, 2)), 1) == 1
	cluster_mean = mean(cell2mat(table2cell(list_of_grms_2{i}(:, 4:size(list_of_grms_2{i}, 2)))),1);
	list_of_cluster_means_2 = [list_of_cluster_means_2; {cluster_mean}];
      else
	cluster_mean = mean(cell2mat(table2cell(list_of_grms_2{i}(:, 4:size(list_of_grms_2{i}, 2)))));
	list_of_cluster_means_2 = [list_of_cluster_means_2; {cluster_mean}];
      end
    end
    
    list_of_cluster_means = cell2mat(list_of_cluster_means);
    list_of_cluster_means_2 = cell2mat(list_of_cluster_means_2);
    
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
      
      highest_correlation = max(z(current_cluster,:));
      lowest_correlation = min(z(current_cluster,:));
      
      
      expression_of_first_subjects_current_cluster = cell2mat(table2cell(list_of_grms{current_cluster}(:, 4:size(list_of_grms{current_cluster}, 2))));
      expression_of_second_subjects_1st_cluster = cell2mat(table2cell(list_of_grms_2{module_with_lowest_correlation}(:, 4:size(list_of_grms_2{module_with_lowest_correlation}, 2))));
      expression_of_second_subjects_2nd_cluster = cell2mat(table2cell(list_of_grms_2{module_with_highest_correlation}(:, 4:size(list_of_grms_2{module_with_highest_correlation}, 2))));
      
      expression_of_first_subjects_cluster = expression_of_first_subjects_current_cluster;
      name_of_first_subjects_cluster = ['M' num2str(current_cluster)];
      mean_of_first_subjects_cluster = list_of_cluster_means(current_cluster,:);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_lowest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_1st_cluster;
      mean_of_second_subjects_cluster = list_of_cluster_means_2(module_with_lowest_correlation,:);
      
      y_limits = [min([min(expression_of_first_subjects_current_cluster), min(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])-.05,max([max(expression_of_first_subjects_current_cluster), max(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])+.05];
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, [name_of_second_subject ' (\rho = ' num2str(lowest_correlation) ')'], name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_limits);
      
      print('-dpdf', [name_of_first_subjects_cluster '_(lowest_correlation).pdf']);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_highest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_2nd_cluster;
      mean_of_second_subjects_cluster = list_of_cluster_means_2(module_with_highest_correlation,:);
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, [name_of_second_subject ' (\rho = ' num2str(highest_correlation) ')'], name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_limits);
      
      print('-dpdf', [name_of_first_subjects_cluster '_(highest_correlation).pdf']);
      
      close all;
    end
    cd('..');
end

function plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_2, y_axis_limits)
			mainFig = figure('units', 'centimeters', 'position', [0, 0, 45, 20]);
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

      handle=title(['Expression of ', name_of_first_subjects_cluster, ' from ', strrep(name_of_first_subject,'_','\_'), '']);

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
      
      handle=title(['Expression of ', name_of_second_subjects_cluster, ' from ', strrep(name_of_second_subject,'_','\_'), '']);

      set(handle,'Position',[2.35 v(4)*1.03 0]);

      hold off;
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
