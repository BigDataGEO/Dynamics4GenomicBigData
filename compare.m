function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;

  fprintf('\n');
  GEO_number = input(['Please enter the accession number of your dataset enclosed in single quotes (e.g., ''GSE52428''): ']);
  
  try
    geo_struct = get_geo_data(GEO_number);
  catch
    display(['Could not retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.']);
    return;
  end
  
  [list_of_genes, raw_gene_expression, raw_time_points, name_of_first_subject, condition, gene_ID_type, number_of_top_DRGs_considered] = step_1(GEO_number);
    
  display(sprintf('\nYou have successfully entered the data for the first subject.'));
    
  prompt = '\nNow you will be required to enter the same information for the second subject (press enter to continue).';
  
  uselessVariable = input(prompt);

  [list_of_genes_2, raw_gene_expression_2, raw_time_points_2, name_of_second_subject, condition_2, gene_ID_type_2, number_of_top_DRGs_considered_2] = step_1(GEO_number);
    
  output_folder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',condition);
    
  mkdir(output_folder);
  cd(output_folder);
    
  % Steps 2, 3, and 4 of the pipeline are run for the first subject.
  [gene_expression, time_points, smooth_gene_trajectories] = step_2(raw_gene_expression, raw_time_points, false);
  [list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, fd_smooth_coefficients] = step_3(list_of_genes, gene_expression, time_points, number_of_top_DRGs_considered, smooth_gene_trajectories, false);
  [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(gene_expression, time_points, list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, number_of_top_DRGs_considered, false);
    
  % Steps 2, 3, and 4 of the pipeline are run for the second subject.    
  [gene_expression_2, time_points_2, smooth_gene_trajectories_2] = step_2(raw_gene_expression_2, raw_time_points_2, false);
  [list_of_DRGs_2, indices_of_DRGs_2, indices_of_genes_sorted_by_F_value_2, smooth_gene_expression_2, fd_smooth_coefficients_2] = step_3(list_of_genes_2, gene_expression_2, time_points_2, number_of_top_DRGs_considered_2, smooth_gene_trajectories_2, false);
  [list_of_gene_clusters_2, gene_expression_by_cluster_2, list_of_cluster_means_2] = step_4(gene_expression_2, time_points_2, list_of_DRGs_2, indices_of_DRGs_2, indices_of_genes_sorted_by_F_value_2, smooth_gene_expression_2, number_of_top_DRGs_considered_2, false);
  
  output_comparison_plots(name_of_first_subject, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_2, indices_of_DRGs);
    
  plot_cluster_matches(name_of_first_subject, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2);
    
  close all;
  cd(Dynamics4GenomicBigData_HOME);
end

function output_comparison_plots(name_of_first_subject, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_2, indices_of_DRGs)

      global Dynamics4GenomicBigData_HOME;
  
      output_folder = 'Gene_expression_of_both_subjects';
      mkdir(output_folder);
      cd(output_folder);
      
      p_values = [];
      alpha_threshold = 0.05;
%        p_value_output_matrix = {'GRM number', 'p-value', ['p-value < ' num2str(alpha_threshold)]};
      
      p_value_output_matrix = {'GRM number', 'p-value (Wilcoxon)', 'p-value (KS)', 'p-value (KW)', 'p-value (correlation with mean)', 'p-value (bootstrap)', 'p-value (permutation)', ['p-value < ' num2str(alpha_threshold)], '', 'Variance is different in all time points?', '# of time points where variance is equal', '# of time points where variance is different'};
      
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
	
	  % THIS MAY NOT BE CORRECT.
	  expression_of_second_subject = gene_expression_2(indices_of_DRGs(list_of_gene_clusters{currentClusterIndex}),:);
	    
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
	  
%  	  plot(list_of_cluster_means(currentClusterIndex,:),'o-g','LineWidth',1.5);
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
	  
	  [p_value_wilcoxon, p_value_ks, p_value_kruskalwallis, p_value_correlation_with_mean, p_value_bootstrap, p_value_permutation] = compare_gene_expression(expression_of_first_subject, expression_of_second_subject, time_points);
	  
	  p_values = [p_values; [p_value_wilcoxon, p_value_ks, p_value_kruskalwallis, p_value_correlation_with_mean, p_value_bootstrap, p_value_permutation]];
	  
%  	  p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(p_value), num2str(p_value < alpha_threshold)}];

	  significant = p_value_wilcoxon < alpha_threshold | p_value_ks < alpha_threshold | p_value_kruskalwallis < alpha_threshold | p_value_correlation_with_mean < alpha_threshold | p_value_bootstrap < alpha_threshold | p_value_permutation < alpha_threshold;
	  
	  variance_p_values = compare_variance_by_time_point(expression_of_first_subject, expression_of_second_subject);
	  
	  number_of_time_points_where_variance_is_different = nnz(variance_p_values < 0.05);
	  
	  number_of_time_points_where_variance_is_equal = min(size(expression_of_first_subject, 2), size(expression_of_second_subject, 2)) - number_of_time_points_where_variance_is_different;
	  
	  variance_is_different_in_all_time_points = (number_of_time_points_where_variance_is_equal == 0);
  
	  p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(p_value_wilcoxon), num2str(p_value_ks), num2str(p_value_kruskalwallis), num2str(p_value_correlation_with_mean), num2str(p_value_bootstrap), num2str(p_value_permutation), num2str(significant), '', num2str(variance_is_different_in_all_time_points), num2str(number_of_time_points_where_variance_is_equal), num2str(number_of_time_points_where_variance_is_different)}];
	  
	  print(h8,'-dpdf', ['M' num2str(currentClusterIndex) '.pdf']);
	  
	  currentClusterIndex = currentClusterIndex + 1;

	  gen = gen + 1;

	end

	close all;
      end
      
      create_exel_file('p-values_control_vs_stimulus.xls',p_value_output_matrix,1,[],Dynamics4GenomicBigData_HOME);
      
      cd('..');
end

function plot_cluster_matches(name_of_first_subject, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2)

    x = list_of_cluster_means;
    y = list_of_cluster_means_2;
    
    z = corr(x', y');
    
    output_folder = 'GRM_matching';
    mkdir(output_folder);
    cd(output_folder);
    
    for current_cluster = 1:size(z,1)
      module_with_highest_correlation = find(z(current_cluster,:) == max(z(current_cluster,:)));
      module_with_lowest_correlation = find(z(current_cluster,:) == min(z(current_cluster,:)));
      
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
