function output_comparison_plots(name_of_first_subject, list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_2, indices_of_DRGs)

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
      
      create_exel_file('p-values_control_vs_stimulus.xls',p_value_output_matrix,1,[],Dynamics4GenomicBigData_HOME);
      
      cd('..');
end
