function plot_cluster_matches(name_of_first_subject, gene_expression_by_cluster, list_of_cluster_means, time_points, name_of_second_subject, gene_expression_by_cluster_2, list_of_cluster_means_2, time_points_2)

    global Dynamics4GenomicBigData_HOME;

    x = list_of_cluster_means;
    y = list_of_cluster_means_2;
    
    z = corr(x', y', 'type', 'Spearman');
    
    output_folder = 'GRM_matching';
    mkdir(output_folder);
    cd(output_folder);
    
    save('Lok.mat')
    
    for current_cluster = 1:size(z,1)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [sorted_values sorted_indices] = sort(z(current_cluster,:), 'descend');
      correlations_of_current_cluster = [{['Modules from ' name_of_second_subject]} {['Spearman correlation with ' name_of_first_subject '''s M' num2str(current_cluster)]}];
      correlations_of_current_cluster = [correlations_of_current_cluster; [strcat('M', strread(num2str(sorted_indices),'%s')')' strread(num2str(sorted_values),'%s')]];
      create_exel_file(['M' num2str(current_cluster) '.xls'],correlations_of_current_cluster,1,[],Dynamics4GenomicBigData_HOME);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
