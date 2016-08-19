function compare()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;

  [~,GEO_num] = xlsread('GEO_list.xlsx');

  for i = 1:size(GEO_num,1)

    GEO_number = char(GEO_num(i));
    Preprocessing_technique = 'Default';
    [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);

    [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs, gene_ID_type] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);
    
    display(sprintf('\nYou have successfully entered the data for the first subject.'));
    
    prompt = '\nNow you will be required to enter the same information for the second subject (press enter to continue).';
    
    tm_ind = input(prompt);
    
    [Data_2, Subject_2, Pos_2, str_ind_2, pr_ind_2, tb_2, Subject_name_2, number_of_top_DRGs_2, gene_ID_type_2] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);  
    
    [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
    con = strrep(con,' ','_');
    con = strrep(con,'/','_');
    con = strrep(con,'.','_');
    
    flder = strcat(Dynamics4GenomicBigData_HOME,'Results/',GEO_number,'/',con,'/Comparison/');
    
    mkdir(flder);
    cd(flder);
    
    options = struct('format','html','outputDir',flder,'showCode',true);
    
    % Steps 2, 3, and 4 of the pipeline are run for the first subject.
    [gexp, gexp2, Time, N, n, subject_name, yCR] = step_2(Data, Subject, Pos, str_ind, false, false);  
    [fdgenens, dfgenens, gcvgenens, lambdagenes, yhat, STDERR, SSE, IND_DRG, GID_DRG, DRG, cutoff, INDF, F, axisLabelFontSize] = step_3(N, Time, yCR, gexp2, n, gid, number_of_top_DRGs, false, false);  
    [std_data, fidxcluster,rmclusters,c,mean_clusters_mat,clusters, n_clusters, Cluster_IDX] = step_4(N, i, yhat, IND_DRG, Time, cutoff, axisLabelFontSize, gexp2, INDF, GID_DRG, false, false, false);
    
    % Steps 2, 3, and 4 of the pipeline are run for the second subject.
    [gexp_2, gexp2_2, Time_2, N_2, n_2, subject_name_2, yCR_2] = step_2(Data_2, Subject_2, Pos_2, str_ind_2, false, false);
    [fdgenens_2, dfgenens_2, gcvgenens_2, lambdagenes_2, yhat_2, STDERR_2, SSE_2, IND_DRG_2, GID_DRG_2, DRG_2, cutoff_2, INDF_2, F_2, axisLabelFontSize_2] = step_3(N_2, Time_2, yCR_2, gexp2_2, n_2, gid, number_of_top_DRGs_2, false, false);
    [std_data_2, fidxcluster_2,rmclusters_2,c_2,mean_clusters_mat_2,clusters_2, n_clusters_2, Cluster_IDX_2] = step_4(N_2, i, yhat_2, IND_DRG_2, Time_2, cutoff_2, axisLabelFontSize_2, gexp2_2, INDF_2, GID_DRG_2, false, false, false);

    % The first subject's clusters are sorted by size, so that the first cluster is the largest one
    % and the last cluster is the smallest one.
    [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster{i}, 1), 'descend');
    
    fidxcluster{i} = fidxcluster{i}(cluster_indexes_by_size);
    clusters{i} = clusters{i}(cluster_indexes_by_size);
    mean_clusters_mat{i} = mean_clusters_mat{i}(cluster_indexes_by_size,:);
    
    % The second subject's clusters are sorted by size, so that the first cluster is the largest one
    % and the last cluster is the smallest one.
    [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster_2{i}, 1), 'descend');
    
    fidxcluster_2{i} = fidxcluster_2{i}(cluster_indexes_by_size);
    clusters_2{i} = clusters_2{i}(cluster_indexes_by_size);
    mean_clusters_mat_2{i} = mean_clusters_mat_2{i}(cluster_indexes_by_size,:);
    
    array_indices_of_first_subjects_DRGs = IND_DRG{i};
    
    first_subjects_clusters = fidxcluster{i};
    
    name_of_first_subject = subject_name;
    name_of_second_subject = subject_name_2;
    
    second_subjects_gene_expression = gexp2_2{i};
    
    first_subjects_expression_by_cluster = clusters{i};
    second_subjects_expression_by_cluster = clusters_2{i};
    
    means_of_first_subjects_clusters = mean_clusters_mat{i};
    means_of_second_subjects_clusters = mean_clusters_mat_2{i};
    
    time_points_of_first_subject = Time{i};
    time_points_of_second_subject = Time_2{i};
    
    output_comparison_plots(name_of_first_subject, first_subjects_clusters, first_subjects_expression_by_cluster, means_of_first_subjects_clusters, time_points_of_first_subject, name_of_second_subject, second_subjects_gene_expression, array_indices_of_first_subjects_DRGs);
    
    plot_cluster_matches(name_of_first_subject, first_subjects_expression_by_cluster, means_of_first_subjects_clusters, time_points_of_first_subject, name_of_second_subject, second_subjects_expression_by_cluster, means_of_second_subjects_clusters, time_points_of_second_subject);
    
    close all;
    cd(Dynamics4GenomicBigData_HOME);
  end
end


function plot_cluster_matches(name_of_first_subject, first_subjects_expression_by_cluster, means_of_first_subjects_clusters, time_points_of_first_subject, name_of_second_subject, second_subjects_expression_by_cluster, means_of_second_subjects_clusters, time_points_of_second_subject)

    x = means_of_first_subjects_clusters;
    y = means_of_second_subjects_clusters;
    
    z = corr(x', y');
    
    output_folder = 'GRM_matching';
    mkdir(output_folder);
    cd(output_folder);
    
    for current_cluster = 1:size(z,1)
      module_with_highest_correlation = find(z(current_cluster,:) == max(z(current_cluster,:)));
      module_with_lowest_correlation = find(z(current_cluster,:) == min(z(current_cluster,:)));
      
      expression_of_first_subjects_current_cluster = first_subjects_expression_by_cluster{current_cluster};      
      expression_of_second_subjects_1st_cluster = second_subjects_expression_by_cluster{module_with_lowest_correlation};      
      expression_of_second_subjects_2nd_cluster = second_subjects_expression_by_cluster{module_with_highest_correlation};
      
      expression_of_first_subjects_cluster = expression_of_first_subjects_current_cluster;
      name_of_first_subjects_cluster = ['M' num2str(current_cluster)];
      mean_of_first_subjects_cluster = means_of_first_subjects_clusters(current_cluster,:);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_lowest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_1st_cluster;
      mean_of_second_subjects_cluster = means_of_second_subjects_clusters(module_with_lowest_correlation,:);
      
      y_limits = [min([min(expression_of_first_subjects_current_cluster), min(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])-.05,max([max(expression_of_first_subjects_current_cluster), max(expression_of_second_subjects_1st_cluster), min(expression_of_second_subjects_2nd_cluster)])+.05];
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points_of_first_subject, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_of_second_subject, y_limits);
      
      print('-dpsc2', '-append', [name_of_first_subjects_cluster '.ps']);
      
      name_of_second_subjects_cluster = ['M' num2str(module_with_highest_correlation)];
      expression_of_second_subjects_cluster = expression_of_second_subjects_2nd_cluster;
      mean_of_second_subjects_cluster = means_of_second_subjects_clusters(module_with_highest_correlation,:);
      
      plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points_of_first_subject, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_of_second_subject, y_limits);
      
      print('-dpsc2', '-append', [name_of_first_subjects_cluster '.ps']);
      
      close all;
    end
    cd('..');
end

function plot_expression_of_two_clusters(name_of_first_subject, name_of_first_subjects_cluster, expression_of_first_subjects_cluster, mean_of_first_subjects_cluster, time_points_of_first_subject, name_of_second_subject, name_of_second_subjects_cluster, expression_of_second_subjects_cluster, mean_of_second_subjects_cluster, time_points_of_second_subject, y_axis_limits)
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
	  
      set(gca,'XTick', 1:size(time_points_of_first_subject));
      set(gca,'XTickLabel', time_points_of_first_subject);

      xlabel('Time', 'FontSize', axisLabelFontSize);
      ylabel('Expression', 'FontSize', axisLabelFontSize);
      
      hold on;

      plot(mean_of_first_subjects_cluster,'o-r','LineWidth',1.5);

      ylim(y_axis_limits);

      v = axis;

      handle=title(['Expression of ', name_of_first_subjects_cluster, ' from ', name_of_first_subject, '']);

      set(handle,'Position',[2.5 v(4)*1. 0]);

      hold off;
      
      % Second subplot
      subplot(1,2,2);
      
      plot(expression_of_second_subjects_cluster','-*b');
	  
      set(gca,'XTick', 1:size(time_points_of_second_subject));
      set(gca,'XTickLabel', time_points_of_second_subject);

      xlabel('Time', 'FontSize', axisLabelFontSize);
      ylabel('Expression', 'FontSize', axisLabelFontSize);
      
      hold on;

      plot(mean_of_second_subjects_cluster,'o-r','LineWidth',1.5);

      ylim(y_axis_limits);

      v = axis;
      
      handle=title(['Expression of ', name_of_second_subjects_cluster, ' from ', name_of_second_subject, '']);

      set(handle,'Position',[2.5 v(4)*1. 0]);

      hold off;
end

function output_comparison_plots(name_of_first_subject, first_subjects_clusters, first_subjects_expression_by_cluster, means_of_first_subjects_clusters, time_points_of_first_subject, name_of_second_subject, second_subjects_gene_expression, array_indices_of_first_subjects_DRGs)

      global Dynamics4GenomicBigData_HOME;
  
      output_folder = 'Gene_expression_of_both_subjects';
      mkdir(output_folder);
      cd(output_folder);
      
      p_values = [];
      alpha_threshold = 0.05;
      p_value_output_matrix = {'GRM number', 'p-value', ['p-value < ' num2str(alpha_threshold)]};
      
      number_of_clusters = size(means_of_first_subjects_clusters,1);

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
	  
	  expression_of_first_subject = first_subjects_expression_by_cluster{currentClusterIndex};
	
	  % THIS MAY NOT BE CORRECT.
	  expression_of_second_subject = second_subjects_gene_expression(array_indices_of_first_subjects_DRGs(first_subjects_clusters{currentClusterIndex}),:);
	    
	  % Plot the expression of the first individual

	  subplot(number_of_rows_per_page,number_of_columns_per_page,gen);

	  plot(expression_of_first_subject','-*b');
	  
	  set(gca,'XTick', 1:size(time_points_of_first_subject));
	  set(gca,'XTickLabel', time_points_of_first_subject);

	  xlabel('Time', 'FontSize', axisLabelFontSize);
	  ylabel('Expression', 'FontSize', axisLabelFontSize);

	  hold on;

	  plot(means_of_first_subjects_clusters(currentClusterIndex,:),'o-r','LineWidth',1.5);

	  y_axis_limits = [min([min(expression_of_first_subject), min(expression_of_second_subject)])-.05,max([max(expression_of_first_subject), max(expression_of_second_subject)])+.05];
	  ylim(y_axis_limits);

	  v = axis;

	  handle=title(['Expression of genes in M',num2str(currentClusterIndex), ' from ', name_of_first_subject, '']);

	  set(handle,'Position',[2.5 v(4)*1. 0]);

	  hold off;
		
		
	  % Plot the expression of the second individual
		
	  gen = gen + 1;

	  subplot(number_of_rows_per_page,number_of_columns_per_page,gen)

	  plot(expression_of_second_subject','-*b');
	  
	  hold on;
	  
%  	  plot(means_of_first_subjects_clusters(currentClusterIndex,:),'o-g','LineWidth',1.5);

	  ylim(y_axis_limits);
	  
	  set(gca,'XTick', 1:size(time_points_of_first_subject));
	  set(gca,'XTickLabel', time_points_of_first_subject);

	  xlabel('Time', 'FontSize', axisLabelFontSize)

	  ylabel('Expression', 'FontSize', axisLabelFontSize)

	  hold on;

	  v = axis;

	  handle=title(['Expression of the same genes in ', name_of_second_subject, '']);

	  set(handle,'Position',[2.5 v(4)*1. 0]);

	  hold off;
	  
	  p_value = wilcoxon_signed_rank_test(expression_of_first_subject, expression_of_second_subject);
	  
	  p_values = [p_values; p_value];
	  
	  p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(p_value), num2str(p_value < alpha_threshold)}];
	  
	  print(h8,'-dpdf', ['M' num2str(currentClusterIndex) '.pdf']);
	  
	  currentClusterIndex = currentClusterIndex + 1;

	  gen = gen + 1;

	end

	close all;
      end
      
      create_exel_file('p-values_control_vs_stimulus.xls',p_value_output_matrix,1,[],Dynamics4GenomicBigData_HOME);
      
      cd('..');
end