function fidxcluster = step_4(yhat, IND_DRG, Time, cutoff, gexp2, INDF, GID_DRG, outputFiles, outputClusterFig, outputClusterByTypeFig)

  flder=pwd;
  
  outputFolder = 'Step_4';
  mkdir(outputFolder);
  
  %  ---------------  Surfaces of Top Genes from F-test  ------------------------
  
  axisLabelFontSize = 30;

  h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

  ind= 0 ;

  surf(yhat(:,IND_DRG),'FaceColor','interp','EdgeColor','none');

  ylim([1,length(Time)]);

  set(gca,'YTick',1:length(Time),'Yticklabel',Time);
  set(gca,'FontSize',11);

  xlim([1,cutoff]);

  zlim([min(min(yhat(:,IND_DRG))),max(max(yhat(:,IND_DRG)))]);
  hold on;
      
  ylabel('Time', 'FontSize', axisLabelFontSize);
  xlabel('Top ranking genes', 'FontSize', axisLabelFontSize);
  zlabel('Expression', 'FontSize', axisLabelFontSize);

  title(['Dynamic Response Genes'], 'FontSize', axisLabelFontSize);
  hold off;

%    print(h,'-dpsc2', '-append', 'Paper_04.ps');
%    movefile('Paper_04.ps', outputFolder);
  
  %  -----------------------------------------------------------------------

  %                       Cluster (IHC)

  %  -----------------------------------------------------------------------



  %Theshold
  alpha = 0.75;

  std_data     = zscore(gexp2(INDF(1:cutoff),:)')';
  [fidxcluster,rmclusters,c,mean_clusters_mat,clusters] = IHC(std_data,alpha);
      
      
  % The following four lines sort the clusters by size.
  [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster, 1), 'descend');
  fidxcluster = fidxcluster(cluster_indexes_by_size);
  clusters = clusters(cluster_indexes_by_size);
  mean_clusters_mat = mean_clusters_mat(cluster_indexes_by_size,:);
  % The previous four lines sort the clusters by size.
      
  n_clusters   = cellfun(@(x) size(x,1),clusters,'UniformOutput', false);

  for l = 1:length(fidxcluster)
    Cluster_IDX(fidxcluster{l}) = l;
  end

  if(outputFiles)
    global Dynamics4GenomicBigData_HOME;
    create_exel_file('Cluster_IDX.xls',Cluster_IDX',1,[],Dynamics4GenomicBigData_HOME);    
    movefile('Cluster_IDX.xls', outputFolder);
  end

  for k=1:1

    sz{k}       = cell2mat(n_clusters);

    ind         = find(sz{k}>99);

    ind1        = find(sz{k}>9 & sz{k}<100);

    ind2        = find(sz{k}>1 & sz{k}<10);

    ind3        = find(sz{k}==1);

    lrg_id{k}   = GID_DRG(vertcat(fidxcluster{ind}));

    med_id{k}   = GID_DRG(vertcat(fidxcluster{ind1}));

    smal_id{k}  = GID_DRG(vertcat(fidxcluster{ind2}));

    sin_id{k}   = GID_DRG(vertcat(fidxcluster{ind3}));

    lrg_ts{k}   = mean_clusters_mat(ind,:);

    med_ts{k}   = mean_clusters_mat(ind1,:);

    smal_ts{k}  = mean_clusters_mat(ind2,:);

    sin_ts{k}   = mean_clusters_mat(ind3,:);

    sizes{k}    = [size(sz{k},1),length(ind),length(ind1),length(ind2),length(ind3)];

  end
  
  if(outputClusterFig)

      [s,ind]=sort(cell2mat(n_clusters),'descend');
      
      number_of_clusters = size(mean_clusters_mat,1);
      number_of_subplots = number_of_clusters;
      
      % Ideally the following two variables should be settable to any values. But for now this works only with 30 and 6.
      number_of_subplots_per_page = 30;
      number_of_columns_per_page = 6;    
      number_of_rows_per_page = number_of_subplots_per_page / number_of_columns_per_page;
      
      number_of_pages = ceil(number_of_subplots / number_of_subplots_per_page);
      number_of_plots_in_last_page = number_of_subplots_per_page;
      if mod(number_of_subplots, number_of_subplots_per_page)~=0
	number_of_plots_in_last_page = mod(number_of_subplots, number_of_subplots_per_page);
      end
      
      cluster_number = 1;

      for b = 1:number_of_pages
      
	  number_of_plots_in_current_page = number_of_subplots_per_page;
	  if(b == number_of_pages) %i.e., if this is the last page
	    number_of_plots_in_current_page = number_of_plots_in_last_page;
	  end

	  h8=figure('units', 'centimeters', 'position', [0, 0, 85, 50]);
	  axisLabelFontSize = 9;
	  
	  set(gcf, 'PaperPositionMode', 'manual');
	  set(gcf, 'PaperUnits', 'centimeters');
	  set(gcf, 'PaperPosition', [0 0 75 50]);
	  set(gcf, 'PaperUnits', 'centimeters');
	  set(gcf, 'PaperSize', [85 50]);

	  for gen = 1:number_of_plots_in_current_page

	      subplot(number_of_rows_per_page,number_of_columns_per_page,gen);

	      plot(clusters{cluster_number}','-*b');
	      
	      set(gca,'XTick', 1:size(Time));
	      set(gca,'XTickLabel', Time);

	      xlabel('Time', 'FontSize', axisLabelFontSize);

	      ylabel('Expression', 'FontSize', axisLabelFontSize);

	      hold on;

	      plot(mean_clusters_mat(cluster_number,:),'o-r','LineWidth',1.5);

	      xlim([0,size(mean_clusters_mat(cluster_number,:),2)]);

	      ylim([min(min(clusters{cluster_number}))-.05,max(max(clusters{cluster_number}))+.05]);

	      v = axis;
	      
	      number_of_genes_in_current_cluster  = s(cluster_number);
	      
	      handle=title(['M' num2str(cluster_number) ' (' num2str(number_of_genes_in_current_cluster) ' genes)' ]);
	      
	      if(number_of_genes_in_current_cluster == 1)
		handle=title(['M' num2str(cluster_number) ' (' num2str(number_of_genes_in_current_cluster) ' gene)' ]);
	      end

	      set(handle,'Position',[2.5 v(4)*1. 0]);

	      hold off;
	      
	      cluster_number = cluster_number + 1;

	  end

	  print(h8,'-dpsc2', '-append', 'Cluster.ps');
      end
    movefile('Cluster.ps', outputFolder);
  end

  if(outputClusterByTypeFig)
    GRMFigure=figure('units', 'centimeters', 'position', [0, 0, 50, 40]);

    axisLabelFontSize = 30;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 -2 50 40]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [65 40]);
    set(gca,'FontSize',11);

    Figure1 = subplot(2,2,1);
    set(gca,'FontSize',11);

    if(~isempty(lrg_ts{1}));

    ribbon(lrg_ts{1}');

    ylim([1,size(lrg_ts{1},2)]);

    if size(lrg_ts{1},1) > 1
      xlim([1,size(lrg_ts{1},1)]);
    end

    zlim([min(min(lrg_ts{1})),max(max(lrg_ts{1}))]);

    ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

    xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

    title('LSM');

    end

    Figure2 = subplot(2,2,2);
    set(gca,'FontSize',11);

    ribbon(med_ts{1}');

    ylim([1,size(med_ts{1},2)]);

    if size(med_ts{1},1) > 1
      xlim([1,size(med_ts{1},1)]);
    end
    zlim([min(min(med_ts{1})),max(max(med_ts{1}))]);

    ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

    xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

    title('MSM');

    Figure3 = subplot(2,2,3);
    set(gca,'FontSize',11);

    ribbon(smal_ts{1}');

    ylim([1,size(smal_ts{1},2)]);

    xlim([1,size(smal_ts{1},1)]);

    zlim([min(min(smal_ts{1})),max(max(smal_ts{1}))]);

    ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

    xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

    title('SSM');

    Figure4 = subplot(2,2,4);
    set(gca,'FontSize',11);

    ribbon(sin_ts{1}');

    ylim([1,size(sin_ts{1},2)]);

    xlim([1,size(sin_ts{1},1)]);

    zlim([min(min(sin_ts{1})),max(max(sin_ts{1}))]);

    ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

    xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

    title('SGM');

    print(GRMFigure,'-dpsc2', '-append', 'GRMs.ps');
    
    movefile('GRMs.ps', outputFolder);

    % *Table 2* provides the number of clusters for each subject and the number of clusters in each category.

    col_hed = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};

    row_hed = strcat(repmat({'Subject '},1,1),cellstr(arrayfun(@num2str, 1:1, 'UniformOutput', false))');

    tmp = round2(vertcat(sizes{:}));

    makeHtmlTable(tmp,[],row_hed,col_hed);
  end
  
  matrix_of_files_descs = [{'File name'} {'Description'}];
  
  matrix_of_files_descs = [matrix_of_files_descs; [{'Cluster_IDX.xls'} {'Cluster indices.'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{'Cluster.ps'} {'Cluster plots.'}]];
  matrix_of_files_descs = [matrix_of_files_descs; [{'GRMs.ps'} {'Clusters plotted by size.'}]];
  
  create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

  movefile('List_and_description_of_output.xls', outputFolder);
  
end