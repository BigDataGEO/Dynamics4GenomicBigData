function [std_data, fidxcluster,rmclusters,c,mean_clusters_mat,clusters, n_clusters, Cluster_IDX] = step_4(N, i, yhat, IND_DRG, Time, cutoff, axisLabelFontSize, gexp2, INDF, path, flder, GID_DRG, outputFiles, outputClusterFig, outputClusterByTypeFig)
  
  %  ---------------  Surfaces of Top Genes from F-test  ------------------------

%    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);
%  
%    ind= 0 ;
%  
%    for sub = 1:N
%  
%      surf(yhat{sub}(:,IND_DRG{i}),'FaceColor','interp','EdgeColor','none');
%  
%      ylim([1,length(Time{sub})]);
%  
%      set(gca,'YTick',1:length(Time{sub}),'Yticklabel',Time{sub});
%      set(gca,'FontSize',11);
%  
%      xlim([1,cutoff]);
%  
%      zlim([min(min(yhat{sub}(:,IND_DRG{i}))),max(max(yhat{sub}(:,IND_DRG{i})))]);
%      hold on;
%        
%      ylabel('Time', 'FontSize', axisLabelFontSize);
%  
%      xlabel('Top ranking genes', 'FontSize', axisLabelFontSize);
%  
%      zlabel('Expression', 'FontSize', axisLabelFontSize);
%  
%      title(['Dynamic Response Genes'], 'FontSize', axisLabelFontSize);
%      hold off;
%  
%    end
%  
%    print('Paper_04.pdf','-dpdf');
  
  %  -----------------------------------------------------------------------

  %                       Cluster (IHC)

  %  -----------------------------------------------------------------------



  %Theshold
  alpha = 0.75;

  for  i = 1:N
      std_data{i}     = zscore(gexp2{i}(INDF{i}(1:cutoff),:)')';
      [fidxcluster{i},rmclusters{i},c{i},mean_clusters_mat{i},clusters{i}] = IHC(std_data{i},alpha);
      n_clusters{i}   = cellfun(@(x) size(x,1),clusters{i},'UniformOutput', false);
  end

  for  i = 1:N
    for l = 1:length(fidxcluster{i})
      Cluster_IDX{i}(fidxcluster{i}{l}) = l;
    end
  end

  if(outputFiles)
    for i = 1:N
      create_exel_file('Cluster_IDX.xls',Cluster_IDX{i}',i,[],path);
    end
    
    disp(strcat('This is a link to the Cluster Indexes <a href="',flder,'/Cluster_IDX.xls">Cluster_Index</a>.'));
  end
  
  for k=1:N

    sz{k}       = cell2mat(n_clusters{k});

    ind         = find(sz{k}>99);

    ind1        = find(sz{k}>9 & sz{k}<100);

    ind2        = find(sz{k}>1 & sz{k}<10);

    ind3        = find(sz{k}==1);

    lrg_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind}));

    med_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind1}));

    smal_id{k}  = GID_DRG{k}(vertcat(fidxcluster{k}{ind2}));

    sin_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind3}));

    lrg_ts{k}   = mean_clusters_mat{k}(ind,:);

    med_ts{k}   = mean_clusters_mat{k}(ind1,:);

    smal_ts{k}  = mean_clusters_mat{k}(ind2,:);

    sin_ts{k}   = mean_clusters_mat{k}(ind3,:);

    sizes{k}    = [size(sz{k},1),length(ind),length(ind1),length(ind2),length(ind3)];

  end
  
  if(outputClusterFig)
    for i=1:N

      [s,ind]=sort(cell2mat(n_clusters{i}),'descend');
      
      number_of_clusters = size(mean_clusters_mat{i},1);
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
  %      for b = 1:1 % For testing purposes, output only the first page.

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

	      subplot(number_of_rows_per_page,number_of_columns_per_page,gen)

	      ind2 = ind(gen+(b-1).*number_of_plots_in_current_page);

	      plot(clusters{i}{ind2}','-*b');
	      
	      set(gca,'XTick', 1:size(Time{i}));
	      set(gca,'XTickLabel', Time{i});

	      xlabel('Time', 'FontSize', axisLabelFontSize);

	      ylabel('Expression', 'FontSize', axisLabelFontSize);

	      hold on;

	      plot(mean_clusters_mat{i}(ind2,:),'o-r','LineWidth',1.5);

	      xlim([0,size(mean_clusters_mat{i}(ind2,:),2)]);

	      ylim([min(min(clusters{i}{ind2}))-.05,max(max(clusters{i}{ind2}))+.05]);

	      v = axis;

	      number_of_genes_in_current_cluster  = s(gen+(b-1).*number_of_plots_in_current_page);
	      
	      handle=title(['M' num2str(cluster_number) ' (' num2str(number_of_genes_in_current_cluster) ' genes)' ]);

	      set(handle,'Position',[2.5 v(4)*1. 0]);

	      hold off;
	      
	      cluster_number = cluster_number + 1;

	  end

	  print(h8,'-dpsc2', '-append', 'Cluster.ps');
      end
    end
    

    disp(strcat('This is a link to the Cluster Plots <a href="',flder,'/Cluster.ps">Cluster_Plots</a>.'));
  
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

    % *Table 2* provides the number of clusters for each subject and the number of clusters in each category.

    col_hed = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};

    row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');

    tmp = round2(vertcat(sizes{:}));

    makeHtmlTable(tmp,[],row_hed,col_hed);
  end