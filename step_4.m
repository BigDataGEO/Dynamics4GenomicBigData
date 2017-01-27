function [list_of_gene_clusters, gene_expression_by_cluster, list_of_cluster_means] = step_4(list_of_probe_ids, list_of_genes, standardized_gene_expression, time_points, list_of_top_DRGs, indices_of_top_DRGs, smooth_gene_expression, output)

  global Dynamics4GenomicBigData_HOME;
  
  flder=pwd;
 
  
  %  -----------------------------------------------------------------------

  %                       Cluster (IHC)

  %  -----------------------------------------------------------------------



  %Theshold
  alpha = 0.75;

  std_data = zscore(standardized_gene_expression(indices_of_top_DRGs,:)')';

  [list_of_gene_clusters, rmclusters, c, list_of_cluster_means, gene_expression_by_cluster] = IHC(std_data, alpha);
      
      
  % The following four lines sort the clusters by size.
  [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', list_of_gene_clusters, 1), 'descend');
  list_of_gene_clusters = list_of_gene_clusters(cluster_indexes_by_size);
  gene_expression_by_cluster = gene_expression_by_cluster(cluster_indexes_by_size);
  list_of_cluster_means = list_of_cluster_means(cluster_indexes_by_size,:);
  % The previous four lines sort the clusters by size.
  
  daSm = 0;
  maxIndex = 0;
  for kkk=1:length(list_of_gene_clusters)
    daSm = daSm + length(list_of_gene_clusters{kkk});
    currentClaster = list_of_gene_clusters{kkk};
    
    for qqq = 1:length(currentClaster)
      if currentClaster(qqq) > maxIndex
	maxIndex = currentClaster(qqq);
      end
    end
  end
      
  n_clusters   = cellfun(@(x) size(x,1),gene_expression_by_cluster,'UniformOutput', false);

  for l = 1:length(list_of_gene_clusters)
    Cluster_IDX(list_of_gene_clusters{l}) = l;
  end
  
  for k=1:1

    sz{k}       = cell2mat(n_clusters);

    ind         = find(sz{k}>99);

    ind1        = find(sz{k}>9 & sz{k}<100);

    ind2        = find(sz{k}>1 & sz{k}<10);

    ind3        = find(sz{k}==1);

    lrg_id{k}   = list_of_top_DRGs(vertcat(list_of_gene_clusters{ind}));

    med_id{k}   = list_of_top_DRGs(vertcat(list_of_gene_clusters{ind1}));

    smal_id{k}  = list_of_top_DRGs(vertcat(list_of_gene_clusters{ind2}));

    sin_id{k}   = list_of_top_DRGs(vertcat(list_of_gene_clusters{ind3}));

    lrg_ts{k}   = list_of_cluster_means(ind,:);

    med_ts{k}   = list_of_cluster_means(ind1,:);

    smal_ts{k}  = list_of_cluster_means(ind2,:);

    sin_ts{k}   = list_of_cluster_means(ind3,:);

    sizes{k}    = [size(sz{k},1),length(ind),length(ind1),length(ind2),length(ind3)];

  end
  
  if(output)
  
    global Dynamics4GenomicBigData_HOME;
    outputFolder = 'Step_4';
    mkdir(outputFolder);
    
    [s,ind]=sort(cell2mat(n_clusters),'descend');

      number_of_clusters = size(list_of_cluster_means,1);
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
	  set(gcf, 'PaperSize', [75 50]);

	  for gen = 1:number_of_plots_in_current_page

	      subplot(number_of_rows_per_page,number_of_columns_per_page,gen);

	      plot(gene_expression_by_cluster{cluster_number}','-*b');
	      
	      set(gca,'XTick', 1:size(time_points));
	      set(gca,'XTickLabel', time_points);

	      xlabel('Time', 'FontSize', axisLabelFontSize);

	      ylabel('Expression', 'FontSize', axisLabelFontSize);

	      hold on;

	      plot(list_of_cluster_means(cluster_number,:),'o-r','LineWidth',1.5);

	      xlim([0,size(list_of_cluster_means(cluster_number,:),2)]);

	      ylim([min(min(gene_expression_by_cluster{cluster_number}))-.05,max(max(gene_expression_by_cluster{cluster_number}))+.05]);

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

	  print(h8,'-dpdf', ['GRMs_' num2str(b) '.pdf']);
	  movefile(['GRMs_' num2str(b) '.pdf'], outputFolder);
      end
      
    GRMFigure=figure('units', 'centimeters', 'position', [0, 0, 50, 40]);

    axisLabelFontSize = 30;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 -2 50 40]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [60 40]);
    set(gca,'FontSize',11);

    Figure1 = subplot(2,2,1);
    set(gca,'FontSize',11);
    
    if(size(lrg_ts{1},1) > 0)

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
    end

    Figure2 = subplot(2,2,2);
    set(gca,'FontSize',11);
    
    if(size(med_ts{1},1) > 0)

      ribbon(med_ts{1}');

      ylim([1,size(med_ts{1},2)]);

      if size(med_ts{1},1) > 1
	xlim([1,size(med_ts{1},1)]);
      end
      zlim([min(min(med_ts{1})),max(max(med_ts{1}))]);

      ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

      xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

      title('MSM');
    end

    Figure3 = subplot(2,2,3);
    set(gca,'FontSize',11);

    if(size(smal_ts{1},1) > 0)
    
      ribbon(smal_ts{1}');

      ylim([1,size(smal_ts{1},2)]);
      
      if size(smal_ts{1},1) > 1
	xlim([1,size(smal_ts{1},1)]);
      end

%        xlim([1,size(smal_ts{1},1)]);

      zlim([min(min(smal_ts{1})),max(max(smal_ts{1}))]);

      ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

      xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

      title('SSM');
    end

    Figure4 = subplot(2,2,4);
    set(gca,'FontSize',11);
    
    if(size(sin_ts{1},1) > 0)

      ribbon(sin_ts{1}');

      ylim([1,size(sin_ts{1},2)]);

      if(size(sin_ts{1},1) > 1)
	xlim([1,size(sin_ts{1},1)]);
      end

      zlim([min(min(sin_ts{1})),max(max(sin_ts{1}))]);

      ylabel('Time (hours)', 'FontSize', axisLabelFontSize);

      xlabel('ith Cluster Center', 'FontSize', axisLabelFontSize);

      title('SGM');
    end
    
    print(GRMFigure,'-dpdf', 'GRMs.pdf');
    
    movefile('GRMs.pdf', outputFolder);
    
    close all;

    % *Table 2* provides the number of clusters for each subject and the number of clusters in each category.

    col_hed = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};

    row_hed = strcat(repmat({'Subject '},1,1),cellstr(arrayfun(@num2str, 1:1, 'UniformOutput', false))');

    tmp = round2(vertcat(sizes{:}));

    matrix_of_files_descs = [{'File name'} {'Description'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Clusters_X.pdf'} {'Cluster plots.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'GRMs.pdf'} {'Clusters plotted by size.'}]];
    
    writetable(cell2table(matrix_of_files_descs), 'List_and_description_of_output.csv', 'WriteVariableNames', false);

    movefile('List_and_description_of_output.csv', outputFolder);
    
    cluster_iteration_ID = 1;
    probe_ids_in_current_cluster = list_of_probe_ids(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
    names_of_genes_in_current_cluster = list_of_genes(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
    
    mkdir([outputFolder '/GRMs']);
    cd([outputFolder '/GRMs'])
    
    for cluster_iteration_ID=1:length(list_of_gene_clusters)
      probe_ids_in_current_cluster = list_of_probe_ids(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
      names_of_genes_in_current_cluster = list_of_genes(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID}));
      
      writetable(cell2table([[{'Row index in GSE matrix'} {'Probe ID'} {'Gene name'} strcat({'t = '}, strtrim(cellstr(strtrim(num2str(time_points)))))']; [num2cell(indices_of_top_DRGs(list_of_gene_clusters{cluster_iteration_ID})) probe_ids_in_current_cluster names_of_genes_in_current_cluster num2cell(gene_expression_by_cluster{cluster_iteration_ID})]]), ['M' num2str(cluster_iteration_ID) '.csv'], 'WriteVariableNames', false);
    end
    
    for cluster_number=1:length(list_of_cluster_means)
      number_of_genes_in_current_cluster  = s(cluster_number);
      gene_expression_plot(gene_expression_by_cluster{cluster_number}, time_points, ['M' num2str(cluster_number) ' (' num2str(number_of_genes_in_current_cluster) ' genes)' ], 'Time', 'Expression', 'Z');
      print(gcf,'-dpdf', ['M' num2str(cluster_number) '.pdf']);      
      close all;
    end
    
    cd('../..');
    
  
  end
  
end


function [rx] = round2(x)
%% Round to two decimal places
 
rx = round(x.*100)./100;
end



function [fidxcluster, rmclusters, c, mean_clusters_mat, clusters] = IHC(data, alpha, dist, link)

tic;

n = size(data,1);

if nargin < 3, dist='Spearman'; end
if nargin < 4, link='average'; end

idxclustertmp = clusterdata(data,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

for j = 1:max(idxclustertmp)
idxcluster{1}{j}=find(idxclustertmp==j);
end

orig_idxcluster{1} = idxcluster{1};

clusters      = cell(length(idxcluster{1}),1);
rmclusters    = cell(length(idxcluster{1}),1);
for j = 1:length(idxcluster{1})
    clusters{j} = data(orig_idxcluster{1}{j},:);
end

mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false); 
mean_clusters_mat = cell2mat(mean_clusters_mat);

k=2;
con =1;
crit = 0;

while(con)
    
     clear('idxclustertmp');
     
     idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

     for j = 1:max(idxclustertmp)
     idxcluster{k}{j}=find(idxclustertmp==j);
     end
    
     clusters      = cell(length(idxcluster{k}),1);
     for j = 1:length(idxcluster{k})
         orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
         clusters{j} = data(orig_idxcluster{k}{j},:);
     end
     
    %% Prune
    for j = 1:length(clusters)
        PCorR = 1-pdist2(clusters{j},mod_mean(clusters{j}),dist);
        if(any(PCorR<(alpha)))
            indrm = find(PCorR<alpha);
            rmclusters{k}{j}= clusters{j}(indrm,:);
            clusters{j}(indrm,:) = [];
            rmclustersidx{k}{j} = orig_idxcluster{k}{j}(indrm);
            orig_idxcluster{k}{j}(indrm) = [];
        end
    end
    
    if(~isempty(rmclusters{k}))      
          orig_idxcluster{k}  = horzcat(orig_idxcluster{k},num2cell(vertcat(rmclustersidx{k}{:}))');
          clusters = vertcat(clusters,num2cell(vertcat(rmclusters{k}{:}),2));
    end

    mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
    mean_clusters_mat = cell2mat(mean_clusters_mat);
    
    %% Stopping Conditions
    
    %1. The indices of the clusters are identical
    %2. The number of iterations is >100
    %3. The number of jumping genes is less than 5%.
    
    for l = 1:length(orig_idxcluster{k})
    markov_ind(orig_idxcluster{k}{l}) = l;
    end

    for p = 1:length(orig_idxcluster{k-1})
    mcl_ind(orig_idxcluster{k-1}{p}) = p;
    end
    
    mem = zeros(n,1);
    mem_a = zeros(n,1);
    D = zeros(n,1);
    for i= 1:n
    ind_g_i = find(markov_ind(i)==markov_ind);
    ind_g_i_2 = find(mcl_ind(i)==mcl_ind);
    mem(ind_g_i) = 1;
    mem_a(ind_g_i_2) = 1;
    D(i) = pdist2(mem',mem_a','hamming');
    W(i) = sum(mem)+sum(mem_a);
    end
    
    c(k) = sum(1 - (D'*n)./W)./n;
         
    if(k>3)
    
    m = mean(abs(diff(c((k-3):k))));
         
    if((k>100) || m<0.02)
         
         con1 = 1;
         while(con1)
             
             k=k+1;
             
             idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance',dist,'linkage',link);

             for j = 1:max(idxclustertmp)
             idxcluster{k}{j}=find(idxclustertmp==j);
             end
             
             clusters      = cell(length(idxcluster{k}),1);
             for j = 1:length(idxcluster{k})
                 orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
                 clusters{j} = data(orig_idxcluster{k}{j},:);
             end
             mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
             mean_clusters_mat = cell2mat(mean_clusters_mat);
             
             C = 1-pdist(mean_clusters_mat,dist);
             
             con1 = max(max(C))>(alpha);
      
         end
         break;
    end
   
    end
     
     %display(['No of cycles: ',num2str(k)]);

     k=k+1;
    
end

fidxcluster = orig_idxcluster{k};

end


function[y] = mod_mean(x)

if(size(x,1)>1) 
    y = mean(x);
else
    y = x ;
end

end
