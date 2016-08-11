close all;
clc;
clear;

path = strcat(pwd,'/');

%if first time running on computer run line 9 and 10 
%cd([path,'SBEToolbox-1.3.3\'])
%install

%Add Paths
addpath(path)
addpath(genpath([path,'fdaM/']))
addpath(genpath([path,'SBEToolbox-1.3.3/']))

[~,GEO_num] = xlsread('GEO_list.xlsx');

py.importlib.import_module('DAVIDWS');

for i = 1:size(GEO_num,1)

  GEO_number = char(GEO_num(i));
  Preprocessing_technique = 'Default';
  [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);

  [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);
  
  display(sprintf('\nYou have successfully entered the data for the first subject.'));
  
  prompt = '\nNow you will be required to enter the same information for the second subject (press enter to continue).';
  
  tm_ind = input(prompt);
  
  [Data_2, Subject_2, Pos_2, str_ind_2, pr_ind_2, tb_2, Subject_name_2] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);  
  
  [~, ~, con] = LCS(char(tb(pr_ind(1),1)),char(tb(pr_ind(end),1)));
  con = strrep(con,' ','_');
  con = strrep(con,'/','_');
  con = strrep(con,'.','_');
  
  flder = strcat(path,'Results/',GEO_number,'/',con);
  
  mkdir(flder)
  cd(flder)
  
  options = struct('format','html','outputDir',flder,'showCode',true);
  
  % This part is extracted from Paper.m.
  [gexp, gexp2, Time, N, n, subject_name] = get_preprocessed_data(Data, Subject, Pos, str_ind);
  [gexp_2, gexp2_2, Time_2, N_2, n_2, subject_name_2] = get_preprocessed_data(Data_2, Subject_2, Pos_2, str_ind_2);

  yCR = Est_Sub_Sel(Time,gexp2,N);

  for i = 1:N
      %  ---------------  set up the b-spline basis  ------------------------
      knots    = Time{i}';
      norder   = 3;
      nbasis   = length(Time{i}) + norder - 1;
      basisobj = create_bspline_basis([min(Time{i}) max(Time{i})], nbasis, norder, knots);

      %  -----------  Otain optimal smoothing parameter  -----------------
      B              = eval_basis(Time{i},basisobj);
      R              = eval_penalty(basisobj,2);
      lambdagenes{i} = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, yCR{i}', R);
      fdParobj       = fdPar(basisobj, 2, lambdagenes{i});
      [fdgenens{i}, dfgenens{i}, gcvgenens{i},~,SSE{i}]  = smooth_basis(Time{i}, gexp2{i}', fdParobj);
      yhat{i}        = eval_fd(Time{i}, fdgenens{i});
      dyhat{i}       = eval_fd(Time{i}, fdgenens{i},1);
      STDERR{i}      = sqrt(sum(SSE{i})/(n*(length(Time{i})-dfgenens{i})));
  end

  F    = cell(N,1);
  INDF = cell(N,1);


  for i = 1:N
      F{i} = Ftest(gexp2{i}, Time{i},  fdgenens{i}, dfgenens{i});
      [SF, INDF{i}] = sort(F{i},'descend');
  end


  for i = 1:N
    cutoff = 3000;
    IND_DRG{i} = INDF{i}(1:cutoff);
    GID_DRG{i} = gid(IND_DRG{i});
    DRG{i}= gexp2{i}(IND_DRG{i},:)';
  end

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

  [uselessVariable, cluster_indexes_by_size] = sort(cellfun('size', fidxcluster{i}, 1), 'descend');
  currentClusterIndex = 1;

  for i=1:N
  
    p_values = [];
    alpha_threshold = 0.05;
    p_value_output_matrix = {'GRM number', 'p-value', ['p-value < ' num2str(alpha_threshold)]};

    [s,ind]=sort(cell2mat(n_clusters{i}),'descend');
    
    number_of_clusters = size(mean_clusters_mat{i},1);

    number_of_subplots = 2 * number_of_clusters;
    
    number_of_subplots_per_page = 4;
    number_of_columns_per_page = 2;    
    number_of_rows_per_page = number_of_subplots_per_page / number_of_columns_per_page;
    
    number_of_pages = ceil(number_of_subplots / number_of_subplots_per_page);
    number_of_plots_in_last_page = number_of_subplots_per_page;
    if mod(number_of_subplots, number_of_subplots_per_page)~=0
      number_of_plots_in_last_page = mod(number_of_subplots, number_of_subplots_per_page);
    end
      
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

      gen = 1;
	  
      while gen <= number_of_plots_in_current_page
      
	ind2 = cluster_indexes_by_size(currentClusterIndex);
	
	expression_of_first_subject = clusters{i}{ind2};
      
	% THIS MAY NOT BE CORRECT.
	expression_of_second_subject = gexp2_2{i}(IND_DRG{i}(fidxcluster{i}{ind2}),:);
	  
	% Plot the expression of the first individual

	subplot(number_of_rows_per_page,number_of_columns_per_page,gen);

	plot(expression_of_first_subject','-*b');
	
	set(gca,'XTick', 1:size(Time{i}));
	set(gca,'XTickLabel', Time{i});

	xlabel('Time', 'FontSize', axisLabelFontSize);
	ylabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	plot(mean_clusters_mat{i}(ind2,:),'o-r','LineWidth',1.5);

	y_axis_limits = [min([min(expression_of_first_subject), min(expression_of_second_subject)])-.05,max([max(expression_of_first_subject), max(expression_of_second_subject)])+.05];
	ylim(y_axis_limits);

	v = axis;

	handle=title(['Gene expression of M',num2str(currentClusterIndex), ' (', subject_name, ')']);

	set(handle,'Position',[2.5 v(4)*1. 0]);

	hold off;
	      
	      
	% Plot the expression of the second individual
	      
	gen = gen + 1;

	subplot(number_of_rows_per_page,number_of_columns_per_page,gen)

	plot(expression_of_second_subject','-*b');
	
	% The following line is intended to plot the mean of the current cluster but it deletes the
	% plotted expression of the second subject in the figure. Commented out for now.
%  	plot(mean_clusters_mat{i}(ind2,:),'o-g','LineWidth',1.5);

	ylim(y_axis_limits);
	
	set(gca,'XTick', 1:size(Time{i}));
	set(gca,'XTickLabel', Time{i});

	xlabel('Time', 'FontSize', axisLabelFontSize)

	ylabel('Expression', 'FontSize', axisLabelFontSize)

	hold on;

	v = axis;

	handle=title(['Gene expression of M',num2str(currentClusterIndex), ' (', subject_name_2, ')']);

	set(handle,'Position',[2.5 v(4)*1. 0]);

	hold off;
	
	p_value = wilcoxon_signed_rank_test(expression_of_first_subject, expression_of_second_subject);
	
	p_values = [p_values; p_value];
	
	p_value_output_matrix = [p_value_output_matrix; {['M' num2str(currentClusterIndex)], num2str(p_value), num2str(p_value < alpha_threshold)}];
	
	currentClusterIndex = currentClusterIndex + 1;

	gen = gen + 1;

      end
	  
      print(h8,'-dpsc2', '-append', 'Comparison.ps');
      
      create_exel_file('p-values_control_vs_stimulus.xls',p_value_output_matrix,i,[],path);

      close all;
    end
  end
  % End of part extracted from Paper.m.
  
  close all;
  cd(path);
end
