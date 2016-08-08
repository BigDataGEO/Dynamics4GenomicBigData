
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

    [s,ind]=sort(cell2mat(n_clusters{i}),'descend');
    number_of_subplots_per_page = 4;    
    number_of_columns_per_page = 2;    
    number_of_rows_per_page = number_of_subplots_per_page / number_of_columns_per_page;
    
    for b = 1:floor(2*size(mean_clusters_mat{i},1)./number_of_subplots_per_page)
%      for b = 1:1 % For testing purposes, output only the first page.

        h8=figure('units', 'centimeters', 'position', [0, 0, 85, 50]);
        axisLabelFontSize = 9;
        
        set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', [0 0 75 50]);
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperSize', [85 50]);

	gen = 1;
	
        while gen <= number_of_subplots_per_page
        
	    % The first individual

            subplot(number_of_rows_per_page,number_of_columns_per_page,gen)

            ind2 = cluster_indexes_by_size(currentClusterIndex);

            plot(clusters{i}{ind2}','-*b')

            xlabel('Time', 'FontSize', axisLabelFontSize)

            ylabel('Expression', 'FontSize', axisLabelFontSize)

            hold on;

            plot(mean_clusters_mat{i}(ind2,:),'o-r','LineWidth',1.5)

            xlim([0,size(mean_clusters_mat{i}(ind2,:),2)])

            ylim([min(min(clusters{i}{ind2}))-.05,max(max(clusters{i}{ind2}))+.05])

            v = axis;

            handle=title(['Gene expression of cluster ',num2str(currentClusterIndex), ' (', subject_name, ')']);

            set(handle,'Position',[2.5 v(4)*1. 0]);

            hold off;
            
            
            % The second individual
            
            gen = gen + 1;
            
            
            % THIS MAY NOT BE CORRECT.
            expression_of_second_individual = gexp2_2{i}(IND_DRG{i}(fidxcluster{i}{ind2}),:);
            
            subplot(number_of_rows_per_page,number_of_columns_per_page,gen)


            plot(expression_of_second_individual','-*b')

            xlabel('Time', 'FontSize', axisLabelFontSize)

            ylabel('Expression', 'FontSize', axisLabelFontSize)

            hold on;

            

            xlim([0,size(mean_clusters_mat{i}(ind2,:),2)])

            ylim([min(min(expression_of_second_individual))-.05,max(max(expression_of_second_individual))+.05])

            v = axis;

            handle=title(['Gene expression of cluster ',num2str(currentClusterIndex), ' (', subject_name_2, ')']);

            set(handle,'Position',[2.5 v(4)*1. 0]);

            hold off;
            
            currentClusterIndex = currentClusterIndex + 1;

            gen = gen + 1;
            
            
        end
        
        print(h8,'-dpsc2', '-append', 'Comparison.ps');

        close all;
    end

end























