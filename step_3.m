function [gene_expression_sorted_by_F_value, standardized_gene_expression_sorted_by_F_value, number_of_statistically_significant_DRGs, smooth_gene_expression, fd_smooth_coefficients] = step_3(gene_expression, time_points, smooth_gene_trajectories, number_of_top_DRGs_considered, output)

  global Dynamics4GenomicBigData_HOME;

  flder = pwd;
  
  gene_expression_mat = cell2mat(table2cell(gene_expression(:,3:size(gene_expression,2))));
  standardized_gene_expression_mat = zscore(gene_expression_mat')';
  
  list_of_probe_ids = strtrim(table2cell(gene_expression(:,1)));
  list_of_genes = strtrim(table2cell(gene_expression(:,2)));
  
  total_number_of_genes_in_geo_record = size(gene_expression_mat,1);
  
  %  -----------------------------------------------------------------------

  %                        FDA

  %  -----------------------------------------------------------------------
  fd_smooth_coefficients    = [];
  degrees_of_freedom    = [];
  gcvgenens   = [];
  lambdagenes = [];
  smooth_gene_expression        = [];
  derivatives_of_smooth_gene_expression_curves       = [];
  STDERR = [];
  SSE = [];

  options = optimset('LargeScale',  'off', 'Display',     'on', 'Diagnostics', 'off', 'GradObj',     'off', 'Hessian',     'off', 'TolFun', 1e-8, 'TolX',        1e-8);
  
  %  ---------------  set up the b-spline basis  ------------------------
  knots    = time_points';
  norder   = 3;
  nbasis   = length(time_points) + norder - 1;
  basisobj = create_bspline_basis([min(time_points) max(time_points)], nbasis, norder, knots);

  %  -----------  Otain optimal smoothing parameter  -----------------
  B = eval_basis(time_points,basisobj);
  R = eval_penalty(basisobj,2);
  lambdagenes = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, smooth_gene_trajectories', R);
  fdParobj = fdPar(basisobj, 2, lambdagenes);
  [fd_smooth_coefficients, degrees_of_freedom, gcvgenens,~,SSE]  = smooth_basis(time_points, gene_expression_mat', fdParobj);
    
  % The F test requires the degrees of freedom as an integer value. And the smooth_basis function does not always return integer values.
  % It is safe to round down to the nearest integer.
  degrees_of_freedom = floor(degrees_of_freedom);
    
  smooth_gene_expression = eval_fd(time_points, fd_smooth_coefficients);
  smooth_gene_expression = smooth_gene_expression';  
  smooth_gene_expression = cell2table([list_of_probe_ids list_of_genes num2cell(smooth_gene_expression)], 'VariableNames', gene_expression.Properties.VariableNames);
  
  
  derivatives_of_smooth_gene_expression_curves = eval_fd(time_points, fd_smooth_coefficients,1);
  derivatives_of_smooth_gene_expression_curves = derivatives_of_smooth_gene_expression_curves';
  STDERR = sqrt(sum(SSE)/(total_number_of_genes_in_geo_record*(length(time_points)-degrees_of_freedom)));

  [F, F_critic] = Ftest(gene_expression_mat, time_points,  fd_smooth_coefficients, degrees_of_freedom);
  
  [SF, indices_of_genes_sorted_by_F_value] = sort(F,'descend');
  indices_of_top_DRGs_in_series_matrix = indices_of_genes_sorted_by_F_value(1:number_of_top_DRGs_considered);
  list_of_top_DRGs = list_of_genes(indices_of_top_DRGs_in_series_matrix);
  list_of_probe_ids_sorted_by_F_value = strtrim(list_of_probe_ids(indices_of_genes_sorted_by_F_value));
  list_of_genes_sorted_by_F_value = strtrim(list_of_genes(indices_of_genes_sorted_by_F_value));
  gene_expression_sorted_by_F_value = gene_expression_mat(indices_of_genes_sorted_by_F_value,:);
  standardized_gene_expression_sorted_by_F_value = standardized_gene_expression_mat(indices_of_genes_sorted_by_F_value,:);
  smooth_gene_expression_sorted_by_F_value = cell2mat(table2cell(smooth_gene_expression(indices_of_genes_sorted_by_F_value,3:size(smooth_gene_expression,2))));

  A = 1:size(time_points,1);
  A = A';
  A = strtrim(cellstr(num2str(A))');
  A = strcat('T', A);
  
  gene_expression_sorted_by_F_value = cell2table([num2cell(indices_of_genes_sorted_by_F_value) list_of_probe_ids_sorted_by_F_value list_of_genes_sorted_by_F_value num2cell(SF) num2cell(gene_expression_sorted_by_F_value)], 'VariableNames', [{'Row_index_in_GSE_matrix'} {'Probe_ID'} {'Gene_name'} {'F_score'} A]);
  
  standardized_gene_expression_sorted_by_F_value = cell2table([num2cell(indices_of_genes_sorted_by_F_value) list_of_probe_ids_sorted_by_F_value list_of_genes_sorted_by_F_value num2cell(SF) num2cell(standardized_gene_expression_sorted_by_F_value)], 'VariableNames', [{'Row_index_in_GSE_matrix'} {'Probe_ID'} {'Gene_name'} {'F_score'} A]);
  
  
  % The DRGs will be determined through an upper one-tailed F test.
  % The DRGs will be those whose F statistics are greater than the F_critic (F_{0.05, numerator_df, denominator_df})
  number_of_statistically_significant_DRGs = length(indices_of_genes_sorted_by_F_value(1:length(find(SF>F_critic))));
  
  if(output)
    outputFolder = 'Step_3';
    mkdir(outputFolder);
    
    axisLabelFontSize = 30;
    
    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    clear title;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 24]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [30 24]);

    ind= 0 ;

    for sub = 1:1

	surf(cell2mat(table2cell(smooth_gene_expression(:,3:size(smooth_gene_expression,2)))),'FaceColor','interp','EdgeColor','none');
	xlim([1,length(time_points)]);
	
	set(gca,'XTick',1:length(time_points),'Xticklabel',time_points);
	set(gca,'FontSize',11);

	ylim([1,size(cell2mat(table2cell(smooth_gene_expression(:,3:size(smooth_gene_expression,2))))',2)]);
	zlim([min(min(cell2mat(table2cell(smooth_gene_expression(:,3:size(smooth_gene_expression,2))))')),max(max(cell2mat(table2cell(smooth_gene_expression(:,3:size(smooth_gene_expression,2))))'))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);
	ylabel('All genes', 'FontSize', axisLabelFontSize);
	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Smooth gene expression of all genes'], 'FontSize', 20);

	hold off;

    end

    print('Smooth_expression_curves.pdf','-dpdf');
    movefile('Smooth_expression_curves.pdf', outputFolder);
    
    axisLabelFontSize = 30;

    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    ind= 0 ;

    surf(smooth_gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs, :)','FaceColor','interp','EdgeColor','none');

    ylim([1,length(time_points)]);

    set(gca,'YTick',1:length(time_points),'Yticklabel',time_points);
    set(gca,'FontSize',11);

    xlim([1,number_of_statistically_significant_DRGs]);

    zlim([min(min(smooth_gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs, :)')),max(max(smooth_gene_expression_sorted_by_F_value(1:number_of_statistically_significant_DRGs, :)'))]);
    hold on;
	
    ylabel('Time', 'FontSize', axisLabelFontSize);
    xlabel('Dynamic response genes', 'FontSize', axisLabelFontSize);
    zlabel('Expression', 'FontSize', axisLabelFontSize);

    title(['Expression of the Dynamic Response Genes'], 'FontSize', axisLabelFontSize);
    hold off;
    
    saveas(gcf, 'Smooth_expression_of_DRGs.png');
    movefile('Smooth_expression_of_DRGs.png', outputFolder);
    
    
    axisLabelFontSize = 30;

    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    ind= 0 ;

    surf(smooth_gene_expression_sorted_by_F_value(1:number_of_top_DRGs_considered, :)','FaceColor','interp','EdgeColor','none');

    ylim([1,length(time_points)]);

    set(gca,'YTick',1:length(time_points),'Yticklabel',time_points);
    set(gca,'FontSize',11);

    xlim([1,number_of_top_DRGs_considered]);

    zlim([min(min(smooth_gene_expression_sorted_by_F_value(1:number_of_top_DRGs_considered, :)')),max(max(smooth_gene_expression_sorted_by_F_value(1:number_of_top_DRGs_considered, :)'))]);
    hold on;
	
    ylabel('Time', 'FontSize', axisLabelFontSize);
    xlabel('Top ranking genes', 'FontSize', axisLabelFontSize);
    zlabel('Expression', 'FontSize', axisLabelFontSize);

    title(['Expression of the top ' num2str(number_of_top_DRGs_considered) ' DRGs'], 'FontSize', axisLabelFontSize);
    hold off;
    
    saveas(gcf, 'Smooth_expression_of_top_DRGs.png');
    movefile('Smooth_expression_of_top_DRGs.png', outputFolder);
    
    close all;
    
    Col = 'A':'X';
    xlRange = [Col(1) '1'];
    
    cd(outputFolder);
    
    writetable(gene_expression_sorted_by_F_value, ['gene_expression_sorted_by_F_value.csv'], 'WriteVariableNames', true);
    
    writetable(standardized_gene_expression_sorted_by_F_value, ['standardized_gene_expression_sorted_by_F_value.csv'], 'WriteVariableNames', true);
    
    writetable(smooth_gene_expression, ['smooth_gene_expression.csv'], 'WriteVariableNames', true);
    writetable(cell2table(num2cell(derivatives_of_smooth_gene_expression_curves)), ['derivatives_of_smooth_gene_expression_curves.csv'], 'WriteVariableNames', false);
    writetable(cell2table({number_of_statistically_significant_DRGs}), ['number_of_statistically_significant_DRGs.csv'], 'WriteVariableNames', false);
    
    cd('..');

    matrix_of_files_descs = [{'File name'} {'Description.'}];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Smooth_expression_curves.pdf'} {'Smooth expression curves.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Smooth_expression_of_DRGs.png'} {'Smooth expression of the statistically significant DRGs.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Smooth_expression_of_top_DRGs.png'} {['Smooth expression of the top ' num2str(number_of_top_DRGs_considered) ' genes, sorted by F-score.']}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'gene_expression_sorted_by_F_value.csv'} {'Gene expression of all the probes, sorted by F-value and listing the associated gene names.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'standardized_gene_expression_sorted_by_F_value.csv'} {'Standardized gene expression of all the probes, sorted by F-value and listing the associated gene names.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'number_of_statistically_significant_DRGs.csv'} {'The number of probes that are statistically-significant DRGs.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'smooth_gene_expression.csv'} {['Smooth expression curves of the top ' num2str(number_of_top_DRGs_considered) ' DRGs.']}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'derivatives_of_smooth_gene_expression_curves.csv'} {['Derivatives of the smooth expression curves of the top ' num2str(number_of_top_DRGs_considered) ' DRGs.']}]];
    
    writetable(cell2table(matrix_of_files_descs), 'List_and_description_of_output.csv', 'WriteVariableNames', false);

    movefile('List_and_description_of_output.csv', outputFolder);
  end
  
end


function [gcv_m]=multiple_GCV_fun(lambda,B,y,R)

    ngenes = size(y,2);
    for g=1:ngenes
    [gcv(g)]=GCV_fun(lambda,B,y(:,g),R);
    end
    
    gcv_m = nansum(gcv);    
    
end

function [F, F_critic] = Ftest(gene_expression, time_points, fd_smooth_coefficients, degrees_of_freedom)

  n = length(time_points);

  % --- Null Hypiothesis is that the gene_t = mu_t
  % as gene_expression is already centered
  %[m,r,n] = size(gene_expression);

  RSS_H0 = 0;
  %for i = 1:r
  %d = reshape(gene_expression(:,i,:),m,n);     
  RSS_H0 = RSS_H0 + nansum(gene_expression.^2,2);
  %end

  % --- Alternative Hypiothesis is that the gene_t = f(t)
  yhat    = eval_fd(time_points, fd_smooth_coefficients);
  RSS_H1  = 0;
  %for i = 1:r
  %d = reshape(gene_expression(:,i,:),m,n);     
  RSS_H1 = RSS_H1 + nansum((gene_expression-yhat').^2,2);
  %end

  % --- F - test adjusted for degrees of freedom
  F = ((RSS_H0 - RSS_H1)/(degrees_of_freedom-1))./(RSS_H1/(n-degrees_of_freedom));
  
  if(any(RSS_H1<1e-8))
    ind = find(RSS_H1<1e-8);     
    F(ind) = min(F);
  end
  
  % The DRGs will be determined through an upper one-tailed F test.
  % The DRGs will be those whose F statistics are greater than the F_{alpha, numerator_df, denominator_df}
  
  d1 = degrees_of_freedom-1;
  d2 = n-degrees_of_freedom;
  
  % The area below the F curve in (-infty, F_critic) is 0.95. Thus the are in (F_critic, infty) is 0.05.
  F_critic = finv(0.95, d1, d2);
end
 
  
