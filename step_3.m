function [list_of_DRGs, indices_of_DRGs, indices_of_genes_sorted_by_F_value, smooth_gene_expression, fd_smooth_coefficients] = step_3(list_of_genes, gene_expression, time_points, number_of_top_DRGs_considered, smooth_gene_trajectories, output)

  flder = pwd;
  
  total_number_of_genes_in_geo_record = size(list_of_genes,1);
  
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
  
  for i = 1:1
    
    %  ---------------  set up the b-spline basis  ------------------------
    knots    = time_points';
    norder   = 3;
    nbasis   = length(time_points) + norder - 1;
    basisobj = create_bspline_basis([min(time_points) max(time_points)], nbasis, norder, knots);

    %  -----------  Otain optimal smoothing parameter  -----------------
    B              = eval_basis(time_points,basisobj);
    R              = eval_penalty(basisobj,2);
    lambdagenes = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, smooth_gene_trajectories', R);
    fdParobj       = fdPar(basisobj, 2, lambdagenes);
    [fd_smooth_coefficients, degrees_of_freedom, gcvgenens,~,SSE]  = smooth_basis(time_points, gene_expression', fdParobj);
    
    % The F test requires the degrees of freedom as an integer value. And the smooth_basis function does not always return integer values.
    % It is safe to round down to the nearest integer.
    degrees_of_freedom = floor(degrees_of_freedom);
    
    smooth_gene_expression        = eval_fd(time_points, fd_smooth_coefficients);
    derivatives_of_smooth_gene_expression_curves       = eval_fd(time_points, fd_smooth_coefficients,1);
    STDERR      = sqrt(sum(SSE)/(total_number_of_genes_in_geo_record*(length(time_points)-degrees_of_freedom)));
  end

  [F, F_critic] = Ftest(gene_expression, time_points,  fd_smooth_coefficients, degrees_of_freedom);
  [SF, indices_of_genes_sorted_by_F_value] = sort(F,'descend');
  
  % The DRGs will be determined through an upper one-tailed F test.
  % The DRGs will be those whose F statistics are greater than the F_critic (F_{0.05, numerator_df, denominator_df})
  indices_of_DRGs = indices_of_genes_sorted_by_F_value(1:length(find(SF>F_critic)));
  
  % However, if a fixed number of DRGs was specified, then this overrides the above.
  if(number_of_top_DRGs_considered > 0)
    indices_of_DRGs = indices_of_genes_sorted_by_F_value(1:number_of_top_DRGs_considered);
  end
  
  list_of_DRGs = list_of_genes(indices_of_DRGs);
  DRG = gene_expression(indices_of_DRGs,:)';
  
  if(output)
  
    global Dynamics4GenomicBigData_HOME;
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

	surf(smooth_gene_expression','FaceColor','interp','EdgeColor','none');
	xlim([1,length(time_points)]);
	
	set(gca,'XTick',1:length(time_points),'Xticklabel',time_points);
	set(gca,'FontSize',11);

	ylim([1,size(smooth_gene_expression,2)]);
	zlim([min(min(smooth_gene_expression)),max(max(smooth_gene_expression))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);
	ylabel('Genes', 'FontSize', axisLabelFontSize);
	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Smooth gene expression curves'], 'FontSize', 20);

	hold off;

    end

    print('Paper_03.pdf','-dpdf');
    movefile('Paper_03.pdf', outputFolder);
  
    col_hed = {'Df','GCV','log10(\lambda)','Std Error'};
    row_hed = strcat(repmat({'Subject '},1,1),cellstr(arrayfun(@num2str, 1:1, 'UniformOutput', false))');
    tmp = round2([degrees_of_freedom, mean(gcvgenens,2),log10(lambdagenes),sum(STDERR,2)]);
%      makeHtmlTable(tmp,[],row_hed,col_hed);
  
    Col = 'A':'X';

    global Dynamics4GenomicBigData_HOME;
    
    for i = 1:1
      xlRange = [Col(i) '1'];
      create_exel_file('F_value.xls',F,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest.xls',indices_of_genes_sorted_by_F_value,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest_DRG.xls',indices_of_DRGs,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Probe_set_ID_Ftest_DRG.xls',list_of_DRGs,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('DRG.xls',DRG',i,[],Dynamics4GenomicBigData_HOME);
      
      movefile('F_value.xls', outputFolder);
      movefile('Index_Ftest.xls', outputFolder);
      movefile('Index_Ftest_DRG.xls', outputFolder);
      movefile('Probe_set_ID_Ftest_DRG.xls', outputFolder);
      movefile('DRG.xls', outputFolder);
    end
    
    for i = 1:1
      create_exel_file('Fitted_curves.xls',smooth_gene_expression',i,[],Dynamics4GenomicBigData_HOME);
      create_exel_file('Derivative_Fitted_Curves.xls',derivatives_of_smooth_gene_expression_curves',i,[],Dynamics4GenomicBigData_HOME);
      
      movefile('Fitted_curves.xls', outputFolder);
      movefile('Derivative_Fitted_Curves.xls', outputFolder);
    end

    matrix_of_files_descs = [{'File name'} {'Description.'}];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Paper_03.pdf'} {'Smooth expression curves.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'F_value.xls'} {'F statistics.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Index_Ftest.xls'} {'Index of F statistics.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Index_Ftest_DRG.xls'} {'Indices of the DRGs.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Probe_set_ID_Ftest_DRG.xls'} {'Probe set Ids for TRGs.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'DRG.xls'} {'DRG values.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Fitted_curves.xls'} {'Fitted Curves.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Derivative_Fitted_Curves.xls'} {'Derivatives of the Fitted Curves.'}]];
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
  
end