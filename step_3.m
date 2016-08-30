function [fdgenens, yhat, IND_DRG, GID_DRG, INDF] = step_3(Time, yCR, gexp2, n, gid, number_of_top_DRGs_considered, outputFig3, outputFiles)

  flder = pwd;
  
  outputFolder = 'Step_3';
  mkdir(outputFolder);
  
  %  -----------------------------------------------------------------------

  %                        FDA

  %  -----------------------------------------------------------------------
  fdgenens    = [];
  dfgenens    = [];
  gcvgenens   = [];
  lambdagenes = [];
  yhat        = [];
  dyhat       = cell(1,1);
  STDERR = [];
  SSE = [];

  options = optimset('LargeScale',  'off', 'Display',     'on', 'Diagnostics', 'off', 'GradObj',     'off', 'Hessian',     'off', 'TolFun', 1e-8, 'TolX',        1e-8);
  
  for i = 1:1
    
    %  ---------------  set up the b-spline basis  ------------------------
    knots    = Time';
    norder   = 3;
    nbasis   = length(Time) + norder - 1;
    basisobj = create_bspline_basis([min(Time) max(Time)], nbasis, norder, knots);

    %  -----------  Otain optimal smoothing parameter  -----------------
    B              = eval_basis(Time,basisobj);
    R              = eval_penalty(basisobj,2);
    lambdagenes = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, yCR', R);
    fdParobj       = fdPar(basisobj, 2, lambdagenes);
    [fdgenens, dfgenens, gcvgenens,~,SSE]  = smooth_basis(Time, gexp2', fdParobj);
    yhat        = eval_fd(Time, fdgenens);
    dyhat{i}       = eval_fd(Time, fdgenens,1);
    STDERR      = sqrt(sum(SSE)/(n*(length(Time)-dfgenens)));
  end
  
  F    = [];
  INDF = [];

  for i = 1:1
    F = Ftest(gexp2, Time,  fdgenens, dfgenens);
    [SF, INDF] = sort(F,'descend');
  end
  
  for i = 1:1
    IND_DRG = INDF(1:number_of_top_DRGs_considered);
    GID_DRG = gid(IND_DRG);
    DRG= gexp2(IND_DRG,:)';
  end
  
  axisLabelFontSize = 30;
  if(outputFig3)
    h=figure('units', 'centimeters', 'position', [0, 0, 30, 24]);

    clear title;

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 24]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [30 24]);

    ind= 0 ;

    for sub = 1:1

	surf(yhat','FaceColor','interp','EdgeColor','none');
	xlim([Time(1),length(Time)]);
	
	set(gca,'XTick',1:length(Time),'Xticklabel',Time);
	set(gca,'FontSize',11);

	ylim([1,size(yhat,2)]);
	zlim([min(min(yhat)),max(max(yhat))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);
	ylabel('Genes', 'FontSize', axisLabelFontSize);
	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Smooth gene expression curves'], 'FontSize', 20);

	hold off;

    end

    print('Paper_03.pdf','-dpdf');
    movefile('Paper_03.pdf', outputFolder);
  end

  if(outputFiles)
  
    col_hed = {'Df','GCV','log10(\lambda)','Std Error'};
    row_hed = strcat(repmat({'Subject '},1,1),cellstr(arrayfun(@num2str, 1:1, 'UniformOutput', false))');
    tmp = round2([dfgenens, mean(gcvgenens,2),log10(lambdagenes),sum(STDERR,2)]);
    makeHtmlTable(tmp,[],row_hed,col_hed);
  
    Col = 'A':'X';

    global Dynamics4GenomicBigData_HOME;
    
    for i = 1:1
      xlRange = [Col(i) '1'];
      create_exel_file('F_value.xls',F,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest.xls',INDF,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest_DRG.xls',IND_DRG,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Probe_set_ID_Ftest_DRG.xls',GID_DRG,1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('DRG.xls',DRG',i,[],Dynamics4GenomicBigData_HOME);
      
      movefile('F_value.xls', outputFolder);
      movefile('Index_Ftest.xls', outputFolder);
      movefile('Index_Ftest_DRG.xls', outputFolder);
      movefile('Probe_set_ID_Ftest_DRG.xls', outputFolder);
      movefile('DRG.xls', outputFolder);
    end
    
    for i = 1:1
      create_exel_file('Fitted_curves.xls',yhat',i,[],Dynamics4GenomicBigData_HOME);
      create_exel_file('Derivative_Fitted_Curves.xls',dyhat{i}',i,[],Dynamics4GenomicBigData_HOME);
      
      movefile('Fitted_curves.xls', outputFolder);
      movefile('Derivative_Fitted_Curves.xls', outputFolder);
    end
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