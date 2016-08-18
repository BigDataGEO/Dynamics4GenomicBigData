function [fdgenens, dfgenens, gcvgenens, lambdagenes, yhat, STDERR, SSE, IND_DRG, GID_DRG, DRG, cutoff, INDF, F, axisLabelFontSize] = step_3(N, Time, yCR, gexp2, n, gid, number_of_top_DRGs_considered, outputFig3, outputFiles)

  flder = pwd;
  %  -----------------------------------------------------------------------

  %                        FDA

  %  -----------------------------------------------------------------------
  fdgenens    = cell(N,1);
  dfgenens    = cell(N,1);
  gcvgenens   = cell(N,1);
  lambdagenes = cell(N,1);
  yhat        = cell(N,1);
  dyhat       = cell(N,1);
  STDERR = cell(N,1);
  SSE = cell(N,1);

  options = optimset('LargeScale',  'off', 'Display',     'on', 'Diagnostics', 'off', 'GradObj',     'off', 'Hessian',     'off', 'TolFun', 1e-8, 'TolX',        1e-8);
  
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
  
  cutoff = size(INDF{i},1);
  if(number_of_top_DRGs_considered>0 && number_of_top_DRGs_considered<=size(INDF{i},1))
    cutoff = number_of_top_DRGs_considered;
  end
  
  for i = 1:N
    IND_DRG{i} = INDF{i}(1:cutoff);
    GID_DRG{i} = gid(IND_DRG{i});
    DRG{i}= gexp2{i}(IND_DRG{i},:)';
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

    for sub = 1:N

	surf(yhat{sub}','FaceColor','interp','EdgeColor','none');
	xlim([Time{sub}(1),length(Time{sub})]);
	
	set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub});
	set(gca,'FontSize',11);

	ylim([1,size(yhat{sub},2)]);
	zlim([min(min(yhat{sub})),max(max(yhat{sub}))]);

	xlabel('Time', 'FontSize', axisLabelFontSize);
	ylabel('Genes', 'FontSize', axisLabelFontSize);
	zlabel('Expression', 'FontSize', axisLabelFontSize);

	hold on;

	title(['Smooth gene expression curves'], 'FontSize', 20);

	hold off;

    end

    print('Paper_03.pdf','-dpdf');
  end

  if(outputFiles)
  
    col_hed = {'Df','GCV','log10(\lambda)','Std Error'};
    row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
    tmp = round2([cell2mat(dfgenens), mean(cell2mat(gcvgenens),2),log10(cell2mat(lambdagenes)),sum(cell2mat(STDERR),2)]);
    makeHtmlTable(tmp,[],row_hed,col_hed);
  
    Col = 'A':'X';

    global Dynamics4GenomicBigData_HOME;
    
    for i = 1:N
      xlRange = [Col(i) '1'];
      create_exel_file('F_value.xls',F{i},1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest.xls',INDF{i},1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Index_Ftest_DRG.xls',IND_DRG{i},1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('Probe_set_ID_Ftest_DRG.xls',GID_DRG{i},1,xlRange,Dynamics4GenomicBigData_HOME);
      create_exel_file('DRG.xls',DRG{i}',i,[],Dynamics4GenomicBigData_HOME);
    end

    disp(strcat('This is a link to the F statistics <a href="',flder,'/F_value.xls">F_value</a>.'));
    disp(strcat('This is a link to the Index F statistics <a href="',flder,'/Index_Ftest.xls">Index_Ftest</a>.'));
    disp(strcat('This is a link to the Index of the DRGs <a href="',flder,'/Index_Ftest_DRG.xls">Index_DRGs</a>.'));
    disp(strcat('This is a link to the Probe set Ids for TRGs <a href="',flder,'/Probe_set_ID_Ftest_DRG.xls">Probe_sets_DRGs</a>.'));
    disp(strcat('This is a link to the DRG values <a href="',flder,'/DRG.xls">Index_DRG</a>.'));
    
    for i = 1:N
      create_exel_file('Fitted_curves.xls',yhat{i}',i,[],Dynamics4GenomicBigData_HOME);
      create_exel_file('Derivative_Fitted_Curves.xls',dyhat{i}',i,[],Dynamics4GenomicBigData_HOME);
    end

    disp(strcat('This is a link to the Fitted Curves <a href="',flder,'/Fitted_curves.xls">Fitted_Curves</a>.'));
    disp(strcat('This is a link to the Derivatives of the Fitted Curves <a href="',flder,'/Derivative_Fitted_Curves.xls">Derivative_Fitted_Curves</a>.'));
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  