function EAS = step_7(clusters_sorted_by_size, IND_DRG, fdgenens, Time, output)
  
  global Dynamics4GenomicBigData_HOME;
  flder = pwd;
  
  

  % %Obtain Smoothed Estimates of the derivative and trajectory.
  for j = 1:length(clusters_sorted_by_size)
    group = IND_DRG(clusters_sorted_by_size{j});
    meanfd = mean_grouped(fdgenens,group);
    TimeEx = linspace(Time(1),Time(end),100)';
    yhatEx(:,j) = eval_fd(TimeEx, meanfd);
    dyhatEx(:,j) = eval_fd(TimeEx, meanfd,1);
  end

  % Obtain LASSO estimate of the parameters.
  EAS   = [];
  Stats =[];

  for j = 1:size(dyhatEx,2)
    [EAS(:,j),Stats{j}] = lasso(yhatEx,dyhatEx(:,j));
  end       
  
  if(output)
    outputFolder = 'Step_7';
    mkdir(outputFolder);
    
    G = (EAS~=0);
    A0 = EAS(G);

    tmp_ind = find(EAS');
    A_tmp   = EAS';
    create_exel_file('Networks.xls',[tmp_ind,A_tmp(tmp_ind)],1,[],Dynamics4GenomicBigData_HOME);
    create_exel_file('Network_matrix.xls',EAS,1,[],Dynamics4GenomicBigData_HOME);
    
    movefile('Networks.xls', outputFolder);
    movefile('Network_matrix.xls', outputFolder);
    
    A = EAS;

    tmp_ind = find(A');
    A_tmp   = A';
    create_exel_file('Networks_Refined.xls',[tmp_ind,A_tmp(tmp_ind)],1,[],Dynamics4GenomicBigData_HOME);
    create_exel_file('Networks_Refined_matrix.xls',A,1,[],Dynamics4GenomicBigData_HOME);
    
    movefile('Networks_Refined.xls', outputFolder);
    movefile('Networks_Refined_matrix.xls', outputFolder);

    matrix_of_files_descs = [{'File name'} {'Description'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Networks.xls'} {'Parameters of the ODE obtained by the two-stage method.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Network_matrix.xls'} {'Full matrix of parameters of the ODE obtained by the two-stage method.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Networks_Refined.xls'} {'Estimated parameters of the ODE.'}]];
    matrix_of_files_descs = [matrix_of_files_descs; [{'Networks_Refined_matrix.xls'} {'Matrix with the estimated parameters of the ODE.'}]];
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
end