function adjacency_matrix_of_gene_regulatory_network = step_7(list_of_gene_clusters, time_points, indices_of_DRGs, fd_smooth_coefficients, output)
  
  global Dynamics4GenomicBigData_HOME;
  flder = pwd;
  
  % %Obtain Smoothed Estimates of the derivative and trajectory.
  for j = 1:length(list_of_gene_clusters)
    group = indices_of_DRGs(list_of_gene_clusters{j});
    meanfd = mean_grouped(fd_smooth_coefficients,group);
    TimeEx = linspace(time_points(1),time_points(end),100)';
    yhatEx(:,j) = eval_fd(TimeEx, meanfd);
    dyhatEx(:,j) = eval_fd(TimeEx, meanfd,1);
  end

  % Obtain LASSO estimate of the parameters.
  adjacency_matrix_of_gene_regulatory_network   = [];
  Stats =[];

  for j = 1:size(dyhatEx,2)
    [adjacency_matrix_of_gene_regulatory_network(:,j),Stats{j}] = lasso(yhatEx,dyhatEx(:,j));
  end       
  
  if(output)
    outputFolder = 'Step_7';
    mkdir(outputFolder);
    
    G = (adjacency_matrix_of_gene_regulatory_network~=0);
    A0 = adjacency_matrix_of_gene_regulatory_network(G);

    tmp_ind = find(adjacency_matrix_of_gene_regulatory_network');
    A_tmp   = adjacency_matrix_of_gene_regulatory_network';
    create_exel_file('Networks.xls',[tmp_ind,A_tmp(tmp_ind)],1,[],Dynamics4GenomicBigData_HOME);
    create_exel_file('Network_matrix.xls',adjacency_matrix_of_gene_regulatory_network,1,[],Dynamics4GenomicBigData_HOME);
    
    movefile('Networks.xls', outputFolder);
    movefile('Network_matrix.xls', outputFolder);
    
    A = adjacency_matrix_of_gene_regulatory_network;

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