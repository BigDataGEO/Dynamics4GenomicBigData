function [coefficients, adjacency_matrix_of_gene_regulatory_network] = step_6(list_of_gene_clusters, time_points, indices_of_DRGs, fd_smooth_coefficients, output)
  
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
  
  coefficients = adjacency_matrix_of_gene_regulatory_network';
  
  if(output)
  
    outputFolder = 'Step_6';
    mkdir(outputFolder);
    
    row_labels = [];
    column_labels = [{''}];
    
    for row_index = 1:size(coefficients, 1)
      row_labels = [row_labels; {strcat('Coefficients of the ODE for dM', num2str(row_index), '/dt')}];
    end
    
    for column_index = 1:size(coefficients, 2)
      column_labels = [column_labels; {strcat('Coefficient of M', num2str(column_index), '')}];
    end
    
    coefficientsToFile = [row_labels num2cell(coefficients)];
    coefficientsToFile = [column_labels'; coefficientsToFile];    
    coefficientsToFile = cell2table(coefficientsToFile);
    
    writetable(coefficientsToFile, 'Coefficients.csv', 'WriteVariableNames', false, 'Delimiter', ',');
    
    movefile('Coefficients.csv', outputFolder);

    matrix_of_files_descs = [{'File name'} {'Description'}];
    
    matrix_of_files_descs = [matrix_of_files_descs; [{'Coefficients.csv'} {'Full matrix of coeffients of the ODE obtained by the two-stage method.'}]];
    
    create_exel_file('List_and_description_of_output.xls', matrix_of_files_descs, 1, [], Dynamics4GenomicBigData_HOME);

    movefile('List_and_description_of_output.xls', outputFolder);
  end
end