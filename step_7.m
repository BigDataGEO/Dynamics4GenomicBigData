function [EAS, Stats, G] = step_7(N, clusters_sorted_by_size, IND_DRG, fdgenens, Time)
  
  global Dynamics4GenomicBigData_HOME;
  flder = pwd;

  % %Obtain Smoothed Estimates of the derivative and trajectory.
  for i = 1:N
    for j = 1:length(clusters_sorted_by_size)
      group               = IND_DRG{i}(clusters_sorted_by_size{j});
      meanfd              = mean_grouped(fdgenens{i},group);
      TimeEx{i}           = linspace(Time{i}(1),Time{i}(end),100)';
      yhatEx{i}(:,j)      = eval_fd(TimeEx{i}, meanfd);
      dyhatEx{i}(:,j)     = eval_fd(TimeEx{i}, meanfd,1);
    end
  end



  % Obtain LASSO estimate of the parameters.
  EAS   = cell(N);
  Stats = cell(N);

  for i = 1:N
    for j = 1:size(dyhatEx{i},2)
      [EAS{i}(:,j),Stats{i}{j}] = lasso(yhatEx{i},dyhatEx{i}(:,j));
    end       
    G{i} = (EAS{i}~=0); 
    A0{i} = EAS{i}(G{i});
  end

  for i = 1:N
    tmp_ind = find(EAS{i}');
    A_tmp   = EAS{i}';
    create_exel_file('Networks.xls',[tmp_ind,A_tmp(tmp_ind)],i,[],Dynamics4GenomicBigData_HOME);
    create_exel_file('Network_matrix.xls',EAS{i},i,[],Dynamics4GenomicBigData_HOME);
    csvwrite('Network_matrix.csv',EAS{i});  
  end
  disp(strcat('This is the parameters of the ODE obtained by two-stage method <a href="',flder,'/Networks.xls">Network</a>.'));
  disp(strcat('This is the full matrix of parameters of the ODE obtained by two-stage method <a href="',flder,'/Networks_matrix.xls">Network</a>.'));
  
  % Obtain Refined estimates of the parameters.

  %  optim_options = optimset('Display', 'iter','Algorithm','levenberg-marquardt','TolFun',1.0000e-08,'TolX', 1.0000e-08);

  % % 

  A = cell(N);
  for i = 1:N 
    A{i} = EAS{i};%lsqnonlin(@rss_sp,A0{i},[],[],optim_options,TimeEx{i},yhatEx{i},G{i});
  end

  for i = 1:N
    tmp_ind = find(A{i}');
    A_tmp   = A{i}';
    create_exel_file('Networks_Refined.xls',[tmp_ind,A_tmp(tmp_ind)],i,[],Dynamics4GenomicBigData_HOME);
    create_exel_file('Networks_Refined_matrix.xls',A{i},i,[],Dynamics4GenomicBigData_HOME);
  end

  disp(strcat('This is the estimated parameters of the ODE <a href="',flder,'/Networks_Refined.xls">Parameters</a>.'));
  disp(strcat('This is the matrix with the estimated parameters of the ODE <a href="',flder,'/Networks_Refined_matrix.xls">Parameters</a>.'));