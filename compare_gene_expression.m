function [p_value_wilcoxon, p_value_ks, p_value_kruskalwallis, p_value_correlation_with_mean, p_value_bootstrap, p_value_permutation] = compare_gene_expression(expression_of_first_subject, expression_of_second_subject, time_points)
  
  first_data_set = mean(expression_of_first_subject, 1)';
  
  second_data_set = mean(expression_of_second_subject, 1)';
  
  p_value_wilcoxon = signrank(first_data_set, second_data_set);%requires column vectors
  
%    p_value = ranksum(first_data_set, second_data_set);%requires column vectors

  [h, p_value_ks] = kstest2(first_data_set',second_data_set');%requires row vectors

  [p_value_kruskalwallis,tbl,stats] = kruskalwallis([first_data_set second_data_set], [], 'off');%requires column vectors

  p_value_correlation_with_mean = compare_curves_correlation_with_mean(first_data_set, second_data_set);%requires column vectors
  
  p_value_bootstrap = compare_curves_bootstrap(first_data_set, second_data_set, time_points);%requires column vectors
  
  p_value_permutation = compare_curves_permutation(first_data_set, second_data_set, time_points);%requires column vectors
end


function p_value = compare_curves_correlation_with_mean(first_data_set, second_data_set)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_CorrelationWithMean.R ' filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);

end

function p_value = compare_curves_bootstrap(first_data_set, second_data_set, time_points)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);
  
  time_points_filename = 'time_points.csv';

  csvwrite(time_points_filename, [time_points]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_Bootstrap.R ' filename ' ' time_points_filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);
  
  delete(time_points_filename);

end


function p_value = compare_curves_permutation(first_data_set, second_data_set, time_points)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);
  
  time_points_filename = 'time_points.csv';

  csvwrite(time_points_filename, [time_points]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves_Permutation.R ' filename ' ' time_points_filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);
  
  delete(time_points_filename);

end