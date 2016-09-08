function p_value = compare_curves(first_data_set, second_data_set, time_points)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);
  
  time_points_filename = 'time_points.csv';

  csvwrite(time_points_filename, [time_points]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves.R ' filename ' ' time_points_filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);
  
  delete(time_points_filename);

end
