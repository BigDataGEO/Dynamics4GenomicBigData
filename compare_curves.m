function p_value = compare_curves(first_data_set, second_data_set)

  global Dynamics4GenomicBigData_HOME;

  filename = 'curves.csv';

  csvwrite(filename, [first_data_set second_data_set]);

  command = ['Rscript ' Dynamics4GenomicBigData_HOME 'TestTwoCurves.R ' filename];

  [status, cmdout] = system(command);

  substrings = strsplit(cmdout,' ');

  p_value=str2num(substrings{2});
  
  delete(filename);

end
