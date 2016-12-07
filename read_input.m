% File must be in .CSV format and each line must be as follows.

% sample,time point

% The name of the file must be formatted as <GSENUMBER>_<CONDITIONNAME>_<NUMBEROFTOPDRGS>.csv.

function [GEO_number, condition, samples, time_points, number_of_top_DRGs] = read_input(filename)
  
  parts = strsplit(filename, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  GEO_number = parts{1};
  
  condition = parts{2};
  
  number_of_top_DRGs = str2num(parts{3});
  
  T = table2cell(readtable(filename, 'ReadVariableNames', false, 'Delimiter', ','));
  
  time_points = T(:,1);
  
  samples = T(:,2);

end
