% File must be in .CSV format and each line must be as follows.

% sample,time point

% The name of the file must be formatted as <GSENUMBER>_<CONDITIONNAME>_<NUMBEROFTOPDRGS>.csv.

function [GEO_number, condition, number_of_top_DRGs] = read_input_from_csv(filename)
  
  parts = strsplit(filename, '.');
  
  parts = strsplit(parts{1}, '_-_');
  
  GEO_number = parts{1};
  
  condition = parts{2};
  
  number_of_top_DRGs = str2num(parts{3});  

end
