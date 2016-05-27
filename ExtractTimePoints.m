% This function expects a column vector of strings.

% The function finds time points as long as they appear as numeric values followed by a delimiter and then by either 'hr' or 'min'. The numeric value must also be preceded by a delimiter.

% The list of delimiters used is given by the variable of the same name.

function timePoints = ExtractTimePoints(Matrix)

  timePoints = [];
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  delimiters = {' ', ',', '_'};

  format bank;
  
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      splitString = strsplit(stringToSearch, delimiters);
      indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hr'));
      indexWhereMinIndicatorIsLocated  = find(strcmp([splitString], 'min'));
      if not(isempty(indexWhereHourIndicatorIsLocated))
	timePoint = splitString{indexWhereHourIndicatorIsLocated-1};
	timePoints = cat(2, timePoints, {timePoint});
      elseif not(isempty(indexWhereMinIndicatorIsLocated))
	% This time point is in minutes and needs to be changed to hours.
	timePoint = num2str(str2num(splitString{indexWhereMinIndicatorIsLocated-1}) / 60);
	timePoints = cat(2, timePoints, {timePoint});
      end
    end
  end
  timePoints = (str2double(timePoints));  
end
