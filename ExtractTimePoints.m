function timePoints = ExtractTimePoints(Matrix)
  timePoints = [];
  
  delimiters = {' ', ',', '_'};
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      splitString = strsplit(stringToSearch, delimiters);
      indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hr'));
      if not(isempty(indexWhereHourIndicatorIsLocated))
	timePoint = splitString{indexWhereHourIndicatorIsLocated-1};
	timePoints = cat(2, timePoints, {timePoint});
      end
    end
  end
  timePoints = sort(str2double(unique(timePoints)));
end
