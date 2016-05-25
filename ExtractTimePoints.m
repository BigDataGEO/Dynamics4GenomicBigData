function timePoints = ExtractTimePoints(Matrix)
  
  
  delimiters = {' ', ',', '_'};
  
  [numberOfRows,numberOfColumns] = size(Matrix);
  
  timePoints    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), Matrix, 'UniformOutput', false));
  
  if(length(timePoints)~=numberOfRows)
  
  timePoints = [];    
      
  for i=1:numberOfRows
      for j=1:numberOfColumns
          stringToSearch = Matrix{i,j};
          splitString = strsplit(stringToSearch, delimiters);
          indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hr'));
          indexWhereMinIndicatorIsLocated  = find(strcmp([splitString], 'min'));
          if not(isempty(indexWhereHourIndicatorIsLocated))
              timePoint = splitString{indexWhereHourIndicatorIsLocated-1};
              timePoints = cat(2, timePoints, {timePoint});
          end
          if not(isempty(indexWhereMinIndicatorIsLocated))
              timePoint = splitString{indexWhereMinIndicatorIsLocated-1};
              timePoints = cat(2, timePoints, {timePoint});
          end
      end
  end
  timePoints = (str2double(unique(timePoints)));
  end
  
  
end
