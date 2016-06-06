% This function expects a column vector of strings.

% The function finds time points as long as they appear as numeric values followed by a delimiter and then by either 'hr' or 'min', for instance.
% The complete list of accepted labels for time units appears in function ExtractTimePoint further below.

% Example input:
% Matrix={'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_3435.3min_DDD_XXXX', 'AAAAA_BBB_CCC_65.3hrs_DDD_XXXX'}'

% Example output:
% [25.3, 57.26, 65.3]

function timePoints = ExtractTimePoints(Matrix)

  timePoints = ExtractTimePointsInHours(Matrix);
  
  timePoints = UniformizeTimePoints(timePoints);
end


% This function expects a column vector of strings.

% The function finds time points as long as they appear as numeric values followed by a delimiter and then by either 'hr' or 'min', for instance.
% The complete list of accepted labels for time units appears in function ExtractTimePoint further below.

% Example input:
% Matrix={'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_3435.3min_DDD_XXXX', 'AAAAA_BBB_CCC_65.3hrs_DDD_XXXX'}'

% Example output:
% [25.3, 57.26, 65.3]

function timePointsInHours = ExtractTimePointsInHours(Matrix)

  timePointsInHours = [];
  
  [numberOfRows,numberOfColumns] = size(Matrix);

  format bank;
  
  % The following block of nested for loops simply goes through the matrix cells looking for the time points and saving them in variable timePointsInHours. In the ends, all values stored in timePointsInHours are measured in hours.
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      timePoint = ExtractTimePoint(stringToSearch);
      timePoint = num2str(timePoint);
      timePointsInHours = cat(2, timePointsInHours, {timePoint});
    end
  end
  timePointsInHours = (str2double(timePointsInHours));  
end


% The following are a few examples that illustrate what this function receives as arguments and what it returns.

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3 hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3mins_DDD_XXXX'
% timeLabelUsed = 'mins'

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3 mins_DDD_XXXX'
% timeLabelUsed = 'mins'

function timePoint = ExtractTimePoint(stringToSearch)

  timePoint = [];

  % The time point may be in days and the units may appear written as 'days' or 'day', for instance. The following section searches the input string to determine whether the time appears in days and the label used.
  dayLabelUsed = '';
  possibleDayLabels = {'day', 'days', 'dd'};
  for i=1:length(possibleDayLabels)
    indexOfDayLabelInString=strfind(stringToSearch, possibleDayLabels{i});
    if not(isempty(indexOfDayLabelInString))
      dayLabelUsed = possibleDayLabels{i};
      break;
    end
  end
  
  % The time point may be in hours and the units may appear written as 'hours' or 'hrs', for instance. The following section searches the input string to determine whether the time appears in hours and the label used.
  hourLabelUsed = '';
  possibleHourLabels = {'hour', 'hours', 'hr', 'hrs'};
  for i=1:length(possibleHourLabels)
    indexOfHourLabelInString=strfind(stringToSearch, possibleHourLabels{i});
    if not(isempty(indexOfHourLabelInString))
      hourLabelUsed = possibleHourLabels{i};
      break;
    end
  end

  % The time point may be in minutes and the units may appear written as 'minutes' or 'mins', for instance. The following section searches the input string to determine whether the time appears in minutes and the label used.
  minuteLabelUsed = '';
  possibleMinuteLabels = {'minute', 'minutes', 'min', 'mins'};
  for i=1:length(possibleMinuteLabels)
    indexOfMinuteLabelInString=strfind(stringToSearch, possibleMinuteLabels{i});
    if not(isempty(indexOfMinuteLabelInString))
      minuteLabelUsed = possibleMinuteLabels{i};
      break;
    end
  end
  

  if not(isempty(dayLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, dayLabelUsed);
    timePoint = timePoint * 24;
  elseif not(isempty(hourLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, hourLabelUsed);
  elseif not(isempty(minuteLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, minuteLabelUsed);
    timePoint = timePoint / 60;
  end
end


% The following are a few examples that illustrate what this function receives as arguments and what it returns.

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3 hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

function numericValueBeforeOrAfterTimeLabel = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, timeLabelUsed)
  
  % First possibility: the time value is before the time unit label. For example, '32 min' or '1 hr'.
  numericValueBeforeOrAfterTimeLabel = ExtractNumericValueBeforeTimeLabel(stringToSearch, timeLabelUsed);
  
  % Second possibility: the time value is after the time unit label. For example, 'Day 3' or 'Min 5'.
  if isempty(numericValueBeforeOrAfterTimeLabel)
    numericValueBeforeOrAfterTimeLabel = ExtractNumericValueAfterTimeLabel(stringToSearch, timeLabelUsed);
  end
end


% The following are a few examples that illustrate what this function receives as arguments and what it returns.

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

% With the following arguments, the function returns 25.3.
% stringToSearch = 'AAAAA_BBB_CCC_25.3 hrs_DDD_XXXX'
% timeLabelUsed = 'hrs'

function numericValueBeforeTimeLabel = ExtractNumericValueBeforeTimeLabel(stringToSearch, timeLabelUsed)

  indexWhereNumericValueShouldEnd = strfind(stringToSearch, timeLabelUsed) - 1;
  indexWhereNumericValueShouldStart = indexWhereNumericValueShouldEnd;
  
  while true
    
    numericValueBeforeTimeLabel = str2num(TrimString(stringToSearch(indexWhereNumericValueShouldStart:indexWhereNumericValueShouldEnd)));
    
    indexWhereNumericValueShouldStart = indexWhereNumericValueShouldStart - 1;
    
    if indexWhereNumericValueShouldStart <= 0
      indexWhereNumericValueShouldStart = 1; % The substring must start at least at the first position.
      break;
    end
    
    % This condition checks if the numeric value read so far is indeed numeric.
    % If true, then the value is numeric. False, otherwise.
    if not(isempty(numericValueBeforeTimeLabel))
      
      % This condition checks if the character immediately before the numeric value that has been read so far is not numeric.
      % If true, then this character is not numeric and the numeric value read so far is THE numeric value.
      % If false, then this character is also numeric and must be added to the numeric value being read.
      if isempty(str2num(strcat('0',strrep(stringToSearch(indexWhereNumericValueShouldStart:indexWhereNumericValueShouldStart), ' ', '_'))))
	indexWhereNumericValueShouldStart = indexWhereNumericValueShouldStart + 1;
	break;
      end
    end    
  end
  
  numericValueBeforeTimeLabel = str2num(TrimString(stringToSearch(indexWhereNumericValueShouldStart:indexWhereNumericValueShouldEnd)));
end


function numericValueAfterTimeLabel = ExtractNumericValueAfterTimeLabel(stringToSearch, timeLabelUsed)

  indexWhereNumericValueShouldStart = strfind(stringToSearch, timeLabelUsed) + length(timeLabelUsed);
  indexWhereNumericValueShouldEnd = indexWhereNumericValueShouldStart;
  
  while true
    
    numericValueAfterTimeLabel = str2num(TrimString(stringToSearch(indexWhereNumericValueShouldStart:indexWhereNumericValueShouldEnd)));
    
    indexWhereNumericValueShouldEnd = indexWhereNumericValueShouldEnd + 1;
    
    if indexWhereNumericValueShouldEnd > length(stringToSearch)
      indexWhereNumericValueShouldEnd = length(stringToSearch); % The substring must end at most at the last position.
      break;
    end
    
    % This condition checks if the numeric value read so far is indeed numeric.
    % If true, then the value is numeric. False, otherwise.
    if not(isempty(numericValueAfterTimeLabel))
      
      % This condition checks if the character immediately after the numeric value that has been read so far is not numeric.
      % If true, then this character is not numeric and the numeric value read so far is THE numeric value.
      % If false, then this character is also numeric and must be added to the numeric value being read.
      if isempty(str2num(strcat('0',strrep(stringToSearch(indexWhereNumericValueShouldEnd:indexWhereNumericValueShouldEnd), ' ', '_'))))
	indexWhereNumericValueShouldEnd = indexWhereNumericValueShouldEnd - 1;
	break;
      end
    end    
  end
  
  numericValueAfterTimeLabel = str2num(TrimString(stringToSearch(indexWhereNumericValueShouldStart:indexWhereNumericValueShouldEnd)));
end


function trimmedString = TrimString(stringToTrim)
  
  charactersToRemove = {'_', ' ', '-', '/'};
  
  trimmedString = stringToTrim;
  
  for i=1:length(charactersToRemove)
    trimmedString=strrep(trimmedString, charactersToRemove{i}, '');
  end  
end

function uniformizedTimePoints = UniformizeTimePoints(timePoints)

  uniformizedTimePoints = timePoints;
end
