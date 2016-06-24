% This function expects a column vector of strings.

% This function returns a Nx2 matrix, where N is the number of time points found (under normal circumstances N should be the number of rows of the input vector). In each row the first column contains the time value (numeric, e.g., 5) as a cell value and the second column contains the time unit used (string, e.g., 'days').

% The returned time points are converted to the most common time unit in the matrix provided.

% The function finds time points as long as they appear as numeric values followed by a delimiter and then by either 'hr' or 'min', for instance.
% The complete list of accepted labels for time units appears in function ExtractTimePoint further below.

% Example input:
% Matrix={'AAAAA_BBB_CCC_25.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_3435.3min_DDD_XXXX', 'AAAAA_BBB_CCC_65.3hrs_DDD_XXXX', 'AAAAA_BBB_CCC_656.3hrs_DDD_XXXX','AAAAA_BBB_CCC_965.3hrs_DDD_XXXX'}'

% Example output:
%      [ 25.30]    'hours'
%      [ 57.26]    'hours'
%      [ 65.30]    'hours'
%      [656.30]    'hours'
%      [965.30]    'hours'

function timePoints = ExtractTimePoints(Matrix)

  timePoints = [];
  
  [numberOfRows,numberOfColumns] = size(Matrix);

  format bank;
  
  % The following block of nested for loops simply goes through the matrix cells looking for the time points and saving them in variable timePoints. In the ends, all values stored in timePoints are measured in hours.
  for i=1:numberOfRows
    for j=1:numberOfColumns
      stringToSearch = Matrix{i,j};
      [timePoint, timeUnit] = ExtractTimePoint(stringToSearch);
      timePointValueAndUnit=cell(1,2);
      timePointValueAndUnit{1} = timePoint;
      timePointValueAndUnit{2} = timeUnit;
      timePoints = [timePoints; timePointValueAndUnit];
    end
  end
  
  if not(isempty(timePoints)) & sum(sum(cellfun(@isempty, timePoints))) == 0
    mostUsedTimeUnit = findMostUsedTimeUnit(timePoints);
    timePoints = changeAllTimePointsToTimeUnit(timePoints, mostUsedTimeUnit);
  end
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

function [timePoint, timeUnit] = ExtractTimePoint(stringToSearch)

  timePoint = [];
  timeUnit = [];
  
  % The time point may be in weeks and the units may appear written as 'weeks' or 'week', for instance. The following section searches the input string to determine whether the time appears in weeks and the label used.
  weekLabelUsed = '';
  possibleWeekLabels = {'week', 'weeks', 'wk'};
  for i=1:length(possibleWeekLabels)
    indexOfWeekLabelInString=strfind(stringToSearch, possibleWeekLabels{i});
    if not(isempty(indexOfWeekLabelInString))
      weekLabelUsed = possibleWeekLabels{i};
      break;
    end
  end

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
  
  if not(isempty(weekLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, weekLabelUsed);
    timeUnit = 'weeks';
  elseif not(isempty(dayLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, dayLabelUsed);
    timeUnit = 'days';
  elseif not(isempty(hourLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, hourLabelUsed);
    timeUnit = 'hours';
  elseif not(isempty(minuteLabelUsed))
    timePoint = ExtractNumericValueBeforeOrAfterTimeLabel(stringToSearch, minuteLabelUsed);
    timeUnit = 'minutes';
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
  if(indexWhereNumericValueShouldEnd<1)
    indexWhereNumericValueShouldEnd = 1;
  end
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


function mostUsedTimeUnit = findMostUsedTimeUnit(timePoints)

  timeLabelsUsed=unique(timePoints(:,2),'stable');
  numberOfTimesEachTimeLabelIsUsed=cellfun(@(x) sum(ismember(timePoints(:,2),x)),timeLabelsUsed,'un',0);
  
  numberOfTimesEachTimeLabelIsUsedAsMatrix= cat(1, numberOfTimesEachTimeLabelIsUsed{:});

  [num] = max(numberOfTimesEachTimeLabelIsUsedAsMatrix(:));
  [x y] = ind2sub(size(numberOfTimesEachTimeLabelIsUsedAsMatrix),find(numberOfTimesEachTimeLabelIsUsedAsMatrix==num));
  
  mostUsedTimeUnit = timeLabelsUsed{x};
end


function updatedTimePoints = changeAllTimePointsToTimeUnit(timePoints, newTimeUnit)

  keySet =   {'hours-to-minutes', 'hours-to-days', 'hours-to-weeks',
	      'minutes-to-hours', 'minutes-to-days', 'minutes-to-weeks',
	      'days-to-minutes',  'days-to-hours', 'days-to-weeks',
	      'weeks-to-days', 'weeks-to-hours', 'weeks-to-minutes'};
  valueSet = [60, 1/24, (1/24)*(1/7),
	      1/60, (1/60)*(1/24), (1/60)*(1/24)*(1/7),
	      24*60, 24, (1/7),
	      7, 7*24, 7*24*60];
  mapObj = containers.Map(keySet,valueSet);

  updatedTimePoints = timePoints;
  
  [numberOfRows,numberOfColumns] = size(updatedTimePoints);

  format bank;
  
  for i=1:numberOfRows
    timePoint = updatedTimePoints{i,1};
    oldTimeUnit = updatedTimePoints{i,2};
    if not(strcmp(oldTimeUnit, newTimeUnit))
      keyString = strcat(oldTimeUnit, '-to-', newTimeUnit);
      updatedTimePoints{i,1} = timePoint*mapObj(keyString);
      updatedTimePoints{i,2} = newTimeUnit;
    end
  end
end


function trimmedString = TrimString(stringToTrim)
  
  charactersToRemove = {'_', ' ', '-', '/'};
  
  trimmedString = stringToTrim;
  
  for i=1:length(charactersToRemove)
    trimmedString=strrep(trimmedString, charactersToRemove{i}, '');
  end  
end
