function trimmedString = TrimString(stringToTrim)
  
  charactersToRemove = {'_', ' ', '-', '/'};
  
  trimmedString = stringToTrim;
  
  for i=1:length(charactersToRemove)
    trimmedString=strrep(trimmedString, charactersToRemove{i}, '');
  end  
end

%  function timePoint = ExtractTimePoint(stringToSearch)
%    delimiters = {' ', ',', '_'};
%    splitString = strsplit(stringToSearch, delimiters);
%        
%    % The following lines attempt to find the time point measured in hours. The time units may appears as 'hr' or 'hours', for instance.
%    indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hr'));
%    if isempty(indexWhereHourIndicatorIsLocated)
%      indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hour'));
%    end
%    if isempty(indexWhereHourIndicatorIsLocated)
%      indexWhereHourIndicatorIsLocated = find(strcmp([splitString], 'hours'));
%    end
%  
%    % The following line attempts to find the time point measured in minutes.
%    indexWhereMinIndicatorIsLocated  = find(strcmp([splitString], 'min'));
%  
%    if not(isempty(indexWhereHourIndicatorIsLocated))
%      timePoint = splitString{indexWhereHourIndicatorIsLocated-1};
%    elseif not(isempty(indexWhereMinIndicatorIsLocated))
%      % The time point found is in minutes and needs to be changed to hours.
%      timePoint = num2str(str2num(splitString{indexWhereMinIndicatorIsLocated-1}) / 60);
%    end
%  end