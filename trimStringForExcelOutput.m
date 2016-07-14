function trimmedString = trimStringForExcelOutput(theString)

  % This is the maximum string length that can be stored in a spreadsheet cell.
  % Any attempt to store a longer string will result in an IllegalArgumentException from Java.
  maxLength = 32767;
  
  [r, c] = size(theString);
  
  if(c <= maxLength)
    trimmedString = theString;
  else
    trimmedString = theString(1:maxLength);
  end
end

%  [success,theMessage] = xlswrite(fileName,cellfun(@trimStringForExcelOutput, xlsData, 'UniformOutput', 0),sheetName,startRange);
%  [success,theMessage] = xlswrite(fileName,cellfun(@trimStringForExcelOutput, cellstr(num2str(xlsData)), 'UniformOutput', 0),sheetName,startRange);