function timePoints = InputTimePointsManually()  
  timePoints = [];
  prompt = 'Please enter the time points as a list of text strings enclosed in curved brackets (e.g {''60 min'', ''2 hr'',''3 hr''}):\n';
  timePoints = ExtractTimePoints(input(prompt).');
end
