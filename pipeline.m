set_paths_and_imports;

cd('Input');
s = dir('*.csv');
file_list = {s.name}';

GEO_number={};
condition = {};
samples = {};
time_points = {};
number_of_top_DRGs = {};

for i=1:length(file_list)
  [GEO_number{i}, condition{i}, samples{i}, time_points{i}, number_of_top_DRGs{i}] = read_input(file_list{i});
end

cd('..');

for i=1:length(file_list)
  run_condition(GEO_number{i}, condition{i}, samples{i}, time_points{i}, number_of_top_DRGs{i});
end

fprintf('\n');
display('The analysis is complete for all the subjects/conditions.');

unique_GEO_numbers = unique(GEO_number);

for i=1:length(unique_GEO_numbers)

  fprintf('\n');
  display(['All results from dataset ' unique_GEO_numbers{i} ' have been output to folder ' Dynamics4GenomicBigData_HOME 'Output/' unique_GEO_numbers{i} '/Conditions/']);

  write_study_report(unique_GEO_numbers{i});

  fprintf('\n');
  display(['A consolidated report on all conditions from dataset ' unique_GEO_numbers{i} ' can be found in ' Dynamics4GenomicBigData_HOME 'Output/' unique_GEO_numbers{i} '/paper.pdf']);

end



s = dir('*.csv');
file_list = {s.name}';

for i=1:length(file_list)
  [GEO_number{i}, condition{i}, samples{i}, time_points{i}, number_of_top_DRGs{i}] = read_input(file_list{i});
  
  
  writetable(cell2table([samples{i} time_points{i}]), file_list{i}, 'WriteVariableNames', false);
  
end
