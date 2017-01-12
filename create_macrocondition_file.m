function create_macrocondition_file()

  set_paths_and_imports;
  
  global Dynamics4GenomicBigData_HOME;

  GEO_number = input(['Enter the GEO series number: ']);
  
  macrocondition_name = 'ALL';
  
  macrocondition_filename = [GEO_number '_-_' macrocondition_name '.txt'];
  
  cd(['Output/' GEO_number '/Conditions']);
  s = dir('*');
  
  file_list = {s.name}';
  
  file_list(1:2) = [];
  
  cd(Dynamics4GenomicBigData_HOME);
  
  cd(['Input']);
  
  T = cell2table(file_list);
  
  writetable(T, macrocondition_filename, 'WriteVariableNames',false);
  
  cd(Dynamics4GenomicBigData_HOME);

end
