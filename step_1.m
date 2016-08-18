function [Data_GEO,gid,titles,Info,PInfo,geoStruct, gene_expression_of_subject_across_time_points, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = step_1(GEO_number)
  
  Preprocessing_technique = 'Default';
  
  [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);
  
  [gene_expression_of_subject_across_time_points, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct);
  
end