function [Data_GEO, gid, titles, Info, PInfo, geoStruct, Data, Pos, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs_considered, gene_ID_type] = step_1(GEO_number)
  
  Preprocessing_technique = 'Default';
  
  [Data_GEO,gid,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);
  
  [Data, Pos, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs_considered, gene_ID_type] = capture_data(GEO_number, Data_GEO, gid, titles, Info, PInfo, geoStruct);
  
end