function [raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered] = step_1(geoStruct)
  
  Preprocessing_technique = 'Default';
  
  [Data_GEO, titles, Info, PInfo] = get_data_from_geo_struct(geoStruct);
  
  [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs_considered, gene_ID_type] = capture_data(geoStruct, Data_GEO, titles, Info, PInfo);
  
  raw_time_points = raw_time_points';
  
end
