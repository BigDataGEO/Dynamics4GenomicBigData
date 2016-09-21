function [list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered, geoStruct] = step_1(GEO_number)
  
  Preprocessing_technique = 'Default';
  
  geoStruct = get_geo_data(GEO_number);
  
  [Data_GEO,list_of_genes,titles,Info,PInfo] = Obtain_data_from_GEO_website_user(geoStruct);
  
  [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs_considered, gene_ID_type] = capture_data(GEO_number, Data_GEO, list_of_genes, titles, Info, PInfo, geoStruct);
  
  raw_time_points = raw_time_points';
  
end