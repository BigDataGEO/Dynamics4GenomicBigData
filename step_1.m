function [list_of_genes, raw_gene_expression, raw_time_points, number_of_top_DRGs_considered, gene_ID_type] = step_1(GEO_number)
  
  Preprocessing_technique = 'Default';
  
  [Data_GEO,list_of_genes,titles,Info,PInfo,geoStruct] = Obtain_data_from_GEO_website_user(GEO_number,Preprocessing_technique);
  
  [raw_gene_expression, raw_time_points, str_ind, pr_ind, tb, Subject_name, number_of_top_DRGs_considered, gene_ID_type] = capture_data(GEO_number, Data_GEO, list_of_genes, titles, Info, PInfo, geoStruct);
  
  raw_time_points = raw_time_points';
  
end