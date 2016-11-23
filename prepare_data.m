function [raw_gene_expression, raw_time_points, subject_name, condition, number_of_top_DRGs_considered, sample_codes] = prepare_data(geoStruct)
  
  [raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo] = get_data_from_geo_struct(geoStruct);
  
  [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs_considered, sample_codes] = capture_data(geoStruct, raw_gene_expression_data_in_all_samples, titles, metadata_for_all_samples, PInfo);
  
  raw_time_points = raw_time_points';
  
end
