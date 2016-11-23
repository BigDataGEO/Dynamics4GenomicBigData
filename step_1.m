function [raw_gene_expression, raw_time_points] = step_1(geoStruct, samples, time_points)
  
  Data_GEO = double(geoStruct.Data);
  
  [indices_of_samples_in_matrix, not_found] = find_in_cell_array_of_strings(geoStruct.Header.Samples.geo_accession, samples);
  
  raw_gene_expression = Data_GEO(:,indices_of_samples_in_matrix);
  
  raw_time_points = ExtractTimePoints(time_points');
  
  raw_time_points = cell2mat(raw_time_points(:,1));
  
end
