function [list_of_genes, raw_gene_expression, raw_time_points, subject_name, condition, gene_ID_type, number_of_top_DRGs_considered, geoStruct] = step_1(GEO_number)
  
  Preprocessing_technique = 'Default';
  
  geoStruct = get_geo_data(GEO_number);
  
  [Data_GEO, row_identifiers_of_gse_matrix, titles, Info, PInfo] = get_data_from_geo_struct(geoStruct);
  
  platform_id = geoStruct.Header.Series.platform_id;
  
  platform_struct = getgeodata(platform_id);
  
  [raw_gene_expression, raw_time_points, condition, subject_name, number_of_top_DRGs_considered, gene_ID_type, index_of_gpl_column_with_gene_ids] = capture_data(GEO_number, Data_GEO, row_identifiers_of_gse_matrix, titles, Info, PInfo, geoStruct, platform_struct);
  
  raw_time_points = raw_time_points';
  
  list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix);
  
end

function list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix)

  row_identifiers_in_platform_record = [];
  if(isnumeric([platform_struct.Data{:, 1}]))
    
    numeric_values=[platform_struct.Data{:, 1}];
    
    for k=1:length(numeric_values)
      row_identifiers_in_platform_record = [row_identifiers_in_platform_record; {num2str(numeric_values(k))}];
    end
  else
    row_identifiers_in_platform_record = platform_struct.Data(:, 1);
  end

  % Under normal circumstances, variable not_found should be EMPTY after the following line.
  [indices, not_found] = find_in_cell_array_of_strings(row_identifiers_in_platform_record, row_identifiers_of_gse_matrix);

  list_of_genes = platform_struct.Data(indices, index_of_gpl_column_with_gene_ids);
end

% The parameters of this function are cell arrays of strings.
% The function returns the indices in the first array where members of the second array were found.
% Example
% cell_array_to_search = {'A', 'B', 'ABC', 'C', 'EFGZ', 'D', 'F'}
% cell_array_to_search_for = {'ABC', 'EFGZ', 'HOHO'}
function [indices, not_found] = find_in_cell_array_of_strings(cell_array_to_search, cell_array_to_search_for)
  indices = [];
  not_found = [];
  for i=1:length(cell_array_to_search_for)
      idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
      if(isempty(idx))
	not_found = [not_found; cell_array_to_search_for{i}];
      else
	indices = [indices; idx];
      end
  end
end