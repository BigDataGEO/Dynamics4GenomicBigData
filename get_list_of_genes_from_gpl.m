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
  
  if ~isempty(not_found)
    msgID = 'MATLAB:rmpath:DirNotFound3';
    msg = ['There are ' num2str(size(not_found,1)) ' genes that could not be mapped from GSE matrix to the GPL record.'];
    baseException = MException(msgID,msg);
    
    dlmwrite('UnmappedGenes.csv', not_found, '');
    
    throw(baseException);
  end

  list_of_genes = platform_struct.Data(indices, index_of_gpl_column_with_gene_ids);
  
end

