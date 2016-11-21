function [list_of_genes, list_of_probe_ids] = get_list_of_gene_ids(geoStruct)

  try
  
    row_identifiers_of_gse_matrix = rownames(geoStruct.Data);
    
    list_of_probe_ids = row_identifiers_of_gse_matrix;
    
    platform_id = geoStruct.Header.Series.platform_id;
    
    platform_struct = get_geo_platfom_data(platform_id);
    
    index_of_gpl_column_with_gene_ids = capture_index_of_gpl_column_with_gene_ids(platform_struct);
    
    fprintf('\n');
    display(['Loading gene IDs/names. This can take some time, please wait...']);
    fprintf('\n');
    
    list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix);
  
  catch causeException
    msgID = 'MATLAB:rmpath:DirNotFound2';
    msg = ['Unable to retrieve genes from GSE record in the Gene Expression Omnibus.'];
    baseException = MException(msgID,msg);
    
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
  end
end
