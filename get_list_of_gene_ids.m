function [geoStruct, list_of_genes] = get_list_of_gene_ids(GEO_number)

%    try

    geoStruct = get_geo_data(GEO_number);
    
    row_identifiers_of_gse_matrix = rownames(geoStruct.Data);
    
    platform_id = geoStruct.Header.Series.platform_id;
    
    platform_struct = getgeodata(platform_id);
    
    index_of_gpl_column_with_gene_ids = capture_index_of_gpl_column_with_gene_ids(platform_struct);
    
    fprintf('\n');
    display(['Loading gene IDs/names. This can take some time, please wait...']);
    fprintf('\n');
    
    list_of_genes = get_list_of_genes_from_gpl(platform_struct, index_of_gpl_column_with_gene_ids, row_identifiers_of_gse_matrix);
  
%    catch
%      msgID = 'get_geo_data:geo_data_not_retrieved';
%      msg = ['Unable to retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.'];
%      baseException = MException(msgID,msg);
%      
%      throw(baseException);
%    end
end