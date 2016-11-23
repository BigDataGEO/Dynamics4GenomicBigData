function [geoStruct, list_of_genes, gene_ID_type, list_of_probe_ids] = get_geo_data(GEO_number)

  cache_folder_name = 'GEO_cache';
  path_to_cached_file = [cache_folder_name '/' GEO_number '.txt'];
  
  try
    if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file, 'file'))
      geoStruct = geoseriesread(path_to_cached_file);
    else
      mkdir(cache_folder_name);
      geoStruct = getgeodata2(GEO_number, 'ToFile', path_to_cached_file);
    end
    
    list_of_genes = [];
    gene_ID_type = [];
    list_of_probe_ids = [];
    
    path_to_cached_file_2 = [cache_folder_name '/' GEO_number '.mat'];
    
    if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file_2, 'file'))
    
      load(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type', 'list_of_probe_ids');
      
    else
      [list_of_genes, list_of_probe_ids] = get_list_of_gene_ids(geoStruct);
      gene_ID_type = capture_type_of_gene_ID(geoStruct.Header.Series.geo_accession);
      
      save(path_to_cached_file_2, 'list_of_genes', 'gene_ID_type', 'list_of_probe_ids');
      
    end
    
  catch causeException
    msgID = 'MATLAB:rmpath:DirNotFound1';
    msg = ['Unable to retrieve dataset ''' GEO_number ''' from the Gene Expression Omnibus.'];
    baseException = MException(msgID,msg);
    
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
  end
end
