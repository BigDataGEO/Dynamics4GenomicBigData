function geo_struct = get_geo_data(GEO_number)

  cache_folder_name = 'GEO_cache';
  path_to_cached_file = [cache_folder_name '/' GEO_number '.soft'];
  
  if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file, 'file'))
    geo_struct = geoseriesread(path_to_cached_file);
  else
    mkdir(cache_folder_name);
    geo_struct = getgeodata(GEO_number, 'ToFile', path_to_cached_file);
  end
end