function platform_struct = get_geo_platfom_data(platform_id)
  cache_folder_name = 'GEO_cache';
  path_to_cached_file = [cache_folder_name '/' platform_id '.txt'];
  
  if(exist(cache_folder_name, 'dir') && exist(path_to_cached_file, 'file'))
    platform_struct = geosoftread(path_to_cached_file);
  else
    mkdir(cache_folder_name);
    platform_struct = getgeodata2(platform_id, 'ToFile', path_to_cached_file);
  end
end

