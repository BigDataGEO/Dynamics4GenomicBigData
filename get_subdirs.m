function subdirs = get_subdirs(folder_name)
  d = dir(folder_name);
  isub = [d(:).isdir];
  subdirs = {d(isub).name}';
  subdirs(ismember(subdirs,{'.','..'})) = [];
end