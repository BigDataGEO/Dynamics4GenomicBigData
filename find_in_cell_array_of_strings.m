% The parameters of this function are cell arrays of strings.
% The function returns the indices in the first array where members of the second array were found.
% Example
% Input
% cell_array_to_search = {'A', 'B', 'ABC', 'C', 'EFGZ', 'D', 'F'}
% cell_array_to_search_for = {'ABC', 'EFGZ', 'HOHO'}
% Output
% indices = [3; 5]
% not_found = ['HOHO']
function [indices, not_found] = find_in_cell_array_of_strings(cell_array_to_search, cell_array_to_search_for)
  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
    if(mod(i,1000)==0)
      display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
    end
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
  display(['Read ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
end
