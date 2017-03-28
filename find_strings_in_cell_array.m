function [indices, not_found] = find_strings_in_cell_array(cell_array_to_search, cell_array_to_search_for)

  indices = [];
  not_found = [];
  not_found_idx = [];
  for i=1:length(cell_array_to_search_for)

    idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
    if(isempty(idx))
      not_found = [not_found; {cell_array_to_search_for{i}}];
      not_found_idx = [not_found_idx; i];
    else
      indices = [indices; idx];
    end
  end
end
