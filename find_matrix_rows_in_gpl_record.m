function [indices, not_found] = find_matrix_rows_in_gpl_record(cell_array_to_search, cell_array_to_search_for)

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
