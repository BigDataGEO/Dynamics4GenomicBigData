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

  list_of_genes = platform_struct.Data(indices, index_of_gpl_column_with_gene_ids);
end

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
  
      % The following line is certain to work. However it is not efficient. Thus it was replaced for the line that follows, which theoretically should be equivalent and more efficient.
      % The rationale is this: if the i-th element of the second array is present in the first array, then it must not be before the i-th position in the first array.
%        idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      
      idx = [];
      if(length(cell_array_to_search) > length(cell_array_to_search_for))
	idx = find(strcmp([cell_array_to_search(i:length(cell_array_to_search))], strtrim(cell_array_to_search_for{i})));
      else
	idx = find(strcmp([cell_array_to_search(:)], strtrim(cell_array_to_search_for{i})));
      end
      
      if(mod(i,1000)==0)
	display(['Loaded ' num2str(i) ' gene IDs/names out of ' num2str(length(cell_array_to_search_for)) '.']);
      end
      
      if(isempty(idx))
	not_found = [not_found; {cell_array_to_search_for{i}}];
	not_found_idx = [not_found_idx; i];
      else
	if(length(cell_array_to_search) > length(cell_array_to_search_for))
	  indices = [indices; (idx+i-1)];
	else
	  indices = [indices; idx];
	end
      end
  end
end